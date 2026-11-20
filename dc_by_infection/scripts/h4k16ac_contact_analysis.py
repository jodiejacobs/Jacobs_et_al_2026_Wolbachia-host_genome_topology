#!/usr/bin/env python3
"""
Analyze H4K16ac enrichment at differential chromatin contacts on the X chromosome.
Compares real differential interactions vs permutation-based null model to assess significance.

FAST VERSION: pure numpy intersection replaces pybedtools in permutation loops.
pybedtools is used only once for the observed overlap calculation.

logFC convention (consistent throughout):
  logFC > 0 → JW18 wMel  (down-regulated contacts in infected)
  logFC < 0 → JW18 uninf. (up-regulated contacts in uninfected)

Based on:
- Fei et al. - NDF associates with MSL complex via H4K16ac
- Larschan et al. - H4K16ac marks active genes for dosage compensation
"""

import pandas as pd
import numpy as np
import pybedtools
from scipy import stats
import argparse
import matplotlib.pyplot as plt
from pathlib import Path
import gzip
from multiprocessing import Pool, cpu_count
from functools import partial


# ---------------------------------------------------------------------------
# logFC convention helpers — single source of truth
# logFC > 0 = wMel, logFC < 0 = uninf.
# ---------------------------------------------------------------------------
WMEL_MASK  = lambda df: df['logFC'] > 0
UNINF_MASK = lambda df: df['logFC'] < 0


# ===========================================================================
# DATA LOADING
# ===========================================================================

def load_h4k16ac_peaks(chip_file):
    """
    Load H4K16ac ChIP-seq peaks (BED format).
    Returns X-chromosome peaks only as a DataFrame.
    """
    print(f"Loading H4K16ac peaks from {chip_file}")
    opener = gzip.open(chip_file, 'rt') if chip_file.endswith('.gz') else open(chip_file, 'r')
    peaks = []
    with opener as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom = parts[0].replace('chr', '')
                peaks.append({
                    'chrom': chrom,
                    'start': int(parts[1]),
                    'end':   int(parts[2]),
                    'score': float(parts[3]) if len(parts) > 3 else 0
                })

    peaks_df = pd.DataFrame(peaks)
    print(f"  Total peaks loaded: {len(peaks_df)}")
    print(f"  Chromosomes: {sorted(peaks_df['chrom'].unique())}")

    x_peaks = peaks_df[peaks_df['chrom'] == 'X'].copy()
    if len(x_peaks) == 0:
        print("  WARNING: No peaks found on X chromosome!")
        return None

    print(f"  X-chromosome peaks: {len(x_peaks)}")
    print(f"  Score range: {x_peaks['score'].min():.2f} – {x_peaks['score'].max():.2f}")
    return x_peaks


def load_differential_interactions(interactions_file, fdr_threshold=0.05):
    """Load and filter significant differential interactions."""
    print(f"Loading interactions from {interactions_file}")
    interactions = pd.read_csv(interactions_file)

    if 'chr1' in interactions.columns:
        interactions['chr1'] = interactions['chr1'].astype(str).str.replace('chr', '')
        interactions['chr2'] = interactions['chr2'].astype(str).str.replace('chr', '')

    print(f"  Total interactions: {len(interactions)}")

    sig = interactions[
        (interactions['FDR'] < fdr_threshold) &
        (abs(interactions['logFC']) > 1)
    ].copy()

    print(f"  Significant (FDR<{fdr_threshold}, |logFC|>1): {len(sig)}")
    if len(sig) > 0:
        print(f"    Cis:          {(sig['chr1'] == sig['chr2']).sum()}")
        print(f"    Trans:        {(sig['chr1'] != sig['chr2']).sum()}")
        print(f"    wMel  (>0):   {WMEL_MASK(sig).sum()}")
        print(f"    uninf. (<0):  {UNINF_MASK(sig).sum()}")

    return sig, interactions


# ===========================================================================
# OVERLAP CALCULATION — pybedtools (used once on real data)
# ===========================================================================

def calculate_h4k16ac_overlap_vectorized(interactions_df, h4k16ac_peaks, window_size=1000):
    """
    Calculate H4K16ac peak overlap for each interaction using pybedtools.
    Called once on the real data; permutations use the fast numpy version.
    """
    if h4k16ac_peaks is None or len(h4k16ac_peaks) == 0:
        return None

    print(f"  Calculating observed overlaps for {len(interactions_df)} interactions "
          f"(window={window_size} bp)...")

    anchor1_rows, anchor2_rows = [], []
    for idx, row in interactions_df.iterrows():
        anchor1_rows.append({'chrom': row['chr1'],
                              'start': max(0, row['start1'] - window_size),
                              'end':   row['end1'] + window_size, 'idx': idx})
        anchor2_rows.append({'chrom': row['chr2'],
                              'start': max(0, row['start2'] - window_size),
                              'end':   row['end2'] + window_size, 'idx': idx})

    a1_df = pd.DataFrame(anchor1_rows)
    a2_df = pd.DataFrame(anchor2_rows)

    a1_bt    = pybedtools.BedTool.from_dataframe(a1_df[['chrom', 'start', 'end', 'idx']])
    a2_bt    = pybedtools.BedTool.from_dataframe(a2_df[['chrom', 'start', 'end', 'idx']])
    peaks_bt = pybedtools.BedTool.from_dataframe(h4k16ac_peaks[['chrom', 'start', 'end']])

    a1_counts = {int(iv.fields[3]): int(iv.fields[4])
                 for iv in a1_bt.intersect(peaks_bt, wa=True, c=True)}
    a2_counts = {int(iv.fields[3]): int(iv.fields[4])
                 for iv in a2_bt.intersect(peaks_bt, wa=True, c=True)}

    results = []
    for idx in interactions_df.index:
        n1 = a1_counts.get(idx, 0)
        n2 = a2_counts.get(idx, 0)
        results.append({
            'interaction_idx':      idx,
            'anchor1_overlap':      n1 > 0,
            'anchor2_overlap':      n2 > 0,
            'both_anchors_overlap': (n1 > 0) and (n2 > 0),
            'any_anchor_overlap':   (n1 > 0) or  (n2 > 0),
            'n_peaks_anchor1':      n1,
            'n_peaks_anchor2':      n2,
            'total_peaks':          n1 + n2
        })

    return pd.DataFrame(results)


# ===========================================================================
# FAST NUMPY PEAK COUNTING — used inside every permutation
# ===========================================================================

def build_peak_index(h4k16ac_peaks):
    """
    Sort peak starts/ends once; passed to every permutation worker.
    Returns (peak_starts, peak_ends) as sorted int64 arrays.
    """
    arr      = h4k16ac_peaks[['start', 'end']].values.astype(np.int64)
    sort_idx = arr[:, 0].argsort()
    return arr[sort_idx, 0], arr[sort_idx, 1]


def _count_peaks_in_windows(window_starts, window_ends, peak_starts, peak_ends):
    """
    For each window [ws, we], count peaks where peak_start < we AND peak_end > ws.
    Uses searchsorted for O(log M) candidate selection, then vectorized filter.

    Parameters
    ----------
    window_starts, window_ends : 1D int arrays  (N windows)
    peak_starts,  peak_ends    : 1D int arrays  (M peaks, sorted by peak_starts)

    Returns
    -------
    counts : 1D int array (N,)
    """
    counts = np.zeros(len(window_starts), dtype=np.int64)
    for i in range(len(window_starts)):
        ws, we = window_starts[i], window_ends[i]
        right  = np.searchsorted(peak_starts, we, side='left')
        if right == 0:
            continue
        counts[i] = np.sum(peak_ends[:right] > ws)
    return counts


def _total_peaks_for_interactions(interactions_df, peak_starts, peak_ends, window_size):
    """
    Return per-interaction total peak count (anchor1 + anchor2) using numpy only.
    Only processes anchors on the X chromosome (peaks are X-only).
    """
    df = interactions_df.reset_index(drop=True)
    total = np.zeros(len(df), dtype=np.int64)

    for chrom_col, start_col, end_col in [('chr1', 'start1', 'end1'),
                                           ('chr2', 'start2', 'end2')]:
        on_x = (df[chrom_col] == 'X').values
        if not on_x.any():
            continue
        ws = np.maximum(0, df.loc[on_x, start_col].values - window_size).astype(np.int64)
        we = (df.loc[on_x, end_col].values + window_size).astype(np.int64)
        total[on_x] += _count_peaks_in_windows(ws, we, peak_starts, peak_ends)

    return total


# ===========================================================================
# PERMUTATION HELPERS
# ===========================================================================

def _permute_positions(interactions_df, x_chrom_end, seed):
    """Shuffle X-chromosome anchor positions; preserve anchor sizes."""
    np.random.seed(seed)
    permuted = interactions_df.copy()
    for chrom_col, start_col, end_col in [('chr1', 'start1', 'end1'),
                                           ('chr2', 'start2', 'end2')]:
        on_x = permuted[chrom_col] == 'X'
        if on_x.any():
            sizes      = permuted.loc[on_x, end_col] - permuted.loc[on_x, start_col]
            max_size   = sizes.max()
            new_starts = np.random.randint(0, max(1, x_chrom_end - max_size + 1),
                                           size=on_x.sum())
            new_starts = np.minimum(new_starts, x_chrom_end - sizes.values)
            permuted.loc[on_x, start_col] = new_starts
            permuted.loc[on_x, end_col]   = new_starts + sizes.values
    return permuted


def _x_chrom_end(interactions_df):
    end1 = interactions_df[interactions_df['chr1'] == 'X']['end1'].max() \
           if (interactions_df['chr1'] == 'X').any() else 0
    end2 = interactions_df[interactions_df['chr2'] == 'X']['end2'].max() \
           if (interactions_df['chr2'] == 'X').any() else 0
    return max(end1, end2)


# ---------------------------------------------------------------------------
# Worker functions (module-level so multiprocessing can pickle them)
# ---------------------------------------------------------------------------

def _worker_overlap(seed, interactions_df, peak_starts, peak_ends, window_size, x_chrom_end):
    """Permutation worker: returns mean binary overlap rate (numpy)."""
    permuted    = _permute_positions(interactions_df, x_chrom_end, seed)
    total       = _total_peaks_for_interactions(permuted, peak_starts, peak_ends, window_size)
    return float((total > 0).mean())


def _worker_density(seed, interactions_df, peak_starts, peak_ends, window_size, x_chrom_end):
    """Permutation worker: returns mean peak count per interaction (numpy)."""
    permuted = _permute_positions(interactions_df, x_chrom_end, seed)
    total    = _total_peaks_for_interactions(permuted, peak_starts, peak_ends, window_size)
    return float(total.mean())


def _run_permutations(worker_fn, interactions_df, peak_starts, peak_ends,
                      window_size, x_chrom_end, n_permutations, n_cores):
    seeds = list(range(n_permutations))
    fn    = partial(worker_fn,
                    interactions_df=interactions_df,
                    peak_starts=peak_starts,
                    peak_ends=peak_ends,
                    window_size=window_size,
                    x_chrom_end=x_chrom_end)
    with Pool(processes=n_cores) as pool:
        results = pool.map(fn, seeds)
    return np.array([x for x in results if not np.isnan(x)])


def _perm_stats(observed, null_arr, n_permutations):
    null_mean = null_arr.mean()
    null_std  = null_arr.std()
    z_score   = (observed - null_mean) / null_std if null_std > 0 else 0.0
    p_value   = max(2 * min(np.mean(null_arr >= observed),
                            np.mean(null_arr <= observed)),
                    1 / n_permutations)
    return null_mean, null_std, z_score, p_value


# ===========================================================================
# PERMUTATION TESTS
# ===========================================================================

def permutation_test_overlap(interactions_df, peak_starts, peak_ends, overlap_results,
                              window_size=1000, n_permutations=1000, n_cores=None):
    """
    Permutation test on binary overlap rate.
    NOTE: uninformative when peaks saturate the chromosome (100% overlap).
    Retained for completeness; primary result is the density test.
    """
    print("\n--- Permutation Test 1: Binary Overlap Rate ---")
    observed  = float(overlap_results['any_anchor_overlap'].mean())
    xend      = _x_chrom_end(interactions_df)
    n_cores   = n_cores or max(1, cpu_count() - 1)

    print(f"  Observed overlap rate: {observed*100:.1f}%")
    print(f"  X chrom extent:        {xend:,} bp")
    print(f"  Running {n_permutations} permutations on {n_cores} cores...")

    null = _run_permutations(_worker_overlap, interactions_df, peak_starts, peak_ends,
                             window_size, xend, n_permutations, n_cores)
    null_mean, null_std, z, p = _perm_stats(observed, null, n_permutations)

    print(f"  Null mean: {null_mean*100:.1f}%  std: {null_std*100:.3f}%")
    print(f"  Z={z:.2f}  p={p:.4e}  {'SIGNIFICANT' if p < 0.05 else 'not significant'}")

    return dict(observed_rate=observed, null_distribution=null,
                null_mean=null_mean, null_std=null_std, z_score=z, p_value=p)


def permutation_test_density(interactions_df, peak_starts, peak_ends, overlap_results,
                              window_size=1000, n_permutations=1000, n_cores=None):
    """
    Permutation test on mean H4K16ac peak COUNT per interaction window.
    Robust to peak saturation — tests whether contacts are positioned at
    higher-density H4K16ac regions than random X-chromosome windows.
    """
    print("\n--- Permutation Test 2: Peak Density (all interactions) ---")
    observed = float(overlap_results['total_peaks'].mean())
    xend     = _x_chrom_end(interactions_df)
    n_cores  = n_cores or max(1, cpu_count() - 1)

    print(f"  Observed mean peaks/interaction: {observed:.2f}")
    print(f"  X chrom extent:                  {xend:,} bp")
    print(f"  Running {n_permutations} permutations on {n_cores} cores...")

    null = _run_permutations(_worker_density, interactions_df, peak_starts, peak_ends,
                             window_size, xend, n_permutations, n_cores)
    null_mean, null_std, z, p = _perm_stats(observed, null, n_permutations)

    print(f"  Null mean: {null_mean:.2f}  std: {null_std:.4f}")
    print(f"  Z={z:.2f}  p={p:.4e}  {'SIGNIFICANT' if p < 0.05 else 'not significant'}")

    return dict(observed_mean=observed, null_distribution=null,
                null_mean=null_mean, null_std=null_std, z_score=z, p_value=p)


def permutation_test_density_by_genotype(interactions_df, peak_starts, peak_ends,
                                          overlap_results, window_size=1000,
                                          n_permutations=1000, n_cores=None):
    """
    Density permutation run separately for wMel (logFC>0) and uninf. (logFC<0).
    Each genotype tested against its own null of same-N random X windows.
    """
    print("\n--- Permutation Test 2b: Peak Density (per genotype) ---")

    # Merge logFC into overlap_results for subsetting
    interactions_r = interactions_df.reset_index(drop=True)
    overlap_r      = overlap_results.reset_index(drop=True)
    merged         = interactions_r.copy()
    for col in overlap_r.columns:
        if col != 'interaction_idx':
            merged[col] = overlap_r[col].values

    xend    = _x_chrom_end(interactions_df)
    n_cores = n_cores or max(1, cpu_count() - 1)
    results = {}

    for label, mask in [('wmel', WMEL_MASK), ('uninf', UNINF_MASK)]:
        subset_int  = interactions_df[mask(interactions_df)].copy()
        subset_mer  = merged[mask(merged)]
        observed    = float(subset_mer['total_peaks'].mean())

        print(f"\n  {label} (n={len(subset_int)}): observed mean = {observed:.2f}")
        print(f"  Running {n_permutations} permutations on {n_cores} cores...")

        null = _run_permutations(_worker_density, subset_int, peak_starts, peak_ends,
                                 window_size, xend, n_permutations, n_cores)
        null_mean, null_std, z, p = _perm_stats(observed, null, n_permutations)

        print(f"    Null mean: {null_mean:.2f}  Z={z:.2f}  p={p:.4e}"
              f"  {'SIGNIFICANT' if p < 0.05 else 'not significant'}")

        results[label] = dict(observed_mean=observed, null_distribution=null,
                              null_mean=null_mean, null_std=null_std, z_score=z, p_value=p)

    return results


# ===========================================================================
# GENOTYPE-SPECIFIC ANALYSIS
# ===========================================================================

def ensure_boolean_columns(df, columns):
    df = df.copy()
    for col in columns:
        if col in df.columns:
            df[col] = df[col].fillna(False).astype(bool)
    return df


def _merge_interactions_overlaps(interactions_df, overlap_results):
    """Merge interactions and overlap results by positional alignment."""
    interactions_r = interactions_df.reset_index(drop=True)
    overlap_r      = overlap_results.reset_index(drop=True)
    assert len(interactions_r) == len(overlap_r), \
        f"Length mismatch: {len(interactions_r)} vs {len(overlap_r)}"
    merged = interactions_r.copy()
    for col in overlap_r.columns:
        if col != 'interaction_idx':
            merged[col] = overlap_r[col].values
    bool_cols = ['anchor1_overlap', 'anchor2_overlap',
                 'both_anchors_overlap', 'any_anchor_overlap']
    return ensure_boolean_columns(merged, bool_cols)


def analyze_by_logfc_direction(interactions_df, overlap_results):
    """
    Compare H4K16ac overlap between wMel (logFC>0) and uninf. (logFC<0).
    Returns direction_results dict and merged DataFrame.
    """
    print("\nAnalyzing by logFC direction (wMel: logFC>0, uninf.: logFC<0)...")
    merged = _merge_interactions_overlaps(interactions_df, overlap_results)
    merged['interaction_type'] = np.where(merged['chr1'] == merged['chr2'], 'cis', 'trans')

    jw18_wmel  = merged[WMEL_MASK(merged)]
    jw18_uninf = merged[UNINF_MASK(merged)]

    def _grp_stats(subset):
        return dict(
            n_interactions    = len(subset),
            overlap_rate      = subset['any_anchor_overlap'].mean(),
            both_anchors_rate = subset['both_anchors_overlap'].mean(),
            mean_peaks        = subset['total_peaks'].mean(),
            median_peaks      = subset['total_peaks'].median()
        )

    results = {'jw18_wmel': _grp_stats(jw18_wmel),
               'jw18_uninf': _grp_stats(jw18_uninf)}

    # Fisher's exact: overlap rate
    contingency = [
        [int(jw18_wmel['any_anchor_overlap'].sum()),
         int((~jw18_wmel['any_anchor_overlap']).sum())],
        [int(jw18_uninf['any_anchor_overlap'].sum()),
         int((~jw18_uninf['any_anchor_overlap']).sum())]
    ]
    odds_ratio, fisher_p = stats.fisher_exact(contingency)

    # Mann-Whitney U: peak counts
    if len(jw18_wmel) > 0 and len(jw18_uninf) > 0:
        u_stat, mw_p = stats.mannwhitneyu(
            jw18_wmel['total_peaks'], jw18_uninf['total_peaks'], alternative='two-sided'
        )
        n1, n2        = len(jw18_wmel), len(jw18_uninf)
        rank_biserial = 1 - (2 * u_stat) / (n1 * n2)
    else:
        mw_p, rank_biserial = np.nan, np.nan

    results['comparison'] = dict(fisher_odds_ratio=odds_ratio,
                                 fisher_p_value=fisher_p,
                                 mannwhitneyu_p_value=mw_p,
                                 effect_size_rank_biserial=rank_biserial)

    for lbl, key in [('JW18 wMel  (logFC>0)', 'jw18_wmel'),
                     ('JW18 uninf.(logFC<0)', 'jw18_uninf')]:
        r = results[key]
        print(f"  {lbl} n={r['n_interactions']}: overlap={r['overlap_rate']*100:.1f}%  "
              f"mean peaks={r['mean_peaks']:.1f}  median={r['median_peaks']:.1f}")
    print(f"  Fisher's p={fisher_p:.2e}  MWU p={mw_p:.2e}  r={rank_biserial:.3f}")

    return results, merged


def analyze_peaks_per_interaction(interactions_df, overlap_results):
    """
    Peaks-per-interaction comparison (wMel vs uninf.) with Mann-Whitney U,
    rank-biserial r, and Cohen's d.
    """
    print("\nAnalyzing peaks per interaction...")
    merged = _merge_interactions_overlaps(interactions_df, overlap_results)

    jw18_wmel  = merged[WMEL_MASK(merged)]
    jw18_uninf = merged[UNINF_MASK(merged)]

    def _stats(subset):
        n = len(subset)
        return dict(
            n_interactions        = n,
            total_peaks           = int(subset['total_peaks'].sum()),
            peaks_per_interaction = subset['total_peaks'].sum() / n if n > 0 else 0.0,
            mean_peaks            = subset['total_peaks'].mean(),
            median_peaks          = subset['total_peaks'].median(),
            std_peaks             = subset['total_peaks'].std()
        )

    results = {'jw18_wmel': _stats(jw18_wmel), 'jw18_uninf': _stats(jw18_uninf)}

    if len(jw18_wmel) > 0 and len(jw18_uninf) > 0:
        u_stat, mw_p  = stats.mannwhitneyu(
            jw18_wmel['total_peaks'], jw18_uninf['total_peaks'], alternative='two-sided'
        )
        n1, n2         = len(jw18_wmel), len(jw18_uninf)
        rank_biserial  = 1 - (2 * u_stat) / (n1 * n2)
        s1 = results['jw18_wmel']['std_peaks']
        s2 = results['jw18_uninf']['std_peaks']
        pooled_std = np.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2))
        cohens_d   = ((results['jw18_wmel']['mean_peaks'] -
                       results['jw18_uninf']['mean_peaks']) / pooled_std
                      if pooled_std > 0 else np.nan)
    else:
        mw_p, rank_biserial, cohens_d = np.nan, np.nan, np.nan

    results['comparison'] = dict(mannwhitneyu_p_value=mw_p,
                                 effect_size_rank_biserial=rank_biserial,
                                 cohens_d=cohens_d)

    for lbl, key in [('JW18 wMel  (logFC>0)', 'jw18_wmel'),
                     ('JW18 uninf.(logFC<0)', 'jw18_uninf')]:
        r = results[key]
        print(f"  {lbl}: peaks/int={r['peaks_per_interaction']:.2f}  "
              f"mean={r['mean_peaks']:.2f}  median={r['median_peaks']:.1f}")
    print(f"  MWU p={mw_p:.2e}  r={rank_biserial:.3f}  d={cohens_d:.3f}")

    return results


# ===========================================================================
# VISUALIZATION
# ===========================================================================

def _boxplot(ax, data_list, labels, colors, ylabel, title, p_value=None, annotation=None):
    """Reusable boxplot with strip overlay."""
    bp = ax.boxplot(data_list, positions=list(range(len(data_list))),
                    widths=0.6, patch_artist=True, showfliers=False,
                    boxprops=dict(facecolor='white', edgecolor='black', linewidth=1),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(color='black', linewidth=1),
                    capprops=dict(color='black', linewidth=1))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color); patch.set_alpha(0.6)
    for d, pos, color in zip(data_list, range(len(data_list)), colors):
        if len(d) > 0:
            x = np.random.normal(pos, 0.08, size=len(d))
            ax.scatter(x, d, alpha=0.4, s=6, color=color,
                       edgecolors='none', zorder=3)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel(ylabel, fontsize=6)
    ax.set_title(title, fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    ax.patch.set_alpha(0)
    if p_value is not None and not np.isnan(p_value):
        ax.text(0.95, 0.95, f'p = {p_value:.2e}',
                transform=ax.transAxes, ha='right', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    if annotation:
        ax.text(0.05, 0.95, annotation,
                transform=ax.transAxes, ha='left', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))


def create_visualization(merged_df, direction_results, peaks_per_int_results, output_prefix):
    """Analysis plots: 2×2 set 1 + 1×2 set 2."""
    print("\nCreating analysis visualizations...")

    bool_cols = ['anchor1_overlap', 'anchor2_overlap',
                 'both_anchors_overlap', 'any_anchor_overlap']
    merged_df = ensure_boolean_columns(merged_df, bool_cols)

    color_wmel  = '#09aa4b'
    color_uninf = '#8fcb84'
    colors      = [color_wmel, color_uninf]
    labels      = ['JW18 wMel', 'JW18 uninf.']

    plt.rcParams.update({'font.size': 6})

    wmel_peaks  = merged_df[WMEL_MASK(merged_df)]['total_peaks'].values
    uninf_peaks = merged_df[UNINF_MASK(merged_df)]['total_peaks'].values

    # ---- Set 1: 2×2 ----
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)

    # [0,0] Peak count boxplot
    mwu_p = direction_results['comparison']['mannwhitneyu_p_value']
    _boxplot(axes[0, 0], [wmel_peaks, uninf_peaks], labels, colors,
             'Number of H4K16ac Peaks', 'H4K16ac Peak Count by Genotype', p_value=mwu_p)

    # [0,1] Peaks per interaction boxplot
    ann = (f"Mean:\n  wMel:  {peaks_per_int_results['jw18_wmel']['mean_peaks']:.1f}\n"
           f"  uninf: {peaks_per_int_results['jw18_uninf']['mean_peaks']:.1f}")
    _boxplot(axes[0, 1], [wmel_peaks, uninf_peaks], labels, colors,
             'H4K16ac Peaks per Interaction', 'Peaks per Interaction',
             p_value=peaks_per_int_results['comparison']['mannwhitneyu_p_value'],
             annotation=ann)

    # [1,0] logFC vs peaks scatter
    ax = axes[1, 0]
    ax.patch.set_alpha(0)
    ax.scatter(merged_df[WMEL_MASK(merged_df)]['logFC'],
               merged_df[WMEL_MASK(merged_df)]['total_peaks'],
               c=color_wmel,  alpha=0.4, s=10, label='JW18 wMel',   edgecolors='none')
    ax.scatter(merged_df[UNINF_MASK(merged_df)]['logFC'],
               merged_df[UNINF_MASK(merged_df)]['total_peaks'],
               c=color_uninf, alpha=0.4, s=10, label='JW18 uninf.', edgecolors='none')
    ax.set_xlabel('log2 Fold Change', fontsize=6)
    ax.set_ylabel('H4K16ac Peaks', fontsize=6)
    ax.set_title('H4K16ac Enrichment vs Effect Size', fontsize=7, fontweight='bold')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.legend(fontsize=5); ax.tick_params(labelsize=6)
    try:
        rho, pval = stats.spearmanr(merged_df['logFC'], merged_df['total_peaks'])
        ax.text(0.05, 0.95, f'ρ={rho:.3f}\np={pval:.2e}',
                transform=ax.transAxes, ha='left', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    except Exception:
        pass

    # [1,1] Overlap pattern bar
    ax = axes[1, 1]
    ax.patch.set_alpha(0)
    pattern_counts = [
        int((~merged_df['any_anchor_overlap']).sum()),
        int((merged_df['any_anchor_overlap'] & ~merged_df['both_anchors_overlap']).sum()),
        int(merged_df['both_anchors_overlap'].sum())
    ]
    bars = ax.bar([0, 1, 2], pattern_counts,
                  color=['#cccccc', '#ff9933', '#cc3333'], alpha=0.8, width=0.6)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['No overlap', 'One anchor', 'Both anchors'], fontsize=6)
    ax.set_ylabel('Number of Interactions', fontsize=6)
    ax.set_title('H4K16ac Overlap Pattern', fontsize=7, fontweight='bold')
    ax.tick_params(axis='x', rotation=30, labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    for i, (bar, cnt) in enumerate(zip(bars, pattern_counts)):
        ax.text(i, bar.get_height(), str(cnt), ha='center', va='bottom', fontsize=5)

    plt.tight_layout()
    out1 = f"{output_prefix}_analysis_set1.pdf"
    plt.savefig(out1, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"  Saved {out1}")

    # ---- Set 2: 1×2 ----
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    fig.patch.set_alpha(0)

    # [0] Overlap rate bar
    ax = axes[0]
    ax.patch.set_alpha(0)
    overlap_rates = [direction_results['jw18_wmel']['overlap_rate'] * 100,
                     direction_results['jw18_uninf']['overlap_rate'] * 100]
    bars = ax.bar([0, 1], overlap_rates, color=colors, alpha=0.8, width=0.6)
    ax.set_xticks([0, 1]); ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel('H4K16ac Overlap Rate (%)', fontsize=6)
    ax.set_title('H4K16ac Overlap Rate by Genotype', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    fp = direction_results['comparison']['fisher_p_value']
    if fp is not None:
        ax.text(0.5, max(overlap_rates) * 0.95, f'Fisher p={fp:.2e}',
                ha='center', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    for i, (bar, rate) in enumerate(zip(bars, overlap_rates)):
        ax.text(i, bar.get_height(), f'{rate:.1f}%', ha='center', va='bottom', fontsize=5)

    # [1] Per-anchor boxplot
    ax = axes[1]
    ax.patch.set_alpha(0)
    wmel_a1  = merged_df[WMEL_MASK(merged_df)]['n_peaks_anchor1'].values
    wmel_a2  = merged_df[WMEL_MASK(merged_df)]['n_peaks_anchor2'].values
    uninf_a1 = merged_df[UNINF_MASK(merged_df)]['n_peaks_anchor1'].values
    uninf_a2 = merged_df[UNINF_MASK(merged_df)]['n_peaks_anchor2'].values
    if all(len(a) > 0 for a in [wmel_a1, wmel_a2, uninf_a1, uninf_a2]):
        bp = ax.boxplot([wmel_a1, wmel_a2, uninf_a1, uninf_a2],
                        positions=[0, 1, 2.5, 3.5], widths=0.6,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor='white', edgecolor='black', linewidth=1),
                        medianprops=dict(color='black', linewidth=1.5),
                        whiskerprops=dict(color='black', linewidth=1),
                        capprops=dict(color='black', linewidth=1))
        anchor_colors = [color_wmel, color_wmel, color_uninf, color_uninf]
        for patch, color in zip(bp['boxes'], anchor_colors):
            patch.set_facecolor(color); patch.set_alpha(0.6)
        for d, pos, color in zip([wmel_a1, wmel_a2, uninf_a1, uninf_a2],
                                  [0, 1, 2.5, 3.5], anchor_colors):
            x = np.random.normal(pos, 0.08, size=len(d))
            ax.scatter(x, d, alpha=0.4, s=6, color=color, edgecolors='none', zorder=3)
        ax.set_xticks([0.5, 3]); ax.set_xticklabels(['JW18 wMel', 'JW18 uninf.'], fontsize=6)
        ax.set_ylabel('H4K16ac Peaks per Anchor', fontsize=6)
        ax.set_title('Per-Anchor H4K16ac Peaks', fontsize=7, fontweight='bold')
        ax.tick_params(labelsize=6)

    plt.tight_layout()
    out2 = f"{output_prefix}_analysis_set2.pdf"
    plt.savefig(out2, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"  Saved {out2}")


def create_permutation_visualization(perm_overlap, perm_density, perm_density_genotype,
                                     output_prefix):
    """Three-panel permutation null distribution plots."""
    print("\nCreating permutation visualizations...")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.patch.set_alpha(0)

    # Panel 1: binary overlap
    ax = axes[0]
    ax.patch.set_alpha(0)
    ax.hist(perm_overlap['null_distribution'] * 100, bins=30,
            color='gray', alpha=0.7, edgecolor='black', linewidth=0.5, label='Null')
    ax.axvline(perm_overlap['observed_rate'] * 100, color='red', linestyle='--',
               linewidth=2, label=f"Obs ({perm_overlap['observed_rate']*100:.1f}%)")
    ax.axvline(perm_overlap['null_mean'] * 100, color='blue', linestyle=':',
               linewidth=1.5, label=f"Null mean ({perm_overlap['null_mean']*100:.1f}%)")
    ax.set_xlabel('H4K16ac Overlap Rate (%)', fontsize=8)
    ax.set_ylabel('Frequency', fontsize=8)
    ax.set_title('Perm Test 1: Binary Overlap Rate', fontsize=9, fontweight='bold')
    ax.legend(fontsize=6)
    ax.text(0.05, 0.95,
            f"Z={perm_overlap['z_score']:.2f}\np={perm_overlap['p_value']:.2e}",
            transform=ax.transAxes, ha='left', va='top', fontsize=7,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Panel 2: density (all)
    ax = axes[1]
    ax.patch.set_alpha(0)
    ax.hist(perm_density['null_distribution'], bins=30,
            color='gray', alpha=0.7, edgecolor='black', linewidth=0.5, label='Null')
    ax.axvline(perm_density['observed_mean'], color='red', linestyle='--',
               linewidth=2, label=f"Obs ({perm_density['observed_mean']:.1f})")
    ax.axvline(perm_density['null_mean'], color='blue', linestyle=':',
               linewidth=1.5, label=f"Null mean ({perm_density['null_mean']:.1f})")
    ax.set_xlabel('Mean H4K16ac Peaks per Interaction', fontsize=8)
    ax.set_ylabel('Frequency', fontsize=8)
    ax.set_title('Perm Test 2: Peak Density (all)', fontsize=9, fontweight='bold')
    ax.legend(fontsize=6)
    ax.text(0.05, 0.95,
            f"Z={perm_density['z_score']:.2f}\np={perm_density['p_value']:.2e}",
            transform=ax.transAxes, ha='left', va='top', fontsize=7,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Panel 3: density per genotype
    ax = axes[2]
    ax.patch.set_alpha(0)
    palette = {'wmel': '#09aa4b', 'uninf': '#8fcb84'}
    for key, color in palette.items():
        r    = perm_density_genotype[key]
        lbl  = 'wMel' if key == 'wmel' else 'uninf.'
        ax.hist(r['null_distribution'], bins=30, alpha=0.5, color=color,
                edgecolor='black', linewidth=0.3, label=f'{lbl} null')
        ax.axvline(r['observed_mean'], color=color, linestyle='--', linewidth=2,
                   label=f'{lbl} obs ({r["observed_mean"]:.0f})')
    ax.set_xlabel('Mean H4K16ac Peaks per Interaction', fontsize=8)
    ax.set_ylabel('Frequency', fontsize=8)
    ax.set_title('Perm Test 2b: Peak Density (per genotype)', fontsize=9, fontweight='bold')
    ax.legend(fontsize=6)
    ptext = (f"wMel:   Z={perm_density_genotype['wmel']['z_score']:.2f}, "
             f"p={perm_density_genotype['wmel']['p_value']:.2e}\n"
             f"uninf.: Z={perm_density_genotype['uninf']['z_score']:.2f}, "
             f"p={perm_density_genotype['uninf']['p_value']:.2e}")
    ax.text(0.05, 0.95, ptext, transform=ax.transAxes, ha='left', va='top', fontsize=6,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    out = f"{output_prefix}_permutation_tests.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"  Saved {out}")


# ===========================================================================
# OUTPUT / SUMMARY
# ===========================================================================

def save_summary(output_prefix, args, all_x_interactions, sig_interactions,
                 h4k16ac_peaks, perm_overlap, perm_density, perm_density_genotype,
                 direction_results, peaks_per_int_results):

    summary = {
        'total_x_interactions':             len(all_x_interactions),
        'significant_x_interactions':       len(sig_interactions),
        'h4k16ac_peaks_on_x':               len(h4k16ac_peaks),
        # perm 1: overlap
        'perm_overlap_observed_rate':        perm_overlap['observed_rate'],
        'perm_overlap_null_mean':            perm_overlap['null_mean'],
        'perm_overlap_null_std':             perm_overlap['null_std'],
        'perm_overlap_z':                    perm_overlap['z_score'],
        'perm_overlap_pvalue':               perm_overlap['p_value'],
        # perm 2: density all
        'perm_density_observed_mean':        perm_density['observed_mean'],
        'perm_density_null_mean':            perm_density['null_mean'],
        'perm_density_null_std':             perm_density['null_std'],
        'perm_density_z':                    perm_density['z_score'],
        'perm_density_pvalue':               perm_density['p_value'],
        # perm 2b: density per genotype
        'perm_density_wmel_observed':        perm_density_genotype['wmel']['observed_mean'],
        'perm_density_wmel_null_mean':       perm_density_genotype['wmel']['null_mean'],
        'perm_density_wmel_z':               perm_density_genotype['wmel']['z_score'],
        'perm_density_wmel_pvalue':          perm_density_genotype['wmel']['p_value'],
        'perm_density_uninf_observed':       perm_density_genotype['uninf']['observed_mean'],
        'perm_density_uninf_null_mean':      perm_density_genotype['uninf']['null_mean'],
        'perm_density_uninf_z':              perm_density_genotype['uninf']['z_score'],
        'perm_density_uninf_pvalue':         perm_density_genotype['uninf']['p_value'],
        # genotype stats
        'jw18_wmel_n':                       direction_results['jw18_wmel']['n_interactions'],
        'jw18_wmel_overlap_rate':            direction_results['jw18_wmel']['overlap_rate'],
        'jw18_wmel_mean_peaks':              direction_results['jw18_wmel']['mean_peaks'],
        'jw18_wmel_median_peaks':            direction_results['jw18_wmel']['median_peaks'],
        'jw18_uninf_n':                      direction_results['jw18_uninf']['n_interactions'],
        'jw18_uninf_overlap_rate':           direction_results['jw18_uninf']['overlap_rate'],
        'jw18_uninf_mean_peaks':             direction_results['jw18_uninf']['mean_peaks'],
        'jw18_uninf_median_peaks':           direction_results['jw18_uninf']['median_peaks'],
        'overlap_fisher_pvalue':             direction_results['comparison']['fisher_p_value'],
        'peak_count_mwu_pvalue':             direction_results['comparison']['mannwhitneyu_p_value'],
        'peak_count_effect_size_rb':         direction_results['comparison']['effect_size_rank_biserial'],
        # peaks per interaction
        'jw18_wmel_peaks_per_interaction':   peaks_per_int_results['jw18_wmel']['peaks_per_interaction'],
        'jw18_uninf_peaks_per_interaction':  peaks_per_int_results['jw18_uninf']['peaks_per_interaction'],
        'peaks_per_int_mwu_pvalue':          peaks_per_int_results['comparison']['mannwhitneyu_p_value'],
        'peaks_per_int_effect_size_rb':      peaks_per_int_results['comparison']['effect_size_rank_biserial'],
        'peaks_per_int_cohens_d':            peaks_per_int_results['comparison']['cohens_d'],
    }

    pd.DataFrame([summary]).to_csv(f"{output_prefix}_summary.csv", index=False)
    print(f"Summary CSV saved to {output_prefix}_summary.csv")

    with open(f"{output_prefix}_summary.txt", 'w') as f:
        f.write("H4K16ac ENRICHMENT ANALYSIS SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        f.write("logFC convention: >0 = JW18 wMel, <0 = JW18 uninf.\n\n")
        f.write(f"Parameters:\n"
                f"  window_size={args.window_size} bp  "
                f"fdr={args.fdr_threshold}  "
                f"n_permutations={args.n_permutations}\n\n")
        f.write(f"Data:\n"
                f"  Total interactions:       {summary['total_x_interactions']}\n"
                f"  Significant interactions: {summary['significant_x_interactions']}\n"
                f"  H4K16ac peaks on X:       {summary['h4k16ac_peaks_on_x']}\n\n")

        f.write("Permutation Test 1: Binary Overlap Rate\n")
        f.write(f"  Observed: {summary['perm_overlap_observed_rate']*100:.1f}%  "
                f"Null: {summary['perm_overlap_null_mean']*100:.1f}% "
                f"(std={summary['perm_overlap_null_std']*100:.3f}%)\n"
                f"  Z={summary['perm_overlap_z']:.2f}  p={summary['perm_overlap_pvalue']:.4e}\n\n")

        f.write("Permutation Test 2: Peak Density (all interactions)\n")
        f.write(f"  Observed mean: {summary['perm_density_observed_mean']:.2f}  "
                f"Null: {summary['perm_density_null_mean']:.2f} "
                f"(std={summary['perm_density_null_std']:.4f})\n"
                f"  Z={summary['perm_density_z']:.2f}  p={summary['perm_density_pvalue']:.4e}\n\n")

        f.write("Permutation Test 2b: Peak Density (per genotype)\n")
        f.write(f"  wMel:   obs={summary['perm_density_wmel_observed']:.2f}  "
                f"null={summary['perm_density_wmel_null_mean']:.2f}  "
                f"Z={summary['perm_density_wmel_z']:.2f}  "
                f"p={summary['perm_density_wmel_pvalue']:.4e}\n"
                f"  uninf.: obs={summary['perm_density_uninf_observed']:.2f}  "
                f"null={summary['perm_density_uninf_null_mean']:.2f}  "
                f"Z={summary['perm_density_uninf_z']:.2f}  "
                f"p={summary['perm_density_uninf_pvalue']:.4e}\n\n")

        f.write("Genotype Analysis:\n")
        for lbl, key in [('JW18 wMel  (logFC>0)', 'jw18_wmel'),
                         ('JW18 uninf.(logFC<0)', 'jw18_uninf')]:
            f.write(f"  {lbl}:\n"
                    f"    N={summary[key+'_n']}  "
                    f"overlap={summary[key+'_overlap_rate']*100:.1f}%  "
                    f"mean peaks={summary[key+'_mean_peaks']:.2f}  "
                    f"median={summary[key+'_median_peaks']:.1f}\n")
        f.write(f"  Fisher p={summary['overlap_fisher_pvalue']:.2e}  "
                f"MWU p={summary['peak_count_mwu_pvalue']:.2e}  "
                f"r={summary['peak_count_effect_size_rb']:.3f}\n\n")

        f.write("Peaks per Interaction:\n"
                f"  wMel:   {summary['jw18_wmel_peaks_per_interaction']:.2f}\n"
                f"  uninf.: {summary['jw18_uninf_peaks_per_interaction']:.2f}\n"
                f"  MWU p={summary['peaks_per_int_mwu_pvalue']:.2e}  "
                f"r={summary['peaks_per_int_effect_size_rb']:.3f}  "
                f"d={summary['peaks_per_int_cohens_d']:.3f}\n")

    print(f"Text summary saved to {output_prefix}_summary.txt")
    return summary


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description='H4K16ac enrichment at differential X-chromosome contacts (fast numpy permutations)'
    )
    parser.add_argument('--interactions',   required=True,
                        help='Differential interactions CSV (from diffHic)')
    parser.add_argument('--h4k16ac_peaks',  required=True,
                        help='H4K16ac ChIP-seq peaks BED file')
    parser.add_argument('--window_size',    type=int,   default=1000,
                        help='Window around anchors in bp (default: 1000)')
    parser.add_argument('--fdr_threshold',  type=float, default=0.05,
                        help='FDR threshold (default: 0.05)')
    parser.add_argument('--n_permutations', type=int,   default=1000,
                        help='Number of permutations (default: 1000)')
    parser.add_argument('--n_cores',        type=int,   default=None,
                        help='CPU cores for parallelism (default: all-1)')
    parser.add_argument('--skip_overlap_perm', action='store_true',
                        help='Skip binary overlap permutation (uninformative when saturated)')
    parser.add_argument('--output_prefix',  required=True,
                        help='Output prefix (including directory)')
    args = parser.parse_args()

    Path(args.output_prefix).parent.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("H4K16ac Enrichment Analysis — Fast Numpy Permutations")
    print("logFC: >0 = wMel, <0 = uninf.")
    print("=" * 60)

    # Load data
    h4k16ac_peaks                     = load_h4k16ac_peaks(args.h4k16ac_peaks)
    sig_interactions, all_interactions = load_differential_interactions(
        args.interactions, args.fdr_threshold
    )
    if h4k16ac_peaks is None or len(sig_interactions) == 0:
        print("ERROR: no data to analyze"); return 1

    # Observed overlaps (pybedtools, run once)
    print(f"\nCalculating observed H4K16ac overlaps (window={args.window_size} bp)...")
    overlap_results = calculate_h4k16ac_overlap_vectorized(
        sig_interactions, h4k16ac_peaks, args.window_size
    )
    if overlap_results is None:
        print("ERROR: overlap calculation failed"); return 1

    print(f"\n  Any anchor overlapping: "
          f"{overlap_results['any_anchor_overlap'].sum()} "
          f"({overlap_results['any_anchor_overlap'].mean()*100:.1f}%)")
    print(f"  Both anchors overlapping: "
          f"{overlap_results['both_anchors_overlap'].sum()} "
          f"({overlap_results['both_anchors_overlap'].mean()*100:.1f}%)")
    print(f"  Mean peaks/interaction: {overlap_results['total_peaks'].mean():.2f}")

    # Build peak index once — shared across all permutation workers
    print("\nBuilding numpy peak index for fast permutations...")
    peak_starts, peak_ends = build_peak_index(h4k16ac_peaks)
    print(f"  Index built: {len(peak_starts):,} peaks")

    # Genotype-level analyses
    direction_results, merged_df = analyze_by_logfc_direction(
        sig_interactions, overlap_results
    )
    peaks_per_int_results = analyze_peaks_per_interaction(
        sig_interactions, overlap_results
    )

    # Permutation tests
    n_cores = args.n_cores or max(1, cpu_count() - 1)
    print(f"\nUsing {n_cores} cores for permutations")

    if args.skip_overlap_perm:
        print("\nSkipping binary overlap permutation (--skip_overlap_perm)")
        perm_overlap = dict(observed_rate=float(overlap_results['any_anchor_overlap'].mean()),
                            null_distribution=np.array([]), null_mean=np.nan,
                            null_std=np.nan, z_score=np.nan, p_value=np.nan)
    else:
        perm_overlap = permutation_test_overlap(
            sig_interactions, peak_starts, peak_ends, overlap_results,
            window_size=args.window_size,
            n_permutations=args.n_permutations,
            n_cores=n_cores
        )

    perm_density = permutation_test_density(
        sig_interactions, peak_starts, peak_ends, overlap_results,
        window_size=args.window_size,
        n_permutations=args.n_permutations,
        n_cores=n_cores
    )

    perm_density_genotype = permutation_test_density_by_genotype(
        sig_interactions, peak_starts, peak_ends, overlap_results,
        window_size=args.window_size,
        n_permutations=args.n_permutations,
        n_cores=n_cores
    )

    # Outputs
    merged_df.to_csv(f"{args.output_prefix}_interactions.csv", index=False)
    print(f"\nDetailed results: {args.output_prefix}_interactions.csv")

    save_summary(args.output_prefix, args, all_interactions, sig_interactions,
                 h4k16ac_peaks, perm_overlap, perm_density, perm_density_genotype,
                 direction_results, peaks_per_int_results)

    create_visualization(merged_df, direction_results, peaks_per_int_results,
                         args.output_prefix)
    create_permutation_visualization(perm_overlap, perm_density, perm_density_genotype,
                                     args.output_prefix)

    print(f"\nDone. All outputs at: {args.output_prefix}")
    return 0


if __name__ == '__main__':
    exit(main())