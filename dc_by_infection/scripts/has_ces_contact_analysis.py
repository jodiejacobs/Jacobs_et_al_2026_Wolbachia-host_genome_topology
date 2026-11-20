#!/usr/bin/env python3
"""
Analyze High-Affinity Sites (HAS) / Chromatin Entry Sites (CES) enrichment 
at differential chromatin contacts on the X chromosome.

HAS/CES are specific DNA sequence motifs where the MSL dosage compensation 
complex initially binds, then spreads along the X chromosome.

Based on:
- Straub et al. (2008) - CES motif discovery
- Alekseyenko et al. (2008) - HAS identification  
- Villa et al. (2016) - HAS network facilitates MSL spreading
"""

import pandas as pd
import numpy as np
import pybedtools
from scipy import stats
from statsmodels.stats.multitest import multipletests
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def load_has_sites(bed_file):
    """Load HAS/CES sites from BED file. Format: chrX start end"""
    print(f"Loading HAS/CES sites from {bed_file}")
    try:
        sites = []
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    chrom = parts[0].replace('chr', '')
                    sites.append({
                        'chrom': chrom,
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'size': int(parts[2]) - int(parts[1])
                    })
        sites_df = pd.DataFrame(sites)
        print(f"Total HAS/CES sites loaded: {len(sites_df)}")
        print(f"Chromosomes present: {sorted(sites_df['chrom'].unique())}")
        x_sites = sites_df[sites_df['chrom'] == 'X'].copy()
        if len(x_sites) == 0:
            print("Warning: No HAS/CES sites found on X chromosome")
            print(f"Available chromosomes: {sites_df['chrom'].unique()}")
            return None
        print(f"X chromosome HAS/CES sites: {len(x_sites)}")
        print(f"Size range: {x_sites['size'].min()}-{x_sites['size'].max()} bp")
        print(f"Genomic span: {x_sites['start'].min():,} - {x_sites['end'].max():,}")
        return x_sites
    except Exception as e:
        print(f"Error loading HAS/CES sites: {e}")
        import traceback
        traceback.print_exc()
        return None


def load_differential_interactions(interactions_file, fdr_threshold=0.05):
    """Load significant differential interactions"""
    print(f"\nLoading differential interactions from {interactions_file}")
    interactions = pd.read_csv(interactions_file)
    if 'chr1' in interactions.columns:
        interactions['chr1'] = interactions['chr1'].astype(str).str.replace('chr', '')
        interactions['chr2'] = interactions['chr2'].astype(str).str.replace('chr', '')
    print(f"Loaded {len(interactions)} total interactions")
    x_interactions = interactions[
        ((interactions['chr1'] == 'X') | (interactions['chr2'] == 'X'))
    ].copy()
    print(f"Found {len(x_interactions)} total X chromosome interactions")
    sig_x = x_interactions[
        (x_interactions['FDR'] < fdr_threshold) &
        (abs(x_interactions['logFC']) > 1)
    ].copy()
    print(f"Found {len(sig_x)} significant X chromosome interactions")
    print(f"  FDR < {fdr_threshold}, |logFC| > 1")
    if len(sig_x) > 0:
        print(f"  Cis interactions: {sum(sig_x['chr1'] == sig_x['chr2'])}")
        print(f"  Trans interactions: {sum(sig_x['chr1'] != sig_x['chr2'])}")
        print(f"  JW18 uninf. (down-regulated): {sum(sig_x['logFC'] > 0)}")
        print(f"  JW18 wMel (up-regulated): {sum(sig_x['logFC'] < 0)}")
    return sig_x, x_interactions


def load_null_model(null_file):
    """Load null model interactions"""
    print(f"\nLoading null model from {null_file}")
    null_data = pd.read_csv(null_file)
    print(f"Loaded {len(null_data)} null model interactions")
    return null_data


def calculate_has_overlap(interactions_df, has_sites, window_size=50000):
    """
    Calculate overlap between interaction anchors and HAS/CES sites.
    Uses larger window (50kb default) as MSL complex spreads from HAS sites.
    """
    print(f"\nAnalyzing HAS/CES overlap (window size: {window_size}bp)")
    print("Note: Large window accounts for MSL spreading from HAS sites")
    if has_sites is None or len(has_sites) == 0:
        print("No HAS/CES sites available")
        return None

    has_extended = has_sites.copy()
    has_extended['start'] = (has_extended['start'] - window_size).clip(lower=0)
    has_extended['end'] = has_extended['end'] + window_size
    has_bt = pybedtools.BedTool.from_dataframe(has_extended[['chrom', 'start', 'end']])

    results = []
    for idx, interaction in interactions_df.iterrows():
        anchor1_chr   = interaction['chr1']
        anchor1_start = interaction['start1']
        anchor1_end   = interaction['end1']
        anchor2_chr   = interaction['chr2']
        anchor2_start = interaction['start2']
        anchor2_end   = interaction['end2']

        try:
            anchor1_bt = pybedtools.BedTool.from_dataframe(pd.DataFrame([{
                'chrom': anchor1_chr, 'start': anchor1_start, 'end': anchor1_end}]))
            anchor2_bt = pybedtools.BedTool.from_dataframe(pd.DataFrame([{
                'chrom': anchor2_chr, 'start': anchor2_start, 'end': anchor2_end}]))

            n_has_anchor1 = len(list(anchor1_bt.intersect(has_bt, wa=True)))
            n_has_anchor2 = len(list(anchor2_bt.intersect(has_bt, wa=True)))

            # Distance to nearest HAS
            if anchor1_chr == 'X':
                d1 = [min(abs(anchor1_start - r['end']), abs(anchor1_end - r['start']))
                      for _, r in has_sites.iterrows()]
                min_dist_anchor1 = min(d1) if d1 else np.inf
            else:
                min_dist_anchor1 = np.inf

            if anchor2_chr == 'X':
                d2 = [min(abs(anchor2_start - r['end']), abs(anchor2_end - r['start']))
                      for _, r in has_sites.iterrows()]
                min_dist_anchor2 = min(d2) if d2 else np.inf
            else:
                min_dist_anchor2 = np.inf

            results.append({
                'interaction_idx': idx,
                'anchor1_overlap': n_has_anchor1 > 0,
                'anchor2_overlap': n_has_anchor2 > 0,
                'any_anchor_overlap': (n_has_anchor1 > 0) or (n_has_anchor2 > 0),
                'both_anchors_overlap': (n_has_anchor1 > 0) and (n_has_anchor2 > 0),
                'n_has_anchor1': n_has_anchor1,
                'n_has_anchor2': n_has_anchor2,
                'total_has': n_has_anchor1 + n_has_anchor2,
                'min_dist_anchor1': min_dist_anchor1,
                'min_dist_anchor2': min_dist_anchor2,
                'min_dist_any': min(min_dist_anchor1, min_dist_anchor2)
            })
        except Exception as e:
            print(f"Warning: Could not process interaction {idx}: {e}")
            continue

    results_df = pd.DataFrame(results)
    print(f"\nHAS/CES Overlap Summary:")
    print(f"  Interactions with any anchor overlap: {results_df['any_anchor_overlap'].sum()} ({results_df['any_anchor_overlap'].mean()*100:.1f}%)")
    print(f"  Interactions with both anchors overlap: {results_df['both_anchors_overlap'].sum()} ({results_df['both_anchors_overlap'].mean()*100:.1f}%)")
    print(f"  Mean HAS sites per interaction: {results_df['total_has'].mean():.2f}")
    print(f"  Median distance to nearest HAS: {results_df['min_dist_any'].median():.0f} bp")
    return results_df


def analyze_distance_distribution(merged_df, has_sites):
    """Analyze the distribution of distances from interactions to HAS sites."""
    print("\nAnalyzing distance distribution to HAS sites...")
    sig_distances = merged_df['min_dist_any'].values
    sig_distances = sig_distances[~np.isinf(sig_distances)]
    if len(sig_distances) == 0:
        print("No valid distances to analyze")
        return None
    quartiles = np.percentile(sig_distances, [25, 50, 75])
    distance_bins = [0, 10000, 25000, 50000, 100000, np.inf]
    bin_labels = ['0-10kb', '10-25kb', '25-50kb', '50-100kb', '>100kb']
    bin_counts = [np.sum((sig_distances >= lo) & (sig_distances < hi))
                  for lo, hi in zip(distance_bins, distance_bins[1:])]
    distance_dist = pd.DataFrame({
        'distance_range': bin_labels,
        'count': bin_counts,
        'percentage': np.array(bin_counts) / len(sig_distances) * 100
    })
    print("Distance distribution:")
    for _, row in distance_dist.iterrows():
        print(f"  {row['distance_range']}: {row['count']} ({row['percentage']:.1f}%)")
    return {
        'quartiles': quartiles,
        'distance_distribution': distance_dist,
        'median_distance': np.median(sig_distances),
        'mean_distance': np.mean(sig_distances)
    }


def compare_to_null_model(real_overlap_rate, null_data, n_interactions, n_overlaps):
    """Compare real HAS/CES overlap rate to null expectation."""
    print(f"\nComparing to null model...")
    null_expectation = 0.075
    enrichment = real_overlap_rate / null_expectation if null_expectation > 0 else np.inf
    from scipy.stats import binom
    log_p_value = binom.logsf(n_overlaps - 1, n_interactions, null_expectation)
    p_value  = np.exp(log_p_value)
    log10_p  = log_p_value / np.log(10)
    print(f"  Real overlap rate: {real_overlap_rate*100:.1f}%")
    print(f"  Expected (null): {null_expectation*100:.1f}%")
    print(f"  Enrichment: {enrichment:.2f}x")
    if p_value == 0.0:
        if np.isinf(log10_p):
            print(f"  P-value: vanishingly small (log10 p = -inf, extreme float underflow)")
        else:
            print(f"  P-value: < 1e{int(np.floor(log10_p))} (log10 p = {log10_p:.1f}, float underflow)")
    else:
        print(f"  P-value: {p_value:.2e}")
    return {
        'enrichment': enrichment,
        'p_value': p_value,
        'log10_p_value': log10_p,
        'real_overlap_rate': real_overlap_rate,
        'null_expectation': null_expectation
    }


def analyze_by_logfc_direction(interactions_df, overlap_df):
    """Analyze HAS/CES enrichment separately for JW18 uninf. and JW18 wMel interactions."""
    print("\nAnalyzing by logFC direction...")

    merged = interactions_df.reset_index(drop=True).copy()
    overlap_clean = overlap_df.reset_index(drop=True)
    for col in overlap_clean.columns:
        if col != 'interaction_idx':
            merged[col] = overlap_clean[col].values

    # Ensure boolean columns
    for col in ['anchor1_overlap', 'anchor2_overlap', 'any_anchor_overlap', 'both_anchors_overlap']:
        if col in merged.columns:
            merged[col] = merged[col].fillna(False).astype(bool)

    # logFC > 0 → more contact in JW18 uninf. (DOX); logFC < 0 → more in JW18 wMel
    jw18_uninf = merged[merged['logFC'] > 0]
    jw18_wmel  = merged[merged['logFC'] < 0]

    def _safe_stat(arr):
        return {
            'mean': arr.mean() if len(arr) > 0 else 0,
            'median': arr.median() if len(arr) > 0 else 0,
            'std': arr.std() if len(arr) > 0 else 0,
        }

    uninf_has_s = _safe_stat(jw18_uninf['total_has'])
    wmel_has_s  = _safe_stat(jw18_wmel['total_has'])

    results = {
        'jw18_uninf': {
            'n_interactions': len(jw18_uninf),
            'overlap_rate': jw18_uninf['any_anchor_overlap'].mean() if len(jw18_uninf) > 0 else 0,
            'mean_has':     uninf_has_s['mean'],
            'median_has':   uninf_has_s['median'],
            'std_has':      uninf_has_s['std'],
            'median_distance': jw18_uninf['min_dist_any'].median() if len(jw18_uninf) > 0 else np.inf
        },
        'jw18_wmel': {
            'n_interactions': len(jw18_wmel),
            'overlap_rate': jw18_wmel['any_anchor_overlap'].mean() if len(jw18_wmel) > 0 else 0,
            'mean_has':     wmel_has_s['mean'],
            'median_has':   wmel_has_s['median'],
            'std_has':      wmel_has_s['std'],
            'median_distance': jw18_wmel['min_dist_any'].median() if len(jw18_wmel) > 0 else np.inf
        }
    }

    if len(jw18_uninf) > 0 and len(jw18_wmel) > 0:
        # Chi-square: overlap rates
        contingency = np.array([
            [jw18_uninf['any_anchor_overlap'].sum(),  (~jw18_uninf['any_anchor_overlap']).sum()],
            [jw18_wmel['any_anchor_overlap'].sum(),   (~jw18_wmel['any_anchor_overlap']).sum()]
        ])
        chi2, chi2_p = stats.chi2_contingency(contingency)[:2]

        # Mann-Whitney: distances
        uninf_dists = jw18_uninf['min_dist_any'].values
        wmel_dists  = jw18_wmel['min_dist_any'].values
        uninf_dists = uninf_dists[~np.isinf(uninf_dists)]
        wmel_dists  = wmel_dists[~np.isinf(wmel_dists)]
        if len(uninf_dists) > 0 and len(wmel_dists) > 0:
            mw_stat, mw_p = stats.mannwhitneyu(uninf_dists, wmel_dists, alternative='two-sided')
        else:
            mw_stat, mw_p = np.nan, 1.0

        # Mann-Whitney: HAS hit counts
        try:
            has_u, has_mw_p = stats.mannwhitneyu(
                jw18_uninf['total_has'].values,
                jw18_wmel['total_has'].values,
                alternative='two-sided'
            )
            n1, n2 = len(jw18_uninf), len(jw18_wmel)
            has_rb = 1 - (2 * has_u) / (n1 * n2)
        except Exception:
            has_mw_p, has_rb = np.nan, np.nan

        # Welch's t-test: HAS hit counts (tests mean difference, robust to unequal N/variance)
        try:
            t_stat, t_p = stats.ttest_ind(
                jw18_uninf['total_has'].values,
                jw18_wmel['total_has'].values,
                equal_var=False,
                alternative='two-sided'
            )
            pooled_std = np.sqrt(
                ((len(jw18_uninf) - 1) * uninf_has_s['std']**2 +
                 (len(jw18_wmel)  - 1) * wmel_has_s['std']**2) /
                (len(jw18_uninf) + len(jw18_wmel) - 2)
            )
            cohens_d = ((uninf_has_s['mean'] - wmel_has_s['mean']) / pooled_std
                        if pooled_std > 0 else np.nan)
        except Exception:
            t_stat, t_p, cohens_d = np.nan, np.nan, np.nan

        # One-sided binomial test: is wMel HAS overlap rate > uninf rate?
        # H0: p_wmel = p_uninf; H1: p_wmel > p_uninf
        try:
            from scipy.stats import binomtest
            p_null_binom = results['jw18_uninf']['overlap_rate']
            n_wmel       = results['jw18_wmel']['n_interactions']
            k_wmel       = int(jw18_wmel['any_anchor_overlap'].sum())
            binom_result = binomtest(k_wmel, n_wmel, p_null_binom, alternative='greater')
            binom_p      = binom_result.pvalue
        except Exception:
            binom_p, k_wmel, n_wmel, p_null_binom = np.nan, np.nan, np.nan, np.nan

        results['comparison'] = {
            'chi2': chi2,
            'chi2_pvalue': chi2_p,
            'mannwhitney_stat': mw_stat,
            'mannwhitney_pvalue': mw_p,
            'has_count_mannwhitney_pvalue': has_mw_p,
            'has_count_effect_size': has_rb,
            'has_count_ttest_pvalue': t_p,
            'has_count_ttest_stat': t_stat,
            'has_count_cohens_d': cohens_d,
            'binom_pvalue': binom_p,
            'binom_k': k_wmel,
            'binom_n': n_wmel,
            'binom_null_prop': p_null_binom
        }

        print(f"\nJW18 uninf. (up-regulated) interactions:")
        print(f"  N = {results['jw18_uninf']['n_interactions']}")
        print(f"  HAS overlap rate: {results['jw18_uninf']['overlap_rate']*100:.1f}%")
        print(f"  Mean HAS hits: {results['jw18_uninf']['mean_has']:.2f}  Median: {results['jw18_uninf']['median_has']:.1f}")
        print(f"  Median distance: {results['jw18_uninf']['median_distance']:.0f} bp")
        print(f"\nJW18 wMel (down-regulated) interactions:")
        print(f"  N = {results['jw18_wmel']['n_interactions']}")
        print(f"  HAS overlap rate: {results['jw18_wmel']['overlap_rate']*100:.1f}%")
        print(f"  Mean HAS hits: {results['jw18_wmel']['mean_has']:.2f}  Median: {results['jw18_wmel']['median_has']:.1f}")
        print(f"  Median distance: {results['jw18_wmel']['median_distance']:.0f} bp")
        print(f"\nStatistical tests:")
        print(f"  Overlap rate difference (chi-square): p = {chi2_p:.2e}")
        print(f"  Distance difference (Mann-Whitney): p = {mw_p:.2e}")
        print(f"  HAS count difference (Mann-Whitney): p = {has_mw_p:.2e}")
        print(f"  HAS count difference (Welch's t-test): t = {t_stat:.3f}, p = {t_p:.2e}, Cohen's d = {cohens_d:.3f}")
        print(f"  One-sided binomial (wMel > uninf): k={k_wmel}/{n_wmel}, null p={p_null_binom:.3f}, p = {binom_p:.2e}")

    return results, merged


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def add_strip_to_boxplot(ax, data, positions, colors, jitter=0.08, alpha=0.5, size=8):
    """Add individual data points with jitter to a boxplot."""
    for d, pos, color in zip(data, positions, colors):
        if len(d) > 0:
            x = np.random.normal(pos, jitter, size=len(d))
            ax.scatter(x, d, alpha=alpha, s=size, color=color,
                       edgecolors='black', linewidths=0.3, zorder=3)


def create_visualization(merged_df, direction_results, distance_analysis, null_comparison, output_prefix):
    """Create two 2x2 figure sets. Set 1 now includes a boxplot+strip of HAS hits per contact."""
    print("\nCreating visualizations...")

    color_jw18_uninf = '#8fcb84'
    color_jw18_wmel  = '#09aa4b'
    colors = [color_jw18_uninf, color_jw18_wmel]
    labels = ['JW18 uninf.', 'JW18 wMel']

    plt.rcParams.update({'font.size': 6})

    # Per-genotype HAS count arrays for boxplots
    uninf_has = merged_df[merged_df['logFC'] > 0]['total_has'].values
    wmel_has  = merged_df[merged_df['logFC'] < 0]['total_has'].values

    # ------------------------------------------------------------------ #
    # Plot Set 1                                                           #
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)

    # --- (0,0): Observed vs Expected bar chart ---
    ax = axes[0, 0]
    ax.patch.set_alpha(0)
    rates = [null_comparison['real_overlap_rate'] * 100,
             null_comparison['null_expectation'] * 100]
    bar_colors = ['#cc3333' if null_comparison['p_value'] < 0.05 else '#999999', '#cccccc']
    bars = ax.bar(['Observed', 'Expected\n(Null)'], rates,
                  color=bar_colors, alpha=0.8, edgecolor='black', linewidth=1)
    ax.set_ylabel('HAS/CES Overlap Rate (%)', fontsize=6)
    ax.set_title(f'Enrichment vs Null\n({null_comparison["enrichment"]:.2f}x)',
                 fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    p_val = null_comparison['p_value']
    log10_p = null_comparison["log10_p_value"]
    if p_val > 0:
        p_str = f'p = {p_val:.2e}'
    elif np.isinf(log10_p):
        p_str = 'p << 0 (underflow)'
    else:
        p_str = f'p < 1e{int(np.floor(log10_p))}'
    ax.text(0.5, max(rates) * 0.95, p_str, ha='center', va='top', fontsize=5,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    for bar, val in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() * 0.5,
                f'{val:.1f}%', ha='center', va='center', fontsize=5)

    # --- (0,1): Direction-specific overlap rates with null ---
    ax = axes[0, 1]
    ax.patch.set_alpha(0)
    overlap_rates_dir = [
        direction_results['jw18_uninf']['overlap_rate'] * 100,
        direction_results['jw18_wmel']['overlap_rate'] * 100,
        null_comparison['null_expectation'] * 100
    ]
    bars = ax.bar(['JW18 uninf.', 'JW18 wMel', 'Expected\n(Null)'],
                  overlap_rates_dir,
                  color=[color_jw18_uninf, color_jw18_wmel, '#cccccc'],
                  alpha=0.8, edgecolor='black', linewidth=1)
    ax.set_ylabel('HAS Overlap Rate (%)', fontsize=6)
    ax.set_title('Genotype-Specific Enrichment', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    if 'comparison' in direction_results:
        pval = direction_results['comparison'].get('chi2_pvalue', 1.0)
        ax.text(0.5, max(overlap_rates_dir) * 0.95, f'p = {pval:.2e}',
                ha='center', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    for bar, val in zip(bars, overlap_rates_dir):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                f'{val:.1f}%', ha='center', va='bottom', fontsize=5)

    # --- (1,0): HAS hits per contact — BOXPLOT + STRIP ---
    # Note: total_has is integer (mostly 0/1/2) so we jitter y slightly for visibility
    ax = axes[1, 0]
    ax.patch.set_alpha(0)

    # Add y-jitter so overlapping integer points are visible
    np.random.seed(42)
    for arr, pos, color in zip([uninf_has, wmel_has], [0, 1], colors):
        if len(arr) > 0:
            x = np.random.normal(pos, 0.12, size=len(arr))
            y = arr + np.random.uniform(-0.15, 0.15, size=len(arr))
            ax.scatter(x, y, alpha=0.25, s=4, color=color,
                       edgecolors='none', zorder=2, rasterized=True)

    bp = ax.boxplot(
        [uninf_has, wmel_has],
        positions=[0, 1],
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(facecolor='white', edgecolor='black', linewidth=1.5),
        medianprops=dict(color='black', linewidth=2),
        whiskerprops=dict(color='black', linewidth=1),
        capprops=dict(color='black', linewidth=1)
    )
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylabel('Number of HAS/CES Sites per Contact', fontsize=6)
    ax.set_title('HAS/CES Hits per Contact', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    ax.set_yticks([0, 1, 2, 3])

    if len(uninf_has) > 0 and len(wmel_has) > 0:
        stat_lines = []
        try:
            _, mw_p = stats.mannwhitneyu(uninf_has, wmel_has, alternative='two-sided')
            stat_lines.append(f'Mann-Whitney p = {mw_p:.2e}')
        except Exception:
            pass
        try:
            _, t_p = stats.ttest_ind(uninf_has, wmel_has, equal_var=False, alternative='two-sided')
            stat_lines.append(f"Welch's t p = {t_p:.2e}")
        except Exception:
            pass
        # One-sided binomial: wMel overlap rate > uninf rate
        try:
            from scipy.stats import binomtest
            p_null_binom = (uninf_has > 0).mean()
            k = int((wmel_has > 0).sum())
            n = len(wmel_has)
            bp_result = binomtest(k, n, p_null_binom, alternative='greater')
            stat_lines.append(f'Binomial (wMel>uninf) p = {bp_result.pvalue:.2e}')
        except Exception:
            pass
        if stat_lines:
            ax.text(0.97, 0.97, '\n'.join(stat_lines),
                    transform=ax.transAxes, ha='right', va='top', fontsize=4.5,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='grey', linewidth=0.5))

    # Annotate medians
    for pos, arr, color in zip([0, 1], [uninf_has, wmel_has], colors):
        med = np.median(arr)
        mn  = np.mean(arr)
        ax.text(pos, -0.35, f'med={med:.1f}\nmean={mn:.2f}',
                ha='center', va='top', fontsize=4.5, color='black')

    # --- (1,1): Overlap type distribution ---
    ax = axes[1, 1]
    ax.patch.set_alpha(0)
    overlap_types = ['No overlap', 'One anchor', 'Both anchors']
    counts = [
        (~merged_df['any_anchor_overlap']).sum(),
        (merged_df['any_anchor_overlap'] & ~merged_df['both_anchors_overlap']).sum(),
        merged_df['both_anchors_overlap'].sum()
    ]
    bars = ax.bar(overlap_types, counts, color=['#cccccc', '#ff9933', '#cc3333'], alpha=0.8)
    ax.set_ylabel('Number of Interactions', fontsize=6)
    ax.set_title('HAS Overlap Pattern', fontsize=7, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    total = sum(counts)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                f'{count}\n({count/total*100:.1f}%)', ha='center', va='bottom', fontsize=5)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}/has_analysis_set1.pdf",
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"Plot set 1 saved to {output_prefix}/has_analysis_set1.pdf")

    # ------------------------------------------------------------------ #
    # Plot Set 2: Distance and effect size analysis                        #
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)

    # --- (0,0): LogFC vs distance to HAS ---
    ax = axes[0, 0]
    ax.patch.set_alpha(0)
    plot_data = merged_df[~np.isinf(merged_df['min_dist_any'])].copy()
    uninf_mask = plot_data['logFC'] > 0
    wmel_mask  = plot_data['logFC'] < 0
    ax.scatter(plot_data[uninf_mask]['logFC'], plot_data[uninf_mask]['min_dist_any'] / 1000,
               c=color_jw18_uninf, alpha=0.6, s=20, label='JW18 uninf.', edgecolors='none')
    ax.scatter(plot_data[wmel_mask]['logFC'],  plot_data[wmel_mask]['min_dist_any'] / 1000,
               c=color_jw18_wmel,  alpha=0.6, s=20, label='JW18 wMel',   edgecolors='none')
    ax.set_xlabel('log2 Fold Change', fontsize=6)
    ax.set_ylabel('Distance to Nearest HAS (kb)', fontsize=6)
    ax.set_title('Effect Size vs HAS Proximity', fontsize=7, fontweight='bold')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.set_yscale('log')
    ax.legend(fontsize=5)
    ax.tick_params(labelsize=6)
    corr, p_val = stats.spearmanr(plot_data['logFC'], plot_data['min_dist_any'])
    ax.text(0.05, 0.95, f'ρ = {corr:.3f}\np = {p_val:.2e}',
            transform=ax.transAxes, ha='left', va='top', fontsize=5,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))

    # --- (0,1): Distance distribution bar chart ---
    ax = axes[0, 1]
    ax.patch.set_alpha(0)
    if distance_analysis is not None:
        dist_data = distance_analysis['distance_distribution']
        bars = ax.bar(range(len(dist_data)), dist_data['percentage'],
                      color='steelblue', alpha=0.8)
        ax.set_xticks(range(len(dist_data)))
        ax.set_xticklabels(dist_data['distance_range'], rotation=45, ha='right', fontsize=6)
        ax.set_ylabel('Percentage of Interactions (%)', fontsize=6)
        ax.set_title('Distance to Nearest HAS', fontsize=7, fontweight='bold')
        ax.tick_params(labelsize=6)
        ax.text(0.95, 0.95, f'Median: {distance_analysis["median_distance"]/1000:.1f}kb',
                transform=ax.transAxes, ha='right', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
        for bar, val in zip(bars, dist_data['percentage']):
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=5)

    # --- (1,0): Median distance by genotype ---
    ax = axes[1, 0]
    ax.patch.set_alpha(0)
    median_dists = [
        direction_results['jw18_uninf']['median_distance'] / 1000,
        direction_results['jw18_wmel']['median_distance'] / 1000
    ]
    bars = ax.bar(labels, median_dists, color=colors, alpha=0.8)
    ax.set_ylabel('Median Distance to HAS (kb)', fontsize=6)
    ax.set_title('HAS Proximity by Genotype', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    if 'comparison' in direction_results and 'mannwhitney_pvalue' in direction_results['comparison']:
        mw_pval = direction_results['comparison']['mannwhitney_pvalue']
        ax.text(0.5, max(median_dists) * 0.95, f'p = {mw_pval:.2e}',
                ha='center', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    for bar, dist in zip(bars, median_dists):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                f'{dist:.1f}kb', ha='center', va='bottom', fontsize=5)

    # --- (1,1): Sample sizes ---
    ax = axes[1, 1]
    ax.patch.set_alpha(0)
    sample_sizes = [
        direction_results['jw18_uninf']['n_interactions'],
        direction_results['jw18_wmel']['n_interactions']
    ]
    bars = ax.bar(labels, sample_sizes, color=colors, alpha=0.8)
    ax.set_ylabel('Number of Interactions', fontsize=6)
    ax.set_title('Sample Sizes by Genotype', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    for bar, size in zip(bars, sample_sizes):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                f'n={size}', ha='center', va='bottom', fontsize=5)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}/has_analysis_set2.pdf",
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"Plot set 2 saved to {output_prefix}/has_analysis_set2.pdf")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze HAS/CES enrichment at differential X chromosome contacts'
    )
    parser.add_argument('--interactions', required=True,
                        help='Significant differential interactions CSV file')
    parser.add_argument('--null_model', required=True,
                        help='Null model interactions CSV file')
    parser.add_argument('--has_sites', required=True,
                        help='HAS/CES sites BED file')
    parser.add_argument('--window_size', type=int, default=50000,
                        help='Window size around HAS sites for MSL spreading (bp)')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                        help='FDR threshold for significant interactions')
    parser.add_argument('--output_prefix', required=True,
                        help='Output directory path')
    args = parser.parse_args()

    Path(args.output_prefix).mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("HAS/CES Enrichment Analysis at X Chromosome Contacts")
    print("="*60)

    has_sites = load_has_sites(args.has_sites)
    sig_interactions, all_x_interactions = load_differential_interactions(
        args.interactions, args.fdr_threshold)
    null_model = load_null_model(args.null_model)

    if has_sites is None or len(sig_interactions) == 0:
        print("Error: No data to analyze")
        return 1

    overlap_results = calculate_has_overlap(sig_interactions, has_sites, args.window_size)
    if overlap_results is None:
        print("Error: Could not calculate overlap")
        return 1

    direction_results, merged_df = analyze_by_logfc_direction(sig_interactions, overlap_results)
    distance_analysis = analyze_distance_distribution(merged_df, has_sites)

    overall_overlap_rate = overlap_results['any_anchor_overlap'].mean()
    n_overlaps     = int(overlap_results['any_anchor_overlap'].sum())
    n_interactions = len(overlap_results)
    null_comparison = compare_to_null_model(
        overall_overlap_rate, null_model, n_interactions, n_overlaps)

    print(f"\n{'='*60}")
    print("SUMMARY")
    print("="*60)
    print(f"  Overlap rate: {null_comparison['real_overlap_rate']*100:.1f}%")
    print(f"  Expected (null): {null_comparison['null_expectation']*100:.1f}%")
    print(f"  Enrichment: {null_comparison['enrichment']:.2f}x")
    p_val   = null_comparison['p_value']
    log10_p = null_comparison['log10_p_value']
    if p_val == 0.0:
        if np.isinf(log10_p):
            print(f"  P-value: vanishingly small (log10 p = -inf, extreme float underflow)")
        else:
            print(f"  P-value: < 1e{int(np.floor(log10_p))}")
    else:
        print(f"  P-value: {p_val:.2e}")
    print(f"  Significant: {'YES' if p_val < 0.05 else 'NO'}")

    create_visualization(merged_df, direction_results, distance_analysis,
                         null_comparison, args.output_prefix)

    merged_df.to_csv(f"{args.output_prefix}/has_interactions.csv", index=False)

    summary = {
        'total_x_interactions': len(all_x_interactions),
        'significant_x_interactions': len(sig_interactions),
        'has_sites_on_x': len(has_sites),
        'window_size': args.window_size,
        'overall_overlap_rate': overall_overlap_rate,
        'enrichment_vs_null': null_comparison['enrichment'],
        'p_value': null_comparison['p_value'],
        'log10_p_value': null_comparison['log10_p_value'],
        'jw18_uninf_n': direction_results['jw18_uninf']['n_interactions'],
        'jw18_uninf_overlap_rate': direction_results['jw18_uninf']['overlap_rate'],
        'jw18_uninf_mean_has': direction_results['jw18_uninf']['mean_has'],
        'jw18_uninf_median_has': direction_results['jw18_uninf']['median_has'],
        'jw18_uninf_median_distance': direction_results['jw18_uninf']['median_distance'],
        'jw18_wmel_n': direction_results['jw18_wmel']['n_interactions'],
        'jw18_wmel_overlap_rate': direction_results['jw18_wmel']['overlap_rate'],
        'jw18_wmel_mean_has': direction_results['jw18_wmel']['mean_has'],
        'jw18_wmel_median_has': direction_results['jw18_wmel']['median_has'],
        'jw18_wmel_median_distance': direction_results['jw18_wmel']['median_distance'],
        'direction_overlap_pvalue': direction_results.get('comparison', {}).get('chi2_pvalue'),
        'direction_distance_pvalue': direction_results.get('comparison', {}).get('mannwhitney_pvalue'),
        'has_count_mannwhitney_pvalue': direction_results.get('comparison', {}).get('has_count_mannwhitney_pvalue'),
        'has_count_effect_size_rb': direction_results.get('comparison', {}).get('has_count_effect_size'),
        'has_count_ttest_pvalue': direction_results.get('comparison', {}).get('has_count_ttest_pvalue'),
        'has_count_ttest_stat': direction_results.get('comparison', {}).get('has_count_ttest_stat'),
        'has_count_cohens_d': direction_results.get('comparison', {}).get('has_count_cohens_d'),
        'binom_pvalue_wmel_gt_uninf': direction_results.get('comparison', {}).get('binom_pvalue'),
        'binom_k': direction_results.get('comparison', {}).get('binom_k'),
        'binom_n': direction_results.get('comparison', {}).get('binom_n'),
        'binom_null_prop': direction_results.get('comparison', {}).get('binom_null_prop'),
        'median_distance_to_has': distance_analysis['median_distance'] if distance_analysis else None,
        'mean_distance_to_has': distance_analysis['mean_distance'] if distance_analysis else None
    }
    pd.DataFrame([summary]).to_csv(f"{args.output_prefix}/has_summary.csv", index=False)
    if distance_analysis is not None:
        distance_analysis['distance_distribution'].to_csv(
            f"{args.output_prefix}/distance_distribution.csv", index=False)

    with open(f"{args.output_prefix}/has_summary.txt", 'w') as f:
        f.write("HAS/CES ENRICHMENT ANALYSIS SUMMARY\n")
        f.write("="*60 + "\n\n")
        f.write(f"Analysis Parameters:\n")
        f.write(f"  Window size: {args.window_size} bp\n")
        f.write(f"  FDR threshold: {args.fdr_threshold}\n")
        f.write(f"  Note: Large window accounts for MSL spreading from HAS\n\n")
        f.write(f"Data Summary:\n")
        f.write(f"  Total X chromosome interactions: {summary['total_x_interactions']}\n")
        f.write(f"  Significant interactions: {summary['significant_x_interactions']}\n")
        f.write(f"  HAS/CES sites on X: {summary['has_sites_on_x']}\n\n")
        f.write(f"HAS/CES Enrichment:\n")
        f.write(f"  Overall overlap rate: {summary['overall_overlap_rate']*100:.1f}%\n")
        f.write(f"  Expected (null): {null_comparison['null_expectation']*100:.1f}%\n")
        f.write(f"  Enrichment: {summary['enrichment_vs_null']:.2f}x\n")
        p_val   = summary['p_value']
        log10_p = summary['log10_p_value']
        if p_val == 0.0:
            p_str = "vanishingly small (log10 p = -inf)" if np.isinf(log10_p) else f"< 1e{int(np.floor(log10_p))}"
        else:
            p_str = f"{p_val:.2e}"
        f.write(f"  P-value: {p_str}\n")
        f.write(f"  Significant: {'YES' if p_val < 0.05 else 'NO'}\n\n")
        if distance_analysis:
            f.write(f"Distance Analysis:\n")
            f.write(f"  Median distance to HAS: {summary['median_distance_to_has']/1000:.1f} kb\n")
            f.write(f"  Mean distance to HAS: {summary['mean_distance_to_has']/1000:.1f} kb\n\n")
        f.write(f"Genotype-specific Analysis:\n")
        for label, key in [('JW18 uninf. (up-regulated)', 'jw18_uninf'),
                            ('JW18 wMel (down-regulated)', 'jw18_wmel')]:
            r = direction_results[key]
            f.write(f"  {label}:\n")
            f.write(f"    N = {r['n_interactions']}\n")
            f.write(f"    HAS overlap: {r['overlap_rate']*100:.1f}%\n")
            f.write(f"    Mean HAS hits/contact: {r['mean_has']:.2f}\n")
            f.write(f"    Median HAS hits/contact: {r['median_has']:.1f}\n")
            f.write(f"    Median distance: {r['median_distance']/1000:.1f} kb\n")
        if summary['direction_overlap_pvalue'] is not None:
            f.write(f"\n  Statistical tests:\n")
            f.write(f"    Overlap rate (chi-square): p = {summary['direction_overlap_pvalue']:.2e}\n")
            f.write(f"    Distance (Mann-Whitney): p = {summary['direction_distance_pvalue']:.2e}\n")
            if summary['has_count_mannwhitney_pvalue'] is not None:
                f.write(f"    HAS hits/contact (Mann-Whitney): p = {summary['has_count_mannwhitney_pvalue']:.2e}\n")
            if summary['has_count_effect_size_rb'] is not None:
                f.write(f"    HAS hits effect size (rank-biserial): {summary['has_count_effect_size_rb']:.3f}\n")
            if summary['has_count_ttest_pvalue'] is not None:
                t_stat = summary['has_count_ttest_stat']
                t_p    = summary['has_count_ttest_pvalue']
                cd     = summary['has_count_cohens_d']
                f.write(f"    HAS hits/contact (Welch's t-test): t = {t_stat:.3f}, p = {t_p:.2e}\n")
                if cd is not None and not (isinstance(cd, float) and cd != cd):
                    f.write(f"    Cohen's d: {cd:.3f}\n")
            if summary['binom_pvalue_wmel_gt_uninf'] is not None:
                bk  = summary['binom_k']
                bn  = summary['binom_n']
                bp  = summary['binom_pvalue_wmel_gt_uninf']
                bnp = summary['binom_null_prop']
                f.write(f"    One-sided binomial (wMel overlap > uninf): k={bk}/{bn}, null p={bnp:.3f}, p = {bp:.2e}\n")
        f.write("\nBiological Interpretation:\n")
        f.write("-"*60 + "\n")
        enr = summary['enrichment_vs_null']
        pv  = summary['p_value']
        if enr > 1.5 and pv < 0.05:
            f.write("STRONG ENRICHMENT: Differential contacts are significantly enriched\n")
            f.write("near HAS/CES sites, suggesting Wolbachia affects chromatin contacts\n")
            f.write("at dosage compensation entry sites.\n")
        elif enr > 1.0 and pv < 0.05:
            f.write("MODERATE ENRICHMENT: Some enrichment detected near HAS sites.\n")
            f.write("Wolbachia may have subtle effects on dosage compensation regions.\n")
        else:
            f.write("NO SIGNIFICANT ENRICHMENT: Differential contacts are not enriched\n")
            f.write("near HAS sites. Changes may be independent of dosage compensation\n")
            f.write("entry sites, possibly affecting MSL spreading regions instead.\n")

    print(f"\nAnalysis complete! Results saved to {args.output_prefix}/")
    return 0


if __name__ == '__main__':
    exit(main())