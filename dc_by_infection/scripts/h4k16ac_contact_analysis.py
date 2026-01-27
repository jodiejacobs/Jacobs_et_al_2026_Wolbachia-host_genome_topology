#!/usr/bin/env python3
"""
Analyze H4K16ac enrichment at differential chromatin contacts on the X chromosome.
Compares real differential interactions vs permutation-based null model to assess significance.

OPTIMIZED VERSION with parallel processing and vectorized operations.

Based on:
- Fei et al. - NDF associates with MSL complex via H4K16ac
- Larschan et al. - H4K16ac marks active genes for dosage compensation
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
import gzip
import sys
from multiprocessing import Pool, cpu_count
from functools import partial

def load_h4k16ac_peaks(chip_file):
    """
    Load H4K16ac ChIP-seq peaks from GSE20784 or similar datasets.
    Expects BED format: chr, start, end, score
    Format: chr2L 67843 71268 1.90337415112343
    """
    print(f"Loading H4K16ac peaks from {chip_file}")
    
    try:
        # Read the file (handles both compressed and uncompressed)
        if chip_file.endswith('.gz'):
            import gzip
            with gzip.open(chip_file, 'rt') as f:
                peaks = []
                for line in f:
                    if line.startswith('#') or line.startswith('track'):
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        # Normalize chromosome name (remove 'chr' prefix if present)
                        chrom = parts[0].replace('chr', '')
                        peaks.append({
                            'chrom': chrom,
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'score': float(parts[3]) if len(parts) > 3 else 0
                        })
        else:
            # Read uncompressed file
            with open(chip_file, 'r') as f:
                peaks = []
                for line in f:
                    if line.startswith('#') or line.startswith('track'):
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        # Normalize chromosome name (remove 'chr' prefix if present)
                        chrom = parts[0].replace('chr', '')
                        peaks.append({
                            'chrom': chrom,
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'score': float(parts[3]) if len(parts) > 3 else 0
                        })
        
        peaks_df = pd.DataFrame(peaks)
        
        print(f"Total peaks loaded: {len(peaks_df)}")
        print(f"Chromosomes present: {sorted(peaks_df['chrom'].unique())}")
        
        # Filter for X chromosome (handle both 'X' and 'chrX')
        x_peaks = peaks_df[peaks_df['chrom'] == 'X'].copy()
        
        if len(x_peaks) == 0:
            print("Warning: No peaks found on X chromosome")
            print(f"Available chromosomes: {peaks_df['chrom'].unique()}")
            return None
        
        print(f"Loaded {len(x_peaks)} H4K16ac peaks on X chromosome")
        print(f"Score range: {x_peaks['score'].min():.2f} - {x_peaks['score'].max():.2f}")
        print(f"Position range: {x_peaks['start'].min():,} - {x_peaks['end'].max():,}")
        
        return x_peaks
        
    except Exception as e:
        print(f"Error loading H4K16ac peaks: {e}")
        import traceback
        traceback.print_exc()
        return None

def load_differential_interactions(interactions_file, fdr_threshold=0.05):
    """Load significant differential interactions"""
    print(f"Loading differential interactions from {interactions_file}")
    
    interactions = pd.read_csv(interactions_file)
    
    # Normalize chromosome names (remove 'chr' prefix if present)
    if 'chr1' in interactions.columns:
        interactions['chr1'] = interactions['chr1'].astype(str).str.replace('chr', '')
        interactions['chr2'] = interactions['chr2'].astype(str).str.replace('chr', '')
    
    print(f"Loaded {len(interactions)} total interactions")
    print(f"Chromosomes in data: {sorted(interactions['chr1'].unique())}")
    
    # Filter for X chromosome interactions
    x_interactions = interactions[
        ((interactions['chr1'] == 'X') | (interactions['chr2'] == 'X'))
    ].copy()
    
    print(f"Found {len(x_interactions)} total X chromosome interactions")
    
    # Filter for significant
    sig_x = x_interactions[
        (x_interactions['FDR'] < fdr_threshold) & 
        (abs(x_interactions['logFC']) > 1)
    ].copy()
    
    print(f"Found {len(sig_x)} significant X chromosome interactions (FDR < {fdr_threshold}, |logFC| > 1)")
    
    if len(sig_x) > 0:
        print(f"  Cis interactions: {sum(sig_x['chr1'] == sig_x['chr2'])}")
        print(f"  Trans interactions: {sum(sig_x['chr1'] != sig_x['chr2'])}")
        print(f"  JW18 uninf. (up-regulated): {sum(sig_x['logFC'] > 0)}")
        print(f"  JW18 wMel (down-regulated): {sum(sig_x['logFC'] < 0)}")
        print(f"  Position range anchor1: {sig_x['start1'].min():,} - {sig_x['end1'].max():,}")
        print(f"  Position range anchor2: {sig_x['start2'].min():,} - {sig_x['end2'].max():,}")
    
    return sig_x, x_interactions

def calculate_h4k16ac_overlap_vectorized(interactions_df, h4k16ac_peaks, window_size=5000):
    """
    OPTIMIZED: Vectorized overlap calculation using single BedTools intersect call.
    Much faster than looping through interactions.
    """
    if h4k16ac_peaks is None or len(h4k16ac_peaks) == 0:
        return None
    
    print(f"  Calculating overlaps for {len(interactions_df)} interactions...")
    
    # Create extended anchor regions for all interactions at once
    anchor1_regions = []
    anchor2_regions = []
    
    for idx, row in interactions_df.iterrows():
        anchor1_regions.append({
            'chrom': row['chr1'],
            'start': max(0, row['start1'] - window_size),
            'end': row['end1'] + window_size,
            'idx': idx
        })
        anchor2_regions.append({
            'chrom': row['chr2'],
            'start': max(0, row['start2'] - window_size),
            'end': row['end2'] + window_size,
            'idx': idx
        })
    
    anchor1_df = pd.DataFrame(anchor1_regions)
    anchor2_df = pd.DataFrame(anchor2_regions)

    # Create BedTools objects
    anchor1_bt = pybedtools.BedTool.from_dataframe(anchor1_df[['chrom', 'start', 'end', 'idx']])
    anchor2_bt = pybedtools.BedTool.from_dataframe(anchor2_df[['chrom', 'start', 'end', 'idx']])
    h4k16_bt = pybedtools.BedTool.from_dataframe(h4k16ac_peaks[['chrom', 'start', 'end']])
    
    # Intersect all at once
    anchor1_overlaps = anchor1_bt.intersect(h4k16_bt, wa=True, c=True)
    anchor2_overlaps = anchor2_bt.intersect(h4k16_bt, wa=True, c=True)
    
    # Parse results
    anchor1_counts = {}
    for interval in anchor1_overlaps:
        idx = int(interval.fields[3])
        count = int(interval.fields[4])
        anchor1_counts[idx] = count
    
    anchor2_counts = {}
    for interval in anchor2_overlaps:
        idx = int(interval.fields[3])
        count = int(interval.fields[4])
        anchor2_counts[idx] = count
    
    # Build results dataframe
    results = []
    for idx in interactions_df.index:
        n_peaks_1 = anchor1_counts.get(idx, 0)
        n_peaks_2 = anchor2_counts.get(idx, 0)
        
        results.append({
            'interaction_idx': idx,
            'anchor1_overlap': n_peaks_1 > 0,
            'anchor2_overlap': n_peaks_2 > 0,
            'both_anchors_overlap': (n_peaks_1 > 0) and (n_peaks_2 > 0),
            'any_anchor_overlap': (n_peaks_1 > 0) or (n_peaks_2 > 0),
            'n_peaks_anchor1': n_peaks_1,
            'n_peaks_anchor2': n_peaks_2,
            'total_peaks': n_peaks_1 + n_peaks_2
        })
    
    return pd.DataFrame(results)

def permute_once(seed, interactions_df, h4k16ac_peaks, window_size, x_chrom_end):
    """
    Single permutation - designed to be called in parallel.
    """
    np.random.seed(seed)
    
    # Shuffle anchor positions - maintain anchor sizes
    permuted = interactions_df.copy()
    
    # Determine which anchors are on X
    anchor1_on_x = permuted['chr1'] == 'X'
    anchor2_on_x = permuted['chr2'] == 'X'
    
    # Permute anchor 1 where it's on X
    if anchor1_on_x.any():
        anchor1_sizes = permuted.loc[anchor1_on_x, 'end1'] - permuted.loc[anchor1_on_x, 'start1']
        new_starts1 = np.random.randint(0, x_chrom_end - anchor1_sizes.max() + 1, size=anchor1_on_x.sum())
        # Ensure no overflow
        new_starts1 = np.minimum(new_starts1, x_chrom_end - anchor1_sizes.values)
        permuted.loc[anchor1_on_x, 'start1'] = new_starts1
        permuted.loc[anchor1_on_x, 'end1'] = new_starts1 + anchor1_sizes.values
    
    # Permute anchor 2 where it's on X
    if anchor2_on_x.any():
        anchor2_sizes = permuted.loc[anchor2_on_x, 'end2'] - permuted.loc[anchor2_on_x, 'start2']
        new_starts2 = np.random.randint(0, x_chrom_end - anchor2_sizes.max() + 1, size=anchor2_on_x.sum())
        # Ensure no overflow
        new_starts2 = np.minimum(new_starts2, x_chrom_end - anchor2_sizes.values)
        permuted.loc[anchor2_on_x, 'start2'] = new_starts2
        permuted.loc[anchor2_on_x, 'end2'] = new_starts2 + anchor2_sizes.values
    
    # Calculate overlap for this permutation
    perm_overlap = calculate_h4k16ac_overlap_vectorized(permuted, h4k16ac_peaks, window_size)
    
    if perm_overlap is not None:
        return perm_overlap['any_anchor_overlap'].mean()
    else:
        return np.nan

def permutation_test_h4k16ac_enrichment(interactions_df, h4k16ac_peaks, 
                                          overlap_results, window_size=5000, 
                                          n_permutations=1000, n_cores=None):
    """
    OPTIMIZED: Parallel permutation test with vectorized overlap calculation.
    """
    print(f"\nPerforming permutation test ({n_permutations} permutations)...")
    
    # Get observed overlap rate
    observed_rate = overlap_results['any_anchor_overlap'].mean()
    
    # Get X chromosome size
    x_chrom_end = max(
        interactions_df[interactions_df['chr1'] == 'X']['end1'].max() if any(interactions_df['chr1'] == 'X') else 0,
        interactions_df[interactions_df['chr2'] == 'X']['end2'].max() if any(interactions_df['chr2'] == 'X') else 0
    )
    
    print(f"X chromosome size for permutations: {x_chrom_end:,} bp")
    print(f"Observed overlap rate: {observed_rate*100:.1f}%")
    
    # Identify interaction types
    cis_x = (interactions_df['chr1'] == 'X') & (interactions_df['chr2'] == 'X')
    trans_x_first = (interactions_df['chr1'] == 'X') & (interactions_df['chr2'] != 'X')
    trans_x_second = (interactions_df['chr1'] != 'X') & (interactions_df['chr2'] == 'X')
    
    print(f"  Cis-X interactions: {cis_x.sum()}")
    print(f"  Trans interactions (X as anchor1): {trans_x_first.sum()}")
    print(f"  Trans interactions (X as anchor2): {trans_x_second.sum()}")
    
    # Determine number of cores
    if n_cores is None:
        n_cores = max(1, cpu_count() - 1)  # Leave one core free
    print(f"  Using {n_cores} cores for parallel processing")
    
    # Create seeds for each permutation (for reproducibility)
    seeds = [42 + i for i in range(n_permutations)]
    
    # Create partial function with fixed arguments
    permute_func = partial(
        permute_once,
        interactions_df=interactions_df,
        h4k16ac_peaks=h4k16ac_peaks,
        window_size=window_size,
        x_chrom_end=x_chrom_end
    )
    
    # Run permutations in parallel
    print(f"  Running {n_permutations} permutations in parallel...")
    with Pool(processes=n_cores) as pool:
        null_overlap_rates = pool.map(permute_func, seeds)
    
    # Remove any NaN values
    null_overlap_rates = np.array([x for x in null_overlap_rates if not np.isnan(x)])
    
    print(f"  Completed {len(null_overlap_rates)} permutations successfully")
    
    # Calculate empirical p-value (two-tailed)
    p_value_right = np.mean(null_overlap_rates >= observed_rate)
    p_value_left = np.mean(null_overlap_rates <= observed_rate)
    p_value = 2 * min(p_value_right, p_value_left)
    
    # Calculate z-score
    null_mean = null_overlap_rates.mean()
    null_std = null_overlap_rates.std()
    z_score = (observed_rate - null_mean) / null_std if null_std > 0 else 0
    
    print(f"\nPermutation Test Results:")
    print(f"  Observed overlap rate: {observed_rate*100:.1f}%")
    print(f"  Null mean: {null_mean*100:.1f}%")
    print(f"  Null std: {null_std*100:.2f}%")
    print(f"  Z-score: {z_score:.2f}")
    print(f"  Empirical p-value: {p_value:.4f}")
    print(f"  Significant: {'YES' if p_value < 0.05 else 'NO'}")
    
    return {
        'observed_rate': observed_rate,
        'null_distribution': null_overlap_rates,
        'null_mean': null_mean,
        'null_std': null_std,
        'z_score': z_score,
        'p_value': p_value
    }

def ensure_boolean_columns(df, columns):
    """
    Ensure specified columns are boolean type and handle NaN values.
    This prevents TypeError when using the ~ operator.
    """
    df = df.copy()
    for col in columns:
        if col in df.columns:
            df[col] = df[col].fillna(False).astype(bool)
    return df

def analyze_by_logfc_direction(interactions_df, overlap_results):
    """Analyze H4K16ac enrichment by logFC direction (JW18 uninf. vs JW18 wMel)"""
    print("\nAnalyzing by logFC direction...")
    
    # FIXED: Proper merge with reset indices
    interactions_clean = interactions_df.reset_index(drop=True)
    overlap_clean = overlap_results.reset_index(drop=True)
    
    # Verify alignment
    assert len(interactions_clean) == len(overlap_clean), \
        f"Mismatch: {len(interactions_clean)} interactions vs {len(overlap_clean)} overlap results"
    
    # Merge by copying columns
    merged = interactions_clean.copy()
    for col in overlap_clean.columns:
        if col != 'interaction_idx':
            merged[col] = overlap_clean[col].values
    
    # Ensure boolean columns are properly typed
    bool_cols = ['anchor1_overlap', 'anchor2_overlap', 'both_anchors_overlap', 'any_anchor_overlap']
    merged = ensure_boolean_columns(merged, bool_cols)
    
    # Add interaction type
    merged['interaction_type'] = np.where(merged['chr1'] == merged['chr2'], 'cis', 'trans')
    
    # Split by direction
    jw18_uninf = merged[merged['logFC'] > 0]
    jw18_wmel = merged[merged['logFC'] < 0]
    
    # Calculate statistics
    results = {
        'jw18_uninf': {
            'n_interactions': len(jw18_uninf),
            'overlap_rate': jw18_uninf['any_anchor_overlap'].mean(),
            'both_anchors_rate': jw18_uninf['both_anchors_overlap'].mean(),
            'mean_peaks': jw18_uninf['total_peaks'].mean(),
            'median_peaks': jw18_uninf['total_peaks'].median()
        },
        'jw18_wmel': {
            'n_interactions': len(jw18_wmel),
            'overlap_rate': jw18_wmel['any_anchor_overlap'].mean(),
            'both_anchors_rate': jw18_wmel['both_anchors_overlap'].mean(),
            'mean_peaks': jw18_wmel['total_peaks'].mean(),
            'median_peaks': jw18_wmel['total_peaks'].median()
        }
    }
    
    # Statistical comparison 1: Fisher's exact for overlap rate
    contingency = [
        [jw18_uninf['any_anchor_overlap'].sum(), 
         (~jw18_uninf['any_anchor_overlap']).sum()],
        [jw18_wmel['any_anchor_overlap'].sum(), 
         (~jw18_wmel['any_anchor_overlap']).sum()]
    ]
    
    odds_ratio, fisher_p = stats.fisher_exact(contingency)
    
    # Statistical comparison 2: Mann-Whitney U for peak counts
    if len(jw18_uninf) > 0 and len(jw18_wmel) > 0:
        try:
            u_stat, mw_p = stats.mannwhitneyu(
                jw18_uninf['total_peaks'], 
                jw18_wmel['total_peaks'], 
                alternative='two-sided'
            )
            
            # Effect size
            n1 = len(jw18_uninf)
            n2 = len(jw18_wmel)
            rank_biserial = 1 - (2*u_stat) / (n1 * n2)
        except Exception as e:
            print(f"Warning: Could not compute Mann-Whitney U test: {e}")
            mw_p = np.nan
            rank_biserial = np.nan
    else:
        mw_p = np.nan
        rank_biserial = np.nan
    
    results['comparison'] = {
        'fisher_odds_ratio': odds_ratio,
        'fisher_p_value': fisher_p,
        'mannwhitneyu_p_value': mw_p,
        'effect_size_rank_biserial': rank_biserial
    }
    
    print(f"\nDirection-specific results:")
    print(f"JW18 uninf. (n={results['jw18_uninf']['n_interactions']}):")
    print(f"  H4K16ac overlap rate: {results['jw18_uninf']['overlap_rate']*100:.1f}%")
    print(f"  Mean peaks: {results['jw18_uninf']['mean_peaks']:.2f}")
    print(f"  Median peaks: {results['jw18_uninf']['median_peaks']:.1f}")
    print(f"JW18 wMel (n={results['jw18_wmel']['n_interactions']}):")
    print(f"  H4K16ac overlap rate: {results['jw18_wmel']['overlap_rate']*100:.1f}%")
    print(f"  Mean peaks: {results['jw18_wmel']['mean_peaks']:.2f}")
    print(f"  Median peaks: {results['jw18_wmel']['median_peaks']:.1f}")
    print(f"\nStatistical comparisons:")
    print(f"  Overlap rate (Fisher's exact): p = {fisher_p:.2e}")
    print(f"  Peak counts (Mann-Whitney U): p = {mw_p:.2e}")
    if not np.isnan(rank_biserial):
        print(f"  Effect size (rank-biserial): {rank_biserial:.3f}")
    
    return results, merged

def create_permutation_visualization(perm_results, output_prefix):
    """Create visualization of permutation test results"""
    print("\nCreating permutation test visualization...")
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    
    # Plot null distribution
    ax.hist(perm_results['null_distribution'] * 100, bins=30, 
            color='gray', alpha=0.7, edgecolor='black', linewidth=0.5,
            label='Permuted null')
    
    # Add observed value
    ax.axvline(perm_results['observed_rate'] * 100, 
               color='red', linestyle='--', linewidth=2,
               label=f"Observed ({perm_results['observed_rate']*100:.1f}%)")
    
    # Add null mean
    ax.axvline(perm_results['null_mean'] * 100,
               color='blue', linestyle=':', linewidth=2,
               label=f"Null mean ({perm_results['null_mean']*100:.1f}%)")
    
    ax.set_xlabel('H4K16ac Overlap Rate (%)', fontsize=8)
    ax.set_ylabel('Frequency', fontsize=8)
    ax.set_title('Permutation Test: H4K16ac Enrichment', fontsize=10, fontweight='bold')
    ax.legend(fontsize=7)
    ax.tick_params(labelsize=7)
    
    # Add statistics text box
    stats_text = (f"Z-score: {perm_results['z_score']:.2f}\n"
                  f"P-value: {perm_results['p_value']:.4f}")
    ax.text(0.05, 0.95, stats_text,
            transform=ax.transAxes, ha='left', va='top', fontsize=7,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_permutation_test.pdf",
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Permutation plot saved to {output_prefix}_permutation_test.pdf")

def add_strip_to_boxplot(ax, data, positions, colors, jitter=0.08, alpha=0.5, size=8):
    """
    Add individual data points to a boxplot with jitter.
    
    Parameters:
    -----------
    ax : matplotlib axis
    data : list of arrays
        Data for each boxplot
    positions : list
        X positions for boxplots
    colors : list
        Colors for each group
    jitter : float
        Amount of horizontal jitter
    alpha : float
        Point transparency
    size : float
        Point size
    """
    for i, (d, pos, color) in enumerate(zip(data, positions, colors)):
        if len(d) > 0:  # Only plot if data exists
            # Add jitter to x positions
            x = np.random.normal(pos, jitter, size=len(d))
            ax.scatter(x, d, alpha=alpha, s=size, color=color, edgecolors='black', 
                      linewidths=0.3, zorder=3)

def create_visualization(merged_df, direction_results, peaks_per_int_results, output_prefix):
    """Create separate 2x2 visualizations with boxplots and individual points"""
    print("\nCreating visualizations...")
    
    # Ensure all boolean columns are properly typed before visualization
    bool_cols = ['anchor1_overlap', 'anchor2_overlap', 'both_anchors_overlap', 'any_anchor_overlap']
    merged_df = ensure_boolean_columns(merged_df, bool_cols)
    
    # Define colors for JW18 uninf. and JW18 wMel
    color_jw18_uninf = '#8fcb84'  # Light green for JW18 uninf. (upregulated)
    color_jw18_wmel = '#09aa4b'   # Dark green for JW18 wMel (downregulated)
    colors = [color_jw18_uninf, color_jw18_wmel]
    
    # Define labels
    labels = ['JW18 uninf.', 'JW18 wMel']
    
    # Set global font size to minimum 6pt
    plt.rcParams.update({'font.size': 6})
    
    # Prepare data for boxplots
    uninf_peaks = merged_df[merged_df['logFC'] > 0]['total_peaks'].values
    wmel_peaks = merged_df[merged_df['logFC'] < 0]['total_peaks'].values
    
    print(f"  JW18 uninf peaks: n={len(uninf_peaks)}, range={uninf_peaks.min()}-{uninf_peaks.max()}")
    print(f"  JW18 wMel peaks: n={len(wmel_peaks)}, range={wmel_peaks.min()}-{wmel_peaks.max()}")
    
    # --- Plot Set 1: Overlap rates and peak counts (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)
    
    # Plot 1: H4K16ac Peak Count by Genotype (BOXPLOT)
    ax = axes[0, 0]
    ax.patch.set_alpha(0)
    
    bp = ax.boxplot([uninf_peaks, wmel_peaks], 
                     positions=[0, 1],
                     widths=0.6,
                     patch_artist=True,
                     showfliers=False,  # We'll add points manually
                     boxprops=dict(facecolor='white', edgecolor='black', linewidth=1),
                     medianprops=dict(color='black', linewidth=1.5),
                     whiskerprops=dict(color='black', linewidth=1),
                     capprops=dict(color='black', linewidth=1))
    
    # Color the boxes
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    
    # Add individual points with jitter
    add_strip_to_boxplot(ax, [uninf_peaks, wmel_peaks], [0, 1], colors, 
                         jitter=0.08, alpha=0.5, size=8)
    
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylabel('Number of H4K16ac Peaks', fontsize=6)
    ax.set_title('H4K16ac Peak Count by Genotype', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add Mann-Whitney U test
    if len(uninf_peaks) > 0 and len(wmel_peaks) > 0:
        try:
            u_stat, p_val = stats.mannwhitneyu(uninf_peaks, wmel_peaks, alternative='two-sided')
            ax.text(0.95, 0.95, f'p = {p_val:.2e}', 
                    transform=ax.transAxes, ha='right', va='top', fontsize=5,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
        except Exception as e:
            print(f"Warning: Could not compute Mann-Whitney U test: {e}")
    
    # Plot 2: NEW - Peaks per Interaction (Bar plot with error bars)
    ax = axes[0, 1]
    ax.patch.set_alpha(0)
    
    peaks_per_int = [
        peaks_per_int_results['jw18_uninf']['peaks_per_interaction'],
        peaks_per_int_results['jw18_wmel']['peaks_per_interaction']
    ]
    
    # Calculate standard errors for error bars
    se_uninf = peaks_per_int_results['jw18_uninf']['std_peaks'] / np.sqrt(peaks_per_int_results['jw18_uninf']['n_interactions'])
    se_wmel = peaks_per_int_results['jw18_wmel']['std_peaks'] / np.sqrt(peaks_per_int_results['jw18_wmel']['n_interactions'])
    errors = [se_uninf, se_wmel]
    
    bars = ax.bar([0, 1], peaks_per_int, yerr=errors, capsize=5,
                   color=colors, alpha=0.8, width=0.6, 
                   error_kw={'linewidth': 1, 'ecolor': 'black'})
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylabel('H4K16ac Peaks per Interaction', fontsize=6)
    ax.set_title('Normalized H4K16ac Enrichment', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add p-value
    if 'comparison' in peaks_per_int_results and 'mannwhitneyu_p_value' in peaks_per_int_results['comparison']:
        p_val = peaks_per_int_results['comparison']['mannwhitneyu_p_value']
        if not np.isnan(p_val):
            y_pos = max(peaks_per_int[0] + errors[0], peaks_per_int[1] + errors[1]) * 1.15
            ax.text(0.5, y_pos, f'p = {p_val:.2e}', 
                    ha='center', va='bottom', fontsize=5, 
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add value labels on bars
    for i, (bar, rate) in enumerate(zip(bars, peaks_per_int)):
        height = bar.get_height()
        ax.text(i, height + errors[i],
                f'{rate:.2f}',
                ha='center', va='bottom', fontsize=5)
    
    # Plot 3: LogFC vs H4K16ac peaks (Scatter)
    ax = axes[1, 0]
    ax.patch.set_alpha(0)
    
    # Color points by genotype
    uninf_mask = merged_df['logFC'] > 0
    wmel_mask = merged_df['logFC'] < 0
    
    ax.scatter(merged_df[uninf_mask]['logFC'], merged_df[uninf_mask]['total_peaks'], 
               c=color_jw18_uninf, alpha=0.5, s=15, label='JW18 uninf.', edgecolors='none')
    ax.scatter(merged_df[wmel_mask]['logFC'], merged_df[wmel_mask]['total_peaks'], 
               c=color_jw18_wmel, alpha=0.5, s=15, label='JW18 wMel', edgecolors='none')
    
    ax.set_xlabel('log2 Fold Change', fontsize=6)
    ax.set_ylabel('Number of H4K16ac Peaks', fontsize=6)
    ax.set_title('H4K16ac Enrichment vs Effect Size', fontsize=7, fontweight='bold')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.legend(fontsize=5)
    ax.tick_params(labelsize=6)
    
    # Add correlation statistics
    try:
        corr, p_val = stats.spearmanr(merged_df['logFC'], merged_df['total_peaks'])
        ax.text(0.05, 0.95, f'ρ = {corr:.3f}\np = {p_val:.2e}', 
                transform=ax.transAxes, ha='left', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    except Exception as e:
        print(f"Warning: Could not compute Spearman correlation: {e}")
    
    # Plot 4: Cis vs Trans interactions (BOXPLOT)
    ax = axes[1, 1]
    ax.patch.set_alpha(0)
    
    if 'interaction_type' in merged_df.columns:
        cis_peaks = merged_df[merged_df['interaction_type'] == 'cis']['total_peaks'].values
        trans_peaks = merged_df[merged_df['interaction_type'] == 'trans']['total_peaks'].values
        
        if len(cis_peaks) > 0 and len(trans_peaks) > 0:
            bp = ax.boxplot([cis_peaks, trans_peaks], 
                             positions=[0, 1],
                             widths=0.6,
                             patch_artist=True,
                             showfliers=False,
                             boxprops=dict(facecolor='white', edgecolor='black', linewidth=1),
                             medianprops=dict(color='black', linewidth=1.5),
                             whiskerprops=dict(color='black', linewidth=1),
                             capprops=dict(color='black', linewidth=1))
            
            # Color the boxes
            box_colors = ['#3366cc', '#33cc66']
            for patch, color in zip(bp['boxes'], box_colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.6)
            
            # Add individual points
            add_strip_to_boxplot(ax, [cis_peaks, trans_peaks], [0, 1], box_colors,
                                 jitter=0.08, alpha=0.5, size=8)
            
            ax.set_xticks([0, 1])
            ax.set_xticklabels(['Cis', 'Trans'])
            ax.set_ylabel('Number of H4K16ac Peaks', fontsize=6)
            ax.set_title('H4K16ac by Interaction Type', fontsize=7, fontweight='bold')
            ax.tick_params(labelsize=6)
            
            # Add Mann-Whitney U test
            try:
                u_stat, p_val = stats.mannwhitneyu(cis_peaks, trans_peaks, alternative='two-sided')
                ax.text(0.95, 0.95, f'p = {p_val:.2e}', 
                        transform=ax.transAxes, ha='right', va='top', fontsize=5,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
            except Exception as e:
                print(f"Warning: Could not compute Mann-Whitney U test for cis/trans: {e}")
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_analysis_set1.pdf", 
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Plot set 1 saved to {output_prefix}_analysis_set1.pdf")
    
    # --- Plot Set 2: Additional analyses (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)
    
    # Plot 1: H4K16ac Overlap Rate (Bar plot)
    ax = axes[0, 0]
    ax.patch.set_alpha(0)
    overlap_rates = [
        direction_results['jw18_uninf']['overlap_rate'] * 100,
        direction_results['jw18_wmel']['overlap_rate'] * 100
    ]
    bars = ax.bar([0, 1], overlap_rates, color=colors, alpha=0.8, width=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylabel('H4K16ac Overlap Rate (%)', fontsize=6)
    ax.set_title('H4K16ac Enrichment by Genotype', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add p-value if available
    if 'comparison' in direction_results and 'fisher_p_value' in direction_results['comparison']:
        p_val = direction_results['comparison']['fisher_p_value']
        ax.text(0.5, max(overlap_rates) * 0.95, f'p = {p_val:.2e}', 
                ha='center', va='top', fontsize=5, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add value labels on bars
    for i, (bar, rate) in enumerate(zip(bars, overlap_rates)):
        height = bar.get_height()
        ax.text(i, height,
                f'{rate:.1f}%',
                ha='center', va='bottom', fontsize=5)
    
    # Plot 2: Anchor-specific enrichment (BOXPLOT)
    ax = axes[0, 1]
    ax.patch.set_alpha(0)
    
    uninf_a1 = merged_df[merged_df['logFC'] > 0]['n_peaks_anchor1'].values
    uninf_a2 = merged_df[merged_df['logFC'] > 0]['n_peaks_anchor2'].values
    wmel_a1 = merged_df[merged_df['logFC'] < 0]['n_peaks_anchor1'].values
    wmel_a2 = merged_df[merged_df['logFC'] < 0]['n_peaks_anchor2'].values
    
    if all(len(arr) > 0 for arr in [uninf_a1, uninf_a2, wmel_a1, wmel_a2]):
        bp = ax.boxplot([uninf_a1, uninf_a2, wmel_a1, wmel_a2],
                         positions=[0, 1, 2.5, 3.5],
                         widths=0.6,
                         patch_artist=True,
                         showfliers=False,
                         boxprops=dict(facecolor='white', edgecolor='black', linewidth=1),
                         medianprops=dict(color='black', linewidth=1.5),
                         whiskerprops=dict(color='black', linewidth=1),
                         capprops=dict(color='black', linewidth=1))
        
        # Color the boxes
        anchor_colors = [color_jw18_uninf, color_jw18_uninf, color_jw18_wmel, color_jw18_wmel]
        for patch, color in zip(bp['boxes'], anchor_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        
        # Add individual points
        add_strip_to_boxplot(ax, [uninf_a1, uninf_a2, wmel_a1, wmel_a2], 
                             [0, 1, 2.5, 3.5], anchor_colors,
                             jitter=0.08, alpha=0.5, size=8)
        
        ax.set_xticks([0.5, 3])
        ax.set_xticklabels(['JW18 uninf.', 'JW18 wMel'])
        ax.set_ylabel('Number of H4K16ac Peaks', fontsize=6)
        ax.set_title('H4K16ac Peaks per Anchor', fontsize=7, fontweight='bold')
        ax.tick_params(labelsize=6)
    
    # Plot 3: Sample size information
    ax = axes[1, 0]
    ax.patch.set_alpha(0)
    
    sample_sizes = [
        direction_results['jw18_uninf']['n_interactions'],
        direction_results['jw18_wmel']['n_interactions']
    ]
    bars = ax.bar([0, 1], sample_sizes, color=colors, alpha=0.8, width=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylabel('Number of Interactions', fontsize=6)
    ax.set_title('Sample Sizes by Genotype', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add value labels on bars
    for i, (bar, size) in enumerate(zip(bars, sample_sizes)):
        height = bar.get_height()
        ax.text(i, height,
                f'n={size}',
                ha='center', va='bottom', fontsize=5)
    
    # Plot 4: Overlap type distribution
    ax = axes[1, 1]
    ax.patch.set_alpha(0)
    overlap_types = ['No overlap', 'One anchor', 'Both anchors']
    counts = [
        (~merged_df['any_anchor_overlap']).sum(),
        (merged_df['any_anchor_overlap'] & ~merged_df['both_anchors_overlap']).sum(),
        merged_df['both_anchors_overlap'].sum()
    ]
    bars = ax.bar([0, 1, 2], counts, color=['#cccccc', '#ff9933', '#cc3333'], alpha=0.8, width=0.6)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(overlap_types)
    ax.set_ylabel('Number of Interactions', fontsize=6)
    ax.set_title('H4K16ac Overlap Pattern', fontsize=7, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    
    # Add value labels on bars
    for i, (bar, count) in enumerate(zip(bars, counts)):
        height = bar.get_height()
        ax.text(i, height,
                f'{count}',
                ha='center', va='bottom', fontsize=5)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_analysis_set2.pdf", 
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Plot set 2 saved to {output_prefix}_analysis_set2.pdf")

def analyze_peaks_per_interaction(interactions_df, overlap_results, h4k16ac_peaks):
    """
    NEW: Calculate peaks per interaction for wMel vs uninf comparison.
    This normalizes by the number of interactions.
    """
    print("\nAnalyzing peaks per interaction (wMel vs uninf)...")
    
    # Merge interactions with overlap results
    interactions_clean = interactions_df.reset_index(drop=True)
    overlap_clean = overlap_results.reset_index(drop=True)
    
    merged = interactions_clean.copy()
    for col in overlap_clean.columns:
        if col != 'interaction_idx':
            merged[col] = overlap_clean[col].values
    
    # Split by direction
    jw18_uninf = merged[merged['logFC'] > 0]
    jw18_wmel = merged[merged['logFC'] < 0]
    
    # Calculate peaks per interaction
    results = {
        'jw18_uninf': {
            'n_interactions': len(jw18_uninf),
            'total_peaks': jw18_uninf['total_peaks'].sum(),
            'peaks_per_interaction': jw18_uninf['total_peaks'].sum() / len(jw18_uninf) if len(jw18_uninf) > 0 else 0,
            'mean_peaks': jw18_uninf['total_peaks'].mean(),
            'median_peaks': jw18_uninf['total_peaks'].median(),
            'std_peaks': jw18_uninf['total_peaks'].std()
        },
        'jw18_wmel': {
            'n_interactions': len(jw18_wmel),
            'total_peaks': jw18_wmel['total_peaks'].sum(),
            'peaks_per_interaction': jw18_wmel['total_peaks'].sum() / len(jw18_wmel) if len(jw18_wmel) > 0 else 0,
            'mean_peaks': jw18_wmel['total_peaks'].mean(),
            'median_peaks': jw18_wmel['total_peaks'].median(),
            'std_peaks': jw18_wmel['total_peaks'].std()
        }
    }
    
    # Statistical test: Mann-Whitney U for peaks per interaction
    if len(jw18_uninf) > 0 and len(jw18_wmel) > 0:
        try:
            u_stat, mw_p = stats.mannwhitneyu(
                jw18_uninf['total_peaks'], 
                jw18_wmel['total_peaks'], 
                alternative='two-sided'
            )
            
            # Effect size (rank-biserial correlation)
            n1 = len(jw18_uninf)
            n2 = len(jw18_wmel)
            rank_biserial = 1 - (2*u_stat) / (n1 * n2)
            
            # Cohen's d (effect size for continuous data)
            pooled_std = np.sqrt(((n1-1)*results['jw18_uninf']['std_peaks']**2 + 
                                   (n2-1)*results['jw18_wmel']['std_peaks']**2) / (n1+n2-2))
            cohens_d = (results['jw18_uninf']['mean_peaks'] - results['jw18_wmel']['mean_peaks']) / pooled_std if pooled_std > 0 else np.nan
            
        except Exception as e:
            print(f"Warning: Could not compute statistics: {e}")
            mw_p = np.nan
            rank_biserial = np.nan
            cohens_d = np.nan
    else:
        mw_p = np.nan
        rank_biserial = np.nan
        cohens_d = np.nan
    
    results['comparison'] = {
        'mannwhitneyu_p_value': mw_p,
        'effect_size_rank_biserial': rank_biserial,
        'cohens_d': cohens_d
    }
    
    print(f"\nPeaks per interaction results:")
    print(f"JW18 uninf.:")
    print(f"  N interactions: {results['jw18_uninf']['n_interactions']}")
    print(f"  Total peaks: {results['jw18_uninf']['total_peaks']}")
    print(f"  Peaks/interaction: {results['jw18_uninf']['peaks_per_interaction']:.2f}")
    print(f"  Mean ± SD: {results['jw18_uninf']['mean_peaks']:.2f} ± {results['jw18_uninf']['std_peaks']:.2f}")
    print(f"  Median: {results['jw18_uninf']['median_peaks']:.1f}")
    print(f"JW18 wMel:")
    print(f"  N interactions: {results['jw18_wmel']['n_interactions']}")
    print(f"  Total peaks: {results['jw18_wmel']['total_peaks']}")
    print(f"  Peaks/interaction: {results['jw18_wmel']['peaks_per_interaction']:.2f}")
    print(f"  Mean ± SD: {results['jw18_wmel']['mean_peaks']:.2f} ± {results['jw18_wmel']['std_peaks']:.2f}")
    print(f"  Median: {results['jw18_wmel']['median_peaks']:.1f}")
    print(f"\nStatistical comparison:")
    print(f"  Mann-Whitney U p-value: {mw_p:.2e}")
    if not np.isnan(rank_biserial):
        print(f"  Effect size (rank-biserial): {rank_biserial:.3f}")
    if not np.isnan(cohens_d):
        print(f"  Effect size (Cohen's d): {cohens_d:.3f}")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Analyze H4K16ac enrichment at differential X chromosome contacts'
    )
    parser.add_argument('--interactions', required=True,
                       help='Significant differential interactions CSV file')
    parser.add_argument('--h4k16ac_peaks', required=True,
                       help='H4K16ac ChIP-seq peaks BED file')
    parser.add_argument('--window_size', type=int, default=1000,
                       help='Window size around interaction anchors (bp)')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significant interactions')
    parser.add_argument('--n_permutations', type=int, default=1000,
                       help='Number of permutations for null model')
    parser.add_argument('--n_cores', type=int, default=None,
                       help='Number of CPU cores to use (default: all available - 1)')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix (including directory path)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("H4K16ac Enrichment Analysis at X Chromosome Contacts")
    print("OPTIMIZED VERSION with parallel processing")
    print("="*60)
    
    # Load data
    h4k16ac_peaks = load_h4k16ac_peaks(args.h4k16ac_peaks)
    sig_interactions, all_x_interactions = load_differential_interactions(
        args.interactions, args.fdr_threshold
    )
    
    if h4k16ac_peaks is None or len(sig_interactions) == 0:
        print("Error: No data to analyze")
        return 1
    
    # Calculate H4K16ac overlap for significant interactions
    print(f"\nAnalyzing H4K16ac overlap (window size: {args.window_size}bp)")
    overlap_results = calculate_h4k16ac_overlap_vectorized(
        sig_interactions, h4k16ac_peaks, args.window_size
    )
    
    if overlap_results is None:
        print("Error: Could not calculate overlap")
        return 1
    
    print(f"\nOverall H4K16ac overlap statistics:")
    print(f"  Interactions with at least one anchor overlapping: {overlap_results['any_anchor_overlap'].sum()} ({overlap_results['any_anchor_overlap'].mean()*100:.1f}%)")
    print(f"  Interactions with both anchors overlapping: {overlap_results['both_anchors_overlap'].sum()} ({overlap_results['both_anchors_overlap'].mean()*100:.1f}%)")
    print(f"  Mean H4K16ac peaks per interaction: {overlap_results['total_peaks'].mean():.2f}")
    
    # Analyze by logFC direction
    direction_results, merged_df = analyze_by_logfc_direction(
        sig_interactions, overlap_results
    )
    
    # Perform permutation test
    perm_results = permutation_test_h4k16ac_enrichment(
        sig_interactions, h4k16ac_peaks, overlap_results,
        window_size=args.window_size,
        n_permutations=args.n_permutations,
        n_cores=args.n_cores
    )
    
    # Create visualizations
    create_visualization(merged_df, direction_results, args.output_prefix)
    create_permutation_visualization(perm_results, args.output_prefix)
    
    # Save results
    output_csv = f"{args.output_prefix}_interactions.csv"
    merged_df.to_csv(output_csv, index=False)
    print(f"\nDetailed results saved to {output_csv}")
    
    # Save summary
    summary = {
        'total_x_interactions': len(all_x_interactions),
        'significant_x_interactions': len(sig_interactions),
        'h4k16ac_peaks_on_x': len(h4k16ac_peaks),
        'observed_overlap_rate': perm_results['observed_rate'],
        'null_mean_overlap_rate': perm_results['null_mean'],
        'null_std_overlap_rate': perm_results['null_std'],
        'z_score': perm_results['z_score'],
        'permutation_pvalue': perm_results['p_value'],
        'jw18_uninf_n': direction_results['jw18_uninf']['n_interactions'],
        'jw18_uninf_overlap_rate': direction_results['jw18_uninf']['overlap_rate'],
        'jw18_uninf_mean_peaks': direction_results['jw18_uninf']['mean_peaks'],
        'jw18_uninf_median_peaks': direction_results['jw18_uninf']['median_peaks'],
        'jw18_wmel_n': direction_results['jw18_wmel']['n_interactions'],
        'jw18_wmel_overlap_rate': direction_results['jw18_wmel']['overlap_rate'],
        'jw18_wmel_mean_peaks': direction_results['jw18_wmel']['mean_peaks'],
        'jw18_wmel_median_peaks': direction_results['jw18_wmel']['median_peaks'],
        'overlap_fisher_pvalue': direction_results.get('comparison', {}).get('fisher_p_value', None),
        'peak_count_mannwhitney_pvalue': direction_results.get('comparison', {}).get('mannwhitneyu_p_value', None),
        'peak_count_effect_size': direction_results.get('comparison', {}).get('effect_size_rank_biserial', None)
    }
    
    summary_df = pd.DataFrame([summary])
    summary_csv = f"{args.output_prefix}_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"Summary statistics saved to {summary_csv}")
    
    # Save text summary
    summary_txt = f"{args.output_prefix}_summary.txt"
    with open(summary_txt, 'w') as f:
        f.write("H4K16ac ENRICHMENT ANALYSIS SUMMARY\n")
        f.write("="*60 + "\n\n")
        f.write(f"Analysis Parameters:\n")
        f.write(f"  Window size: {args.window_size} bp\n")
        f.write(f"  FDR threshold: {args.fdr_threshold}\n")
        f.write(f"  Number of permutations: {args.n_permutations}\n\n")
        
        f.write(f"Data Summary:\n")
        f.write(f"  Total X chromosome interactions: {summary['total_x_interactions']}\n")
        f.write(f"  Significant interactions: {summary['significant_x_interactions']}\n")
        f.write(f"  H4K16ac peaks on X: {summary['h4k16ac_peaks_on_x']}\n\n")
        
        f.write(f"Permutation Test Results:\n")
        f.write(f"  Observed overlap rate: {summary['observed_overlap_rate']*100:.1f}%\n")
        f.write(f"  Null mean overlap rate: {summary['null_mean_overlap_rate']*100:.1f}%\n")
        f.write(f"  Null std: {summary['null_std_overlap_rate']*100:.2f}%\n")
        f.write(f"  Z-score: {summary['z_score']:.2f}\n")
        f.write(f"  P-value: {summary['permutation_pvalue']:.4f}\n")
        f.write(f"  Significant: {'YES' if summary['permutation_pvalue'] < 0.05 else 'NO'}\n\n")
        
        f.write(f"Genotype-specific Analysis:\n")
        f.write(f"  JW18 uninf. (up-regulated) interactions:\n")
        f.write(f"    N = {summary['jw18_uninf_n']}\n")
        f.write(f"    H4K16ac overlap: {summary['jw18_uninf_overlap_rate']*100:.1f}%\n")
        f.write(f"    Mean H4K16ac peaks: {summary['jw18_uninf_mean_peaks']:.2f}\n")
        f.write(f"    Median H4K16ac peaks: {summary['jw18_uninf_median_peaks']:.1f}\n")
        f.write(f"  JW18 wMel (down-regulated) interactions:\n")
        f.write(f"    N = {summary['jw18_wmel_n']}\n")
        f.write(f"    H4K16ac overlap: {summary['jw18_wmel_overlap_rate']*100:.1f}%\n")
        f.write(f"    Mean H4K16ac peaks: {summary['jw18_wmel_mean_peaks']:.2f}\n")
        f.write(f"    Median H4K16ac peaks: {summary['jw18_wmel_median_peaks']:.1f}\n")
        f.write(f"\n")
        f.write(f"  Statistical tests:\n")
        if summary['overlap_fisher_pvalue'] is not None:
            f.write(f"    Overlap rate comparison (Fisher's): p = {summary['overlap_fisher_pvalue']:.2e}\n")
        if summary['peak_count_mannwhitney_pvalue'] is not None:
            f.write(f"    Peak count comparison (Mann-Whitney U): p = {summary['peak_count_mannwhitney_pvalue']:.2e}\n")
        if summary['peak_count_effect_size'] is not None:
            f.write(f"    Effect size (rank-biserial): {summary['peak_count_effect_size']:.3f}\n")
    
    print(f"Text summary saved to {summary_txt}")
    print(f"\nAnalysis complete! All results saved with prefix: {args.output_prefix}")
    
    return 0

if __name__ == '__main__':
    exit(main())