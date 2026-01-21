#!/usr/bin/env python3
"""
Analyze H3K36me3 enrichment at differential chromatin contacts on the X chromosome.
Compares real differential interactions vs permutation-based null model to assess significance.

Based on:
- Fei et al. - NDF associates with MSL complex via H3K36me3
- Larschan et al. - H3K36me3 marks active genes for dosage compensation
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

def load_h3k36me3_peaks(chip_file):
    """
    Load H3K36me3 ChIP-seq peaks from GSE20784 or similar datasets.
    Expects BED format: chr, start, end, score
    Format: chr2L 67843 71268 1.90337415112343
    """
    print(f"Loading H3K36me3 peaks from {chip_file}")
    
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
        
        print(f"Loaded {len(x_peaks)} H3K36me3 peaks on X chromosome")
        print(f"Score range: {x_peaks['score'].min():.2f} - {x_peaks['score'].max():.2f}")
        
        return x_peaks
        
    except Exception as e:
        print(f"Error loading H3K36me3 peaks: {e}")
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
    
    return sig_x, x_interactions

def calculate_h3k36me3_overlap(interactions_df, h3k36me3_peaks, window_size=5000):
    """
    Calculate overlap between interaction anchors and H3K36me3 peaks.
    Returns both anchor-level and interaction-level statistics.
    """
    if h3k36me3_peaks is None or len(h3k36me3_peaks) == 0:
        return None
    
    # Create BedTool for H3K36me3 peaks
    h3k36_bt = pybedtools.BedTool.from_dataframe(
        h3k36me3_peaks[['chrom', 'start', 'end']]
    )
    
    results = []
    
    for idx, interaction in interactions_df.iterrows():
        # Extract anchor coordinates
        anchor1_chr = interaction['chr1']
        anchor1_start = interaction['start1']
        anchor1_end = interaction['end1']
        
        anchor2_chr = interaction['chr2']
        anchor2_start = interaction['start2']
        anchor2_end = interaction['end2']
        
        # Create windows around anchors
        anchor1_window = pd.DataFrame([{
            'chrom': anchor1_chr,
            'start': max(0, anchor1_start - window_size),
            'end': anchor1_end + window_size
        }])
        
        anchor2_window = pd.DataFrame([{
            'chrom': anchor2_chr,
            'start': max(0, anchor2_start - window_size),
            'end': anchor2_end + window_size
        }])
        
        # Check for overlap
        try:
            anchor1_bt = pybedtools.BedTool.from_dataframe(anchor1_window)
            anchor2_bt = pybedtools.BedTool.from_dataframe(anchor2_window)
            
            anchor1_overlap = len(anchor1_bt.intersect(h3k36_bt)) > 0
            anchor2_overlap = len(anchor2_bt.intersect(h3k36_bt)) > 0
            
            # Count number of H3K36me3 peaks overlapping
            n_peaks_anchor1 = len(anchor1_bt.intersect(h3k36_bt, wa=True))
            n_peaks_anchor2 = len(anchor2_bt.intersect(h3k36_bt, wa=True))
            
            results.append({
                'interaction_idx': idx,
                'anchor1_overlap': anchor1_overlap,
                'anchor2_overlap': anchor2_overlap,
                'both_anchors_overlap': anchor1_overlap and anchor2_overlap,
                'any_anchor_overlap': anchor1_overlap or anchor2_overlap,
                'n_peaks_anchor1': n_peaks_anchor1,
                'n_peaks_anchor2': n_peaks_anchor2,
                'total_peaks': n_peaks_anchor1 + n_peaks_anchor2
            })
        except Exception as e:
            results.append({
                'interaction_idx': idx,
                'anchor1_overlap': False,
                'anchor2_overlap': False,
                'both_anchors_overlap': False,
                'any_anchor_overlap': False,
                'n_peaks_anchor1': 0,
                'n_peaks_anchor2': 0,
                'total_peaks': 0
            })
    
    results_df = pd.DataFrame(results)
    
    return results_df

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
    """Analyze H3K36me3 enrichment by logFC direction (JW18 uninf. vs JW18 wMel)"""
    print("\nAnalyzing by logFC direction...")
    
    # Merge interactions with overlap results
    merged = interactions_df.copy()
    merged['overlap_idx'] = range(len(merged))
    merged = merged.merge(overlap_results, left_on='overlap_idx', right_on='interaction_idx', how='left')
    
    # Ensure boolean columns are properly typed
    bool_cols = ['anchor1_overlap', 'anchor2_overlap', 'both_anchors_overlap', 'any_anchor_overlap']
    merged = ensure_boolean_columns(merged, bool_cols)
    
    # Add interaction type
    merged['interaction_type'] = np.where(merged['chr1'] == merged['chr2'], 'cis', 'trans')
    
    # Split by direction - JW18 uninf. (positive logFC) vs JW18 wMel (negative logFC)
    jw18_uninf = merged[merged['logFC'] > 0]  # Up-regulated = JW18 uninf.
    jw18_wmel = merged[merged['logFC'] < 0]   # Down-regulated = JW18 wMel
    
    # Calculate statistics
    results = {
        'jw18_uninf': {
            'n_interactions': len(jw18_uninf),
            'overlap_rate': jw18_uninf['any_anchor_overlap'].mean(),
            'both_anchors_rate': jw18_uninf['both_anchors_overlap'].mean(),
            'mean_peaks': jw18_uninf['total_peaks'].mean()
        },
        'jw18_wmel': {
            'n_interactions': len(jw18_wmel),
            'overlap_rate': jw18_wmel['any_anchor_overlap'].mean(),
            'both_anchors_rate': jw18_wmel['both_anchors_overlap'].mean(),
            'mean_peaks': jw18_wmel['total_peaks'].mean()
        }
    }
    
    # Statistical comparison using Fisher's exact test
    contingency = [
        [jw18_uninf['any_anchor_overlap'].sum(), 
         (~jw18_uninf['any_anchor_overlap']).sum()],
        [jw18_wmel['any_anchor_overlap'].sum(), 
         (~jw18_wmel['any_anchor_overlap']).sum()]
    ]
    
    odds_ratio, p_value = stats.fisher_exact(contingency)
    
    results['comparison'] = {
        'odds_ratio': odds_ratio,
        'p_value': p_value
    }
    
    print(f"\nDirection-specific results:")
    print(f"JW18 uninf. (n={results['jw18_uninf']['n_interactions']}):")
    print(f"  H3K36me3 overlap rate: {results['jw18_uninf']['overlap_rate']*100:.1f}%")
    print(f"  Mean peaks: {results['jw18_uninf']['mean_peaks']:.2f}")
    print(f"JW18 wMel (n={results['jw18_wmel']['n_interactions']}):")
    print(f"  H3K36me3 overlap rate: {results['jw18_wmel']['overlap_rate']*100:.1f}%")
    print(f"  Mean peaks: {results['jw18_wmel']['mean_peaks']:.2f}")
    print(f"Fisher's exact test p-value: {p_value:.2e}")
    
    return results, merged

def permutation_test_h3k36me3_enrichment(interactions_df, h3k36me3_peaks, 
                                          overlap_results, window_size=5000, 
                                          n_permutations=1000):
    """
    Perform permutation test by shuffling interaction anchor positions
    on the X chromosome and recalculating H3K36me3 overlap.
    
    Returns:
        dict with observed overlap rate, null distribution, p-value, and z-score
    """
    print(f"\nPerforming permutation test ({n_permutations} permutations)...")
    
    # Get observed overlap rate
    observed_rate = overlap_results['any_anchor_overlap'].mean()
    
    # Get X chromosome size from interactions
    x_chrom_end = max(
        interactions_df['end1'].max(),
        interactions_df['end2'].max()
    )
    
    print(f"X chromosome size for permutations: {x_chrom_end:,} bp")
    print(f"Observed overlap rate: {observed_rate*100:.1f}%")
    
    # Store null distribution
    null_overlap_rates = []
    
    for perm_idx in range(n_permutations):
        if (perm_idx + 1) % 100 == 0:
            print(f"  Permutation {perm_idx + 1}/{n_permutations}...")
        
        # Shuffle anchor positions - maintain anchor sizes
        permuted_interactions = interactions_df.copy()
        
        for idx in permuted_interactions.index:
            # Anchor 1
            anchor1_size = permuted_interactions.loc[idx, 'end1'] - permuted_interactions.loc[idx, 'start1']
            new_start1 = np.random.randint(0, max(1, x_chrom_end - anchor1_size))
            permuted_interactions.loc[idx, 'start1'] = new_start1
            permuted_interactions.loc[idx, 'end1'] = new_start1 + anchor1_size
            
            # Anchor 2
            anchor2_size = permuted_interactions.loc[idx, 'end2'] - permuted_interactions.loc[idx, 'start2']
            new_start2 = np.random.randint(0, max(1, x_chrom_end - anchor2_size))
            permuted_interactions.loc[idx, 'start2'] = new_start2
            permuted_interactions.loc[idx, 'end2'] = new_start2 + anchor2_size
        
        # Calculate overlap for permuted data
        perm_overlap = calculate_h3k36me3_overlap(
            permuted_interactions, h3k36me3_peaks, window_size
        )
        
        if perm_overlap is not None:
            null_overlap_rates.append(perm_overlap['any_anchor_overlap'].mean())
    
    null_overlap_rates = np.array(null_overlap_rates)
    
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
    
    ax.set_xlabel('H3K36me3 Overlap Rate (%)', fontsize=8)
    ax.set_ylabel('Frequency', fontsize=8)
    ax.set_title('Permutation Test: H3K36me3 Enrichment', fontsize=10, fontweight='bold')
    ax.legend(fontsize=7)
    ax.tick_params(labelsize=7)
    
    # Add statistics text box
    stats_text = (f"Z-score: {perm_results['z_score']:.2f}\n"
                  f"P-value: {perm_results['p_value']:.4f}")
    ax.text(0.05, 0.95, stats_text,
            transform=ax.transAxes, ha='left', va='top', fontsize=7,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/h3k36me3_permutation_test.pdf",
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Permutation plot saved to {output_prefix}/h3k36me3_permutation_test.pdf")

def add_strip_to_boxplot(ax, data, positions, colors, jitter=0.08, alpha=0.4, size=3):
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
        # Add jitter to x positions
        x = np.random.normal(pos, jitter, size=len(d))
        ax.scatter(x, d, alpha=alpha, s=size, color=color, edgecolors='none', zorder=3)

def create_visualization(merged_df, direction_results, output_prefix):
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
    
    # --- Plot Set 1: Overlap rates and peak counts (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)
    
    # Plot 1: H3K36me3 Peak Count by Genotype (BOXPLOT)
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
                         jitter=0.08, alpha=0.4, size=5)
    
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylabel('Number of H3K36me3 Peaks', fontsize=6)
    ax.set_title('H3K36me3 Peak Count by Genotype', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add Mann-Whitney U test
    if len(uninf_peaks) > 0 and len(wmel_peaks) > 0:
        u_stat, p_val = stats.mannwhitneyu(uninf_peaks, wmel_peaks, alternative='two-sided')
        ax.text(0.95, 0.95, f'p = {p_val:.2e}', 
                transform=ax.transAxes, ha='right', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Plot 2: H3K36me3 Overlap Rate (Binary - Bar plot makes more sense)
    ax = axes[0, 1]
    ax.patch.set_alpha(0)
    overlap_rates = [
        direction_results['jw18_uninf']['overlap_rate'] * 100,
        direction_results['jw18_wmel']['overlap_rate'] * 100
    ]
    bars = ax.bar([0, 1], overlap_rates, color=colors, alpha=0.8, width=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylabel('H3K36me3 Overlap Rate (%)', fontsize=6)
    ax.set_title('H3K36me3 Enrichment by Genotype', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add p-value if available
    if 'comparison' in direction_results and 'p_value' in direction_results['comparison']:
        p_val = direction_results['comparison']['p_value']
        ax.text(0.5, max(overlap_rates) * 0.95, f'p = {p_val:.2e}', 
                ha='center', va='top', fontsize=5, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add value labels on bars
    for i, (bar, rate) in enumerate(zip(bars, overlap_rates)):
        height = bar.get_height()
        ax.text(i, height,
                f'{rate:.1f}%',
                ha='center', va='bottom', fontsize=5)
    
    # Plot 3: LogFC vs H3K36me3 peaks (Scatter with boxplots on margins)
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
    ax.set_ylabel('Number of H3K36me3 Peaks', fontsize=6)
    ax.set_title('H3K36me3 Enrichment vs Effect Size', fontsize=7, fontweight='bold')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.legend(fontsize=5)
    ax.tick_params(labelsize=6)
    
    # Add correlation statistics
    corr, p_val = stats.spearmanr(merged_df['logFC'], merged_df['total_peaks'])
    ax.text(0.05, 0.95, f'ρ = {corr:.3f}\np = {p_val:.2e}', 
            transform=ax.transAxes, ha='left', va='top', fontsize=5,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Plot 4: Cis vs Trans interactions (BOXPLOT)
    ax = axes[1, 1]
    ax.patch.set_alpha(0)
    
    if 'interaction_type' in merged_df.columns:
        cis_peaks = merged_df[merged_df['interaction_type'] == 'cis']['total_peaks'].values
        trans_peaks = merged_df[merged_df['interaction_type'] == 'trans']['total_peaks'].values
        
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
                             jitter=0.08, alpha=0.4, size=5)
        
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Cis', 'Trans'])
        ax.set_ylabel('Number of H3K36me3 Peaks', fontsize=6)
        ax.set_title('H3K36me3 by Interaction Type', fontsize=7, fontweight='bold')
        ax.tick_params(labelsize=6)
        
        # Add Mann-Whitney U test
        if len(cis_peaks) > 0 and len(trans_peaks) > 0:
            u_stat, p_val = stats.mannwhitneyu(cis_peaks, trans_peaks, alternative='two-sided')
            ax.text(0.95, 0.95, f'p = {p_val:.2e}', 
                    transform=ax.transAxes, ha='right', va='top', fontsize=5,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/h3k36me3_analysis_set1.pdf", 
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Plot set 1 saved to {output_prefix}/h3k36me3_analysis_set1.pdf")
    
    # --- Plot Set 2: Additional analyses (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)
    
    # Plot 1: Both anchors overlap comparison (BOXPLOT)
    ax = axes[0, 0]
    ax.patch.set_alpha(0)
    
    uninf_both = merged_df[merged_df['logFC'] > 0]['both_anchors_overlap'].astype(int).values
    wmel_both = merged_df[merged_df['logFC'] < 0]['both_anchors_overlap'].astype(int).values
    
    # For binary data, show as bar plot with individual points
    both_anchors_rates = [
        direction_results['jw18_uninf']['both_anchors_rate'] * 100,
        direction_results['jw18_wmel']['both_anchors_rate'] * 100
    ]
    bars = ax.bar([0, 1], both_anchors_rates, color=colors, alpha=0.8, width=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylabel('Both Anchors Overlap Rate (%)', fontsize=6)
    ax.set_title('H3K36me3 at Both Interaction Anchors', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add value labels on bars
    for i, (bar, rate) in enumerate(zip(bars, both_anchors_rates)):
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
                         jitter=0.08, alpha=0.4, size=5)
    
    ax.set_xticks([0.5, 3])
    ax.set_xticklabels(['JW18 uninf.', 'JW18 wMel'])
    ax.set_ylabel('Number of H3K36me3 Peaks', fontsize=6)
    ax.set_title('H3K36me3 Peaks per Anchor', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add legend for anchors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='lightgray', edgecolor='black', label='Anchor 1'),
                      Patch(facecolor='darkgray', edgecolor='black', label='Anchor 2')]
    
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
    ax.set_title('H3K36me3 Overlap Pattern', fontsize=7, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    
    # Add value labels on bars
    for i, (bar, count) in enumerate(zip(bars, counts)):
        height = bar.get_height()
        ax.text(i, height,
                f'{count}',
                ha='center', va='bottom', fontsize=5)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/h3k36me3_analysis_set2.pdf", 
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Plot set 2 saved to {output_prefix}/h3k36me3_analysis_set2.pdf")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze H3K36me3 enrichment at differential X chromosome contacts'
    )
    parser.add_argument('--interactions', required=True,
                       help='Significant differential interactions CSV file')
    parser.add_argument('--h3k36me3_peaks', required=True,
                       help='H3K36me3 ChIP-seq peaks BED file')
    parser.add_argument('--window_size', type=int, default=5000,
                       help='Window size around interaction anchors (bp)')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significant interactions')
    parser.add_argument('--n_permutations', type=int, default=1000,
                       help='Number of permutations for null model')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("H3K36me3 Enrichment Analysis at X Chromosome Contacts")
    print("="*60)
    
    # Load data
    h3k36me3_peaks = load_h3k36me3_peaks(args.h3k36me3_peaks)
    sig_interactions, all_x_interactions = load_differential_interactions(
        args.interactions, args.fdr_threshold
    )
    
    if h3k36me3_peaks is None or len(sig_interactions) == 0:
        print("Error: No data to analyze")
        return 1
    
    # Calculate H3K36me3 overlap for significant interactions
    print(f"\nAnalyzing H3K36me3 overlap (window size: {args.window_size}bp)")
    overlap_results = calculate_h3k36me3_overlap(
        sig_interactions, h3k36me3_peaks, args.window_size
    )
    
    if overlap_results is None:
        print("Error: Could not calculate overlap")
        return 1
    
    print(f"\nOverall H3K36me3 overlap statistics:")
    print(f"  Interactions with at least one anchor overlapping: {overlap_results['any_anchor_overlap'].sum()} ({overlap_results['any_anchor_overlap'].mean()*100:.1f}%)")
    print(f"  Interactions with both anchors overlapping: {overlap_results['both_anchors_overlap'].sum()} ({overlap_results['both_anchors_overlap'].mean()*100:.1f}%)")
    print(f"  Mean H3K36me3 peaks per interaction: {overlap_results['total_peaks'].mean():.2f}")
    
    # Analyze by logFC direction
    direction_results, merged_df = analyze_by_logfc_direction(
        sig_interactions, overlap_results
    )
    
    # Perform permutation test
    perm_results = permutation_test_h3k36me3_enrichment(
        sig_interactions, h3k36me3_peaks, overlap_results,
        window_size=args.window_size,
        n_permutations=args.n_permutations
    )
    
    # Create visualizations
    create_visualization(merged_df, direction_results, args.output_prefix)
    create_permutation_visualization(perm_results, args.output_prefix)
    
    # Save detailed results
    merged_df.to_csv(f"{args.output_prefix}/h3k36me3_interactions.csv", index=False)
    
    # Save summary statistics
    summary = {
        'total_x_interactions': len(all_x_interactions),
        'significant_x_interactions': len(sig_interactions),
        'h3k36me3_peaks_on_x': len(h3k36me3_peaks),
        'observed_overlap_rate': perm_results['observed_rate'],
        'null_mean_overlap_rate': perm_results['null_mean'],
        'null_std_overlap_rate': perm_results['null_std'],
        'z_score': perm_results['z_score'],
        'permutation_pvalue': perm_results['p_value'],
        'jw18_uninf_n': direction_results['jw18_uninf']['n_interactions'],
        'jw18_uninf_overlap_rate': direction_results['jw18_uninf']['overlap_rate'],
        'jw18_wmel_n': direction_results['jw18_wmel']['n_interactions'],
        'jw18_wmel_overlap_rate': direction_results['jw18_wmel']['overlap_rate'],
        'direction_comparison_pvalue': direction_results.get('comparison', {}).get('p_value', None)
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(f"{args.output_prefix}/h3k36me3_summary.csv", index=False)
    
    # Create text summary
    with open(f"{args.output_prefix}/h3k36me3_summary.txt", 'w') as f:
        f.write("H3K36me3 ENRICHMENT ANALYSIS SUMMARY\n")
        f.write("="*60 + "\n\n")
        f.write(f"Analysis Parameters:\n")
        f.write(f"  Window size: {args.window_size} bp\n")
        f.write(f"  FDR threshold: {args.fdr_threshold}\n")
        f.write(f"  Number of permutations: {args.n_permutations}\n\n")
        
        f.write(f"Data Summary:\n")
        f.write(f"  Total X chromosome interactions: {summary['total_x_interactions']}\n")
        f.write(f"  Significant interactions: {summary['significant_x_interactions']}\n")
        f.write(f"  H3K36me3 peaks on X: {summary['h3k36me3_peaks_on_x']}\n\n")
        
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
        f.write(f"    H3K36me3 overlap: {summary['jw18_uninf_overlap_rate']*100:.1f}%\n")
        f.write(f"  JW18 wMel (down-regulated) interactions:\n")
        f.write(f"    N = {summary['jw18_wmel_n']}\n")
        f.write(f"    H3K36me3 overlap: {summary['jw18_wmel_overlap_rate']*100:.1f}%\n")
        
        if summary['direction_comparison_pvalue'] is not None:
            f.write(f"  Genotype comparison p-value: {summary['direction_comparison_pvalue']:.2e}\n")
    
    print(f"\nAnalysis complete! Results saved to {args.output_prefix}/*")
    
    return 0

if __name__ == '__main__':
    exit(main())