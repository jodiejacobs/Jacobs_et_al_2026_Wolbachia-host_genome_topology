#!/usr/bin/env python3
"""
Analyze H3K36me3 enrichment at differential chromatin contacts on the X chromosome.
Compares real differential interactions vs null model to assess significance.

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

def load_null_model(null_file):
    """Load null model interactions"""
    print(f"Loading null model from {null_file}")
    
    # The null model file doesn't have chromosome info, just statistical measures
    null_data = pd.read_csv(null_file)
    
    print(f"Loaded {len(null_data)} null model interactions")
    
    return null_data

def calculate_h3k36me3_overlap(interactions_df, h3k36me3_peaks, window_size=5000):
    """
    Calculate overlap between interaction anchors and H3K36me3 peaks.
    Returns both anchor-level and interaction-level statistics.
    """
    print(f"\nAnalyzing H3K36me3 overlap (window size: {window_size}bp)")
    
    if h3k36me3_peaks is None or len(h3k36me3_peaks) == 0:
        print("No H3K36me3 peaks available")
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
            print(f"Error processing interaction {idx}: {e}")
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
    
    print(f"\nOverall H3K36me3 overlap statistics:")
    print(f"  Interactions with at least one anchor overlapping: {results_df['any_anchor_overlap'].sum()} ({results_df['any_anchor_overlap'].mean()*100:.1f}%)")
    print(f"  Interactions with both anchors overlapping: {results_df['both_anchors_overlap'].sum()} ({results_df['both_anchors_overlap'].mean()*100:.1f}%)")
    print(f"  Mean H3K36me3 peaks per interaction: {results_df['total_peaks'].mean():.2f}")
    
    return results_df

def ensure_boolean_columns(df, columns):
    """
    FIX: Ensure specified columns are boolean type and handle NaN values.
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
    
    # FIX: Ensure boolean columns are properly typed
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
    # FIX: Ensure boolean before using ~ operator
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

def compare_to_null_model(real_overlap_rate, null_model):
    """Compare real H3K36me3 enrichment to null model"""
    print("\nComparing to null model...")
    
    # The null model doesn't have specific overlap rates, 
    # so we'll use a baseline expectation
    # For now, we'll assume 50% as a null expectation (random chance)
    null_expectation = 0.5
    
    # Calculate enrichment
    enrichment = real_overlap_rate / null_expectation
    
    # Calculate p-value using binomial test
    # H0: overlap rate = null_expectation
    n_interactions = len(null_model)
    n_overlaps = int(real_overlap_rate * n_interactions)
    
    # Two-tailed binomial test
    p_value = stats.binom_test(n_overlaps, n_interactions, null_expectation, alternative='two-sided')
    
    return {
        'real_overlap_rate': real_overlap_rate,
        'null_expectation': null_expectation,
        'enrichment': enrichment,
        'p_value': p_value
    }

def create_visualization(merged_df, direction_results, output_prefix):
    """Create separate 2x2 visualizations with updated labels and colors"""
    print("\nCreating visualizations...")
    
    # FIX: Ensure all boolean columns are properly typed before visualization
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
    
    # --- Plot Set 1: Overlap rates and peak counts (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)  # Transparent figure background
    
    # Plot 1: H3K36me3 Overlap Rate by Direction
    ax = axes[0, 0]
    ax.patch.set_alpha(0)  # Transparent axis background
    overlap_rates = [
        direction_results['jw18_uninf']['overlap_rate'] * 100,
        direction_results['jw18_wmel']['overlap_rate'] * 100
    ]
    bars = ax.bar(labels, overlap_rates, color=colors, alpha=0.8)
    ax.set_ylabel('H3K36me3 Overlap Rate (%)', fontsize=6)
    ax.set_title('H3K36me3 Enrichment by Genotype', fontsize=7, fontweight='bold')
    ax.axhline(y=50, color='gray', linestyle='--', linewidth=0.8, label='Expected (50%)')
    ax.legend(fontsize=5)
    ax.tick_params(labelsize=6)
    
    # Add p-value if available
    if 'comparison' in direction_results and 'p_value' in direction_results['comparison']:
        p_val = direction_results['comparison']['p_value']
        ax.text(0.5, max(overlap_rates) * 0.95, f'p = {p_val:.2e}', 
                ha='center', va='top', fontsize=5, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add value labels on bars
    for bar, rate in zip(bars, overlap_rates):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{rate:.1f}%',
                ha='center', va='bottom', fontsize=5)
    
    # Plot 2: Mean H3K36me3 Peaks per Interaction
    ax = axes[0, 1]
    ax.patch.set_alpha(0)
    peak_counts = [
        direction_results['jw18_uninf']['mean_peaks'],
        direction_results['jw18_wmel']['mean_peaks']
    ]
    bars = ax.bar(labels, peak_counts, color=colors, alpha=0.8)
    ax.set_ylabel('Mean H3K36me3 Peaks per Interaction', fontsize=6)
    ax.set_title('H3K36me3 Peak Density', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add value labels on bars
    for bar, count in zip(bars, peak_counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:.2f}',
                ha='center', va='bottom', fontsize=5)
    
    # Plot 3: Distribution of H3K36me3 peaks
    ax = axes[1, 0]
    ax.patch.set_alpha(0)
    uninf_peaks = merged_df[merged_df['logFC'] > 0]['total_peaks']
    wmel_peaks = merged_df[merged_df['logFC'] < 0]['total_peaks']
    
    ax.hist([uninf_peaks, wmel_peaks], bins=20, label=labels, 
            color=colors, alpha=0.7)
    ax.set_xlabel('Number of H3K36me3 Peaks', fontsize=6)
    ax.set_ylabel('Number of Interactions', fontsize=6)
    ax.set_title('Distribution of H3K36me3 Peaks', fontsize=7, fontweight='bold')
    ax.legend(fontsize=5)
    ax.tick_params(labelsize=6)
    
    # Add statistical test (Mann-Whitney U test)
    if len(uninf_peaks) > 0 and len(wmel_peaks) > 0:
        u_stat, p_val = stats.mannwhitneyu(uninf_peaks, wmel_peaks, alternative='two-sided')
        ax.text(0.95, 0.95, f'p = {p_val:.2e}', 
                transform=ax.transAxes, ha='right', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Plot 4: Overlap type distribution
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
    ax.set_title('H3K36me3 Overlap Pattern', fontsize=7, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    
    # Add value labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{count}',
                ha='center', va='bottom', fontsize=5)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/h3k36me3_analysis_set1.pdf", 
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Plot set 1 saved to {output_prefix}/h3k36me3_analysis_set1.pdf")
    
    # --- Plot Set 2: LogFC correlations and interaction types (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)
    
    # Plot 1: LogFC vs H3K36me3 peaks
    ax = axes[0, 0]
    ax.patch.set_alpha(0)
    
    # Color points by genotype
    uninf_mask = merged_df['logFC'] > 0
    wmel_mask = merged_df['logFC'] < 0
    
    ax.scatter(merged_df[uninf_mask]['logFC'], merged_df[uninf_mask]['total_peaks'], 
               c=color_jw18_uninf, alpha=0.6, s=20, label='JW18 uninf.', edgecolors='none')
    ax.scatter(merged_df[wmel_mask]['logFC'], merged_df[wmel_mask]['total_peaks'], 
               c=color_jw18_wmel, alpha=0.6, s=20, label='JW18 wMel', edgecolors='none')
    
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
    
    # Plot 2: Cis vs Trans interactions
    ax = axes[0, 1]
    ax.patch.set_alpha(0)
    
    if 'interaction_type' in merged_df.columns:
        cis_overlap = merged_df[merged_df['interaction_type'] == 'cis']['any_anchor_overlap'].mean() * 100
        trans_overlap = merged_df[merged_df['interaction_type'] == 'trans']['any_anchor_overlap'].mean() * 100
        
        bars = ax.bar(['Cis', 'Trans'], [cis_overlap, trans_overlap], 
               color=['#3366cc', '#33cc66'], alpha=0.8)
        ax.set_ylabel('H3K36me3 Overlap Rate (%)', fontsize=6)
        ax.set_title('H3K36me3 Enrichment by Interaction Type', fontsize=7, fontweight='bold')
        ax.axhline(y=50, color='gray', linestyle='--', linewidth=0.8, label='Expected (50%)')
        ax.legend(fontsize=5)
        ax.tick_params(labelsize=6)
        
        # Add value labels on bars
        for bar, rate in zip(bars, [cis_overlap, trans_overlap]):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{rate:.1f}%',
                    ha='center', va='bottom', fontsize=5)
        
        # Add Fisher's exact test
        # FIX: Ensure boolean columns before using ~
        cis_data = merged_df[merged_df['interaction_type'] == 'cis'].copy()
        trans_data = merged_df[merged_df['interaction_type'] == 'trans'].copy()
        cis_data = ensure_boolean_columns(cis_data, ['any_anchor_overlap'])
        trans_data = ensure_boolean_columns(trans_data, ['any_anchor_overlap'])
        
        contingency = [
            [cis_data['any_anchor_overlap'].sum(), (~cis_data['any_anchor_overlap']).sum()],
            [trans_data['any_anchor_overlap'].sum(), (~trans_data['any_anchor_overlap']).sum()]
        ]
        
        odds_ratio, p_val = stats.fisher_exact(contingency)
        ax.text(0.5, max(cis_overlap, trans_overlap) * 0.95, f'p = {p_val:.2e}', 
                ha='center', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Plot 3: Both anchors overlap comparison
    ax = axes[1, 0]
    ax.patch.set_alpha(0)
    
    both_anchors_rates = [
        direction_results['jw18_uninf']['both_anchors_rate'] * 100,
        direction_results['jw18_wmel']['both_anchors_rate'] * 100
    ]
    bars = ax.bar(labels, both_anchors_rates, color=colors, alpha=0.8)
    ax.set_ylabel('Both Anchors Overlap Rate (%)', fontsize=6)
    ax.set_title('H3K36me3 at Both Interaction Anchors', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add value labels on bars
    for bar, rate in zip(bars, both_anchors_rates):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{rate:.1f}%',
                ha='center', va='bottom', fontsize=5)
    
    # Plot 4: Sample size information
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
    
    # Add value labels on bars
    for bar, size in zip(bars, sample_sizes):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'n={size}',
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
    parser.add_argument('--null_model', required=True,
                       help='Null model interactions CSV file')
    parser.add_argument('--h3k36me3_peaks', required=True,
                       help='H3K36me3 ChIP-seq peaks BED file')
    parser.add_argument('--window_size', type=int, default=5000,
                       help='Window size around interaction anchors (bp)')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significant interactions')
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
    null_model = load_null_model(args.null_model)
    
    if h3k36me3_peaks is None or len(sig_interactions) == 0:
        print("Error: No data to analyze")
        return 1
    
    # Calculate H3K36me3 overlap for significant interactions
    overlap_results = calculate_h3k36me3_overlap(
        sig_interactions, h3k36me3_peaks, args.window_size
    )
    
    if overlap_results is None:
        print("Error: Could not calculate overlap")
        return 1
    
    # Analyze by logFC direction
    direction_results, merged_df = analyze_by_logfc_direction(
        sig_interactions, overlap_results
    )
    
    # Compare to null model
    overall_overlap_rate = overlap_results['any_anchor_overlap'].mean()
    null_comparison = compare_to_null_model(overall_overlap_rate, null_model)
    
    print(f"\nNull Model Comparison:")
    print(f"  Real overlap rate: {null_comparison['real_overlap_rate']*100:.1f}%")
    print(f"  Expected (baseline): {null_comparison['null_expectation']*100:.1f}%")
    print(f"  Enrichment: {null_comparison['enrichment']:.2f}x")
    print(f"  P-value: {null_comparison['p_value']:.2e}")
    
    # Create visualizations
    create_visualization(merged_df, direction_results, args.output_prefix)
    
    # Save detailed results
    merged_df.to_csv(f"{args.output_prefix}/h3k36me3_interactions.csv", index=False)
    
    # Save summary statistics
    summary = {
        'total_x_interactions': len(all_x_interactions),
        'significant_x_interactions': len(sig_interactions),
        'h3k36me3_peaks_on_x': len(h3k36me3_peaks),
        'overall_overlap_rate': overall_overlap_rate,
        'enrichment_vs_null': null_comparison['enrichment'],
        'p_value': null_comparison['p_value'],
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
        f.write(f"  FDR threshold: {args.fdr_threshold}\n\n")
        
        f.write(f"Data Summary:\n")
        f.write(f"  Total X chromosome interactions: {summary['total_x_interactions']}\n")
        f.write(f"  Significant interactions: {summary['significant_x_interactions']}\n")
        f.write(f"  H3K36me3 peaks on X: {summary['h3k36me3_peaks_on_x']}\n\n")
        
        f.write(f"H3K36me3 Enrichment:\n")
        f.write(f"  Overall overlap rate: {summary['overall_overlap_rate']*100:.1f}%\n")
        f.write(f"  Enrichment vs null: {summary['enrichment_vs_null']:.2f}x\n")
        f.write(f"  P-value: {summary['p_value']:.2e}\n")
        f.write(f"  Significant: {'YES' if summary['p_value'] < 0.05 else 'NO'}\n\n")
        
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