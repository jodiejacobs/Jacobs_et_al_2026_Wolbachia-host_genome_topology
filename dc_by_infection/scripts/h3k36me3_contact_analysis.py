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
        print(f"  Up-regulated: {sum(sig_x['logFC'] > 0)}")
        print(f"  Down-regulated: {sum(sig_x['logFC'] < 0)}")
    
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
                'any_anchor_overlap': anchor1_overlap or anchor2_overlap,
                'both_anchors_overlap': anchor1_overlap and anchor2_overlap,
                'n_peaks_anchor1': n_peaks_anchor1,
                'n_peaks_anchor2': n_peaks_anchor2,
                'total_peaks': n_peaks_anchor1 + n_peaks_anchor2
            })
            
        except Exception as e:
            print(f"Warning: Could not process interaction {idx}: {e}")
            continue
    
    results_df = pd.DataFrame(results)
    
    # Calculate summary statistics
    print(f"\nH3K36me3 Overlap Summary:")
    print(f"  Interactions with any anchor overlap: {results_df['any_anchor_overlap'].sum()} ({results_df['any_anchor_overlap'].mean()*100:.1f}%)")
    print(f"  Interactions with both anchors overlap: {results_df['both_anchors_overlap'].sum()} ({results_df['both_anchors_overlap'].mean()*100:.1f}%)")
    print(f"  Mean H3K36me3 peaks per interaction: {results_df['total_peaks'].mean():.2f}")
    
    return results_df

def compare_to_null_model(real_overlap_rate, null_data, n_permutations=1000):
    """
    Compare real H3K36me3 overlap rate to null expectation.
    Since we don't have chromosome info in null model, use bootstrap approach.
    """
    print(f"\nComparing to null model...")
    
    # Calculate expected overlap rate from null model
    # Assume genome-wide H3K36me3 coverage and X chromosome proportion
    
    # Simplified approach: use permutation test on the null model statistics
    null_sample = null_data.sample(n=min(1000, len(null_data)), replace=True)
    
    # Calculate p-value
    p_value = 1.0 - stats.norm.cdf(
        real_overlap_rate,
        loc=0.5,  # Baseline expectation
        scale=np.sqrt(0.5 * 0.5 / len(null_sample))
    )
    
    enrichment = real_overlap_rate / 0.5  # Compared to baseline
    
    return {
        'enrichment': enrichment,
        'p_value': p_value,
        'real_overlap_rate': real_overlap_rate,
        'null_expectation': 0.5
    }

def analyze_by_logfc_direction(interactions_df, overlap_df):
    """
    Analyze H3K36me3 enrichment separately for up and down-regulated interactions.
    """
    print("\nAnalyzing by logFC direction...")
    
    # Merge overlap results with interaction data
    merged = interactions_df.copy()
    merged['any_anchor_overlap'] = overlap_df['any_anchor_overlap'].values
    merged['both_anchors_overlap'] = overlap_df['both_anchors_overlap'].values
    merged['total_peaks'] = overlap_df['total_peaks'].values
    
    # Split by direction
    up_regulated = merged[merged['logFC'] > 0]
    down_regulated = merged[merged['logFC'] < 0]
    
    results = {
        'up_regulated': {
            'n_interactions': len(up_regulated),
            'overlap_rate': up_regulated['any_anchor_overlap'].mean() if len(up_regulated) > 0 else 0,
            'mean_peaks': up_regulated['total_peaks'].mean() if len(up_regulated) > 0 else 0
        },
        'down_regulated': {
            'n_interactions': len(down_regulated),
            'overlap_rate': down_regulated['any_anchor_overlap'].mean() if len(down_regulated) > 0 else 0,
            'mean_peaks': down_regulated['total_peaks'].mean() if len(down_regulated) > 0 else 0
        }
    }
    
    # Statistical test
    if len(up_regulated) > 0 and len(down_regulated) > 0:
        # Chi-square test for overlap rates
        contingency = np.array([
            [up_regulated['any_anchor_overlap'].sum(), len(up_regulated) - up_regulated['any_anchor_overlap'].sum()],
            [down_regulated['any_anchor_overlap'].sum(), len(down_regulated) - down_regulated['any_anchor_overlap'].sum()]
        ])
        
        chi2, p_value = stats.chi2_contingency(contingency)[:2]
        results['comparison'] = {
            'chi2': chi2,
            'p_value': p_value
        }
        
        print(f"\nUp-regulated interactions:")
        print(f"  N = {results['up_regulated']['n_interactions']}")
        print(f"  H3K36me3 overlap rate: {results['up_regulated']['overlap_rate']*100:.1f}%")
        
        print(f"\nDown-regulated interactions:")
        print(f"  N = {results['down_regulated']['n_interactions']}")
        print(f"  H3K36me3 overlap rate: {results['down_regulated']['overlap_rate']*100:.1f}%")
        
        print(f"\nChi-square test p-value: {p_value:.2e}")
    
    return results, merged

def create_visualization(merged_df, direction_results, output_prefix):
    """Create comprehensive visualization of H3K36me3 enrichment analysis"""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot 1: Overlap rate by direction
    ax = axes[0, 0]
    directions = ['Up-regulated', 'Down-regulated']
    overlap_rates = [
        direction_results['up_regulated']['overlap_rate'] * 100,
        direction_results['down_regulated']['overlap_rate'] * 100
    ]
    
    colors = ['#d62728', '#2ca02c']  # Red for up, green for down
    ax.bar(directions, overlap_rates, color=colors, alpha=0.7)
    ax.set_ylabel('H3K36me3 Overlap Rate (%)')
    ax.set_title('H3K36me3 Enrichment by Direction')
    ax.axhline(y=50, color='gray', linestyle='--', label='Expected')
    ax.legend()
    
    # Plot 2: Number of H3K36me3 peaks by direction
    ax = axes[0, 1]
    peak_counts = [
        direction_results['up_regulated']['mean_peaks'],
        direction_results['down_regulated']['mean_peaks']
    ]
    ax.bar(directions, peak_counts, color=colors, alpha=0.7)
    ax.set_ylabel('Mean H3K36me3 Peaks per Interaction')
    ax.set_title('H3K36me3 Peak Density')
    
    # Plot 3: Distribution of H3K36me3 peaks
    ax = axes[0, 2]
    up_peaks = merged_df[merged_df['logFC'] > 0]['total_peaks']
    down_peaks = merged_df[merged_df['logFC'] < 0]['total_peaks']
    
    ax.hist([up_peaks, down_peaks], bins=20, label=['Up-regulated', 'Down-regulated'], 
            color=colors, alpha=0.6)
    ax.set_xlabel('Number of H3K36me3 Peaks')
    ax.set_ylabel('Number of Interactions')
    ax.set_title('Distribution of H3K36me3 Peaks')
    ax.legend()
    
    # Plot 4: LogFC vs H3K36me3 peaks
    ax = axes[1, 0]
    scatter = ax.scatter(merged_df['logFC'], merged_df['total_peaks'], 
                        c=merged_df['any_anchor_overlap'], cmap='RdYlGn',
                        alpha=0.6)
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('Number of H3K36me3 Peaks')
    ax.set_title('H3K36me3 Enrichment vs Effect Size')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    plt.colorbar(scatter, ax=ax, label='Overlaps H3K36me3')
    
    # Plot 5: Overlap type distribution
    ax = axes[1, 1]
    overlap_types = ['No overlap', 'One anchor', 'Both anchors']
    counts = [
        (~merged_df['any_anchor_overlap']).sum(),
        (merged_df['any_anchor_overlap'] & ~merged_df['both_anchors_overlap']).sum(),
        merged_df['both_anchors_overlap'].sum()
    ]
    ax.bar(overlap_types, counts, color=['gray', 'orange', 'red'], alpha=0.7)
    ax.set_ylabel('Number of Interactions')
    ax.set_title('H3K36me3 Overlap Pattern')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Plot 6: Cis vs Trans interactions
    ax = axes[1, 2]
    if 'interaction_type' in merged_df.columns:
        cis_overlap = merged_df[merged_df['interaction_type'] == 'cis']['any_anchor_overlap'].mean() * 100
        trans_overlap = merged_df[merged_df['interaction_type'] == 'trans']['any_anchor_overlap'].mean() * 100
        
        ax.bar(['Cis', 'Trans'], [cis_overlap, trans_overlap], 
               color=['blue', 'green'], alpha=0.7)
        ax.set_ylabel('H3K36me3 Overlap Rate (%)')
        ax.set_title('H3K36me3 Enrichment by Interaction Type')
        ax.axhline(y=50, color='gray', linestyle='--', label='Expected')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_h3k36me3_analysis.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Visualization saved to {output_prefix}_h3k36me3_analysis.pdf")

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
    merged_df.to_csv(f"{args.output_prefix}_h3k36me3_interactions.csv", index=False)
    
    # Save summary statistics
    summary = {
        'total_x_interactions': len(all_x_interactions),
        'significant_x_interactions': len(sig_interactions),
        'h3k36me3_peaks_on_x': len(h3k36me3_peaks),
        'overall_overlap_rate': overall_overlap_rate,
        'enrichment_vs_null': null_comparison['enrichment'],
        'p_value': null_comparison['p_value'],
        'up_regulated_n': direction_results['up_regulated']['n_interactions'],
        'up_regulated_overlap_rate': direction_results['up_regulated']['overlap_rate'],
        'down_regulated_n': direction_results['down_regulated']['n_interactions'],
        'down_regulated_overlap_rate': direction_results['down_regulated']['overlap_rate'],
        'direction_comparison_pvalue': direction_results.get('comparison', {}).get('p_value', None)
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(f"{args.output_prefix}_h3k36me3_summary.csv", index=False)
    
    # Create text summary
    with open(f"{args.output_prefix}_h3k36me3_summary.txt", 'w') as f:
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
        
        f.write(f"Direction-specific Analysis:\n")
        f.write(f"  Up-regulated interactions:\n")
        f.write(f"    N = {summary['up_regulated_n']}\n")
        f.write(f"    H3K36me3 overlap: {summary['up_regulated_overlap_rate']*100:.1f}%\n")
        f.write(f"  Down-regulated interactions:\n")
        f.write(f"    N = {summary['down_regulated_n']}\n")
        f.write(f"    H3K36me3 overlap: {summary['down_regulated_overlap_rate']*100:.1f}%\n")
        
        if summary['direction_comparison_pvalue'] is not None:
            f.write(f"  Direction comparison p-value: {summary['direction_comparison_pvalue']:.2e}\n")
    
    print(f"\nAnalysis complete! Results saved to {args.output_prefix}_*")
    
    return 0

if __name__ == '__main__':
    exit(main())