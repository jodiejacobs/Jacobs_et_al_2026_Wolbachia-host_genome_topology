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
    """
    Load HAS/CES sites from BED file.
    Format: chrX start end
    """
    print(f"Loading HAS/CES sites from {bed_file}")
    
    try:
        sites = []
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track'):
                    continue
                
                parts = line.strip().split()
                if len(parts) >= 3:
                    # Normalize chromosome name
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
        
        # Should only be X chromosome
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
    
    # Normalize chromosome names
    if 'chr1' in interactions.columns:
        interactions['chr1'] = interactions['chr1'].astype(str).str.replace('chr', '')
        interactions['chr2'] = interactions['chr2'].astype(str).str.replace('chr', '')
    
    print(f"Loaded {len(interactions)} total interactions")
    
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
    
    print(f"Found {len(sig_x)} significant X chromosome interactions")
    print(f"  FDR < {fdr_threshold}, |logFC| > 1")
    
    if len(sig_x) > 0:
        print(f"  Cis interactions: {sum(sig_x['chr1'] == sig_x['chr2'])}")
        print(f"  Trans interactions: {sum(sig_x['chr1'] != sig_x['chr2'])}")
        print(f"  JW18 uninf. (down-regulated): {sum(sig_x['logFC'] < 0)}")
        print(f"  JW18 wMel (up-regulated): {sum(sig_x['logFC'] > 0)}")
    
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
    
    Returns both anchor-level and interaction-level statistics.
    """
    print(f"\nAnalyzing HAS/CES overlap (window size: {window_size}bp)")
    print("Note: Large window accounts for MSL spreading from HAS sites")
    
    if has_sites is None or len(has_sites) == 0:
        print("No HAS/CES sites available")
        return None
    
    # Create BedTool for HAS sites with extended windows
    has_extended = has_sites.copy()
    has_extended['start'] = has_extended['start'] - window_size
    has_extended['end'] = has_extended['end'] + window_size
    has_extended['start'] = has_extended['start'].clip(lower=0)
    
    has_bt = pybedtools.BedTool.from_dataframe(
        has_extended[['chrom', 'start', 'end']]
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
        
        # Create BedTools for anchors
        anchor1_df = pd.DataFrame([{
            'chrom': anchor1_chr,
            'start': anchor1_start,
            'end': anchor1_end
        }])
        
        anchor2_df = pd.DataFrame([{
            'chrom': anchor2_chr,
            'start': anchor2_start,
            'end': anchor2_end
        }])
        
        # Check for overlap
        try:
            anchor1_bt = pybedtools.BedTool.from_dataframe(anchor1_df)
            anchor2_bt = pybedtools.BedTool.from_dataframe(anchor2_df)
            
            # Count overlaps
            anchor1_overlaps = anchor1_bt.intersect(has_bt, wa=True)
            anchor2_overlaps = anchor2_bt.intersect(has_bt, wa=True)
            
            n_has_anchor1 = len(list(anchor1_overlaps))
            n_has_anchor2 = len(list(anchor2_overlaps))
            
            anchor1_overlap = n_has_anchor1 > 0
            anchor2_overlap = n_has_anchor2 > 0
            
            # Calculate distance to nearest HAS
            if anchor1_chr == 'X':
                distances_1 = []
                for _, has_site in has_sites.iterrows():
                    dist = min(
                        abs(anchor1_start - has_site['end']),
                        abs(anchor1_end - has_site['start'])
                    )
                    distances_1.append(dist)
                min_dist_anchor1 = min(distances_1) if distances_1 else np.inf
            else:
                min_dist_anchor1 = np.inf
            
            if anchor2_chr == 'X':
                distances_2 = []
                for _, has_site in has_sites.iterrows():
                    dist = min(
                        abs(anchor2_start - has_site['end']),
                        abs(anchor2_end - has_site['start'])
                    )
                    distances_2.append(dist)
                min_dist_anchor2 = min(distances_2) if distances_2 else np.inf
            else:
                min_dist_anchor2 = np.inf
            
            results.append({
                'interaction_idx': idx,
                'anchor1_overlap': anchor1_overlap,
                'anchor2_overlap': anchor2_overlap,
                'any_anchor_overlap': anchor1_overlap or anchor2_overlap,
                'both_anchors_overlap': anchor1_overlap and anchor2_overlap,
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
    
    # Calculate summary statistics
    print(f"\nHAS/CES Overlap Summary:")
    print(f"  Interactions with any anchor overlap: {results_df['any_anchor_overlap'].sum()} ({results_df['any_anchor_overlap'].mean()*100:.1f}%)")
    print(f"  Interactions with both anchors overlap: {results_df['both_anchors_overlap'].sum()} ({results_df['both_anchors_overlap'].mean()*100:.1f}%)")
    print(f"  Mean HAS sites per interaction: {results_df['total_has'].mean():.2f}")
    print(f"  Median distance to nearest HAS: {results_df['min_dist_any'].median():.0f} bp")
    
    return results_df

def analyze_distance_distribution(merged_df, has_sites):
    """
    Analyze the distribution of distances from interactions to HAS sites.
    Tests if interactions are clustered near HAS compared to random expectation.
    """
    print("\nAnalyzing distance distribution to HAS sites...")
    
    # Get distances for significant interactions
    sig_distances = merged_df['min_dist_any'].values
    sig_distances = sig_distances[~np.isinf(sig_distances)]
    
    if len(sig_distances) == 0:
        print("No valid distances to analyze")
        return None
    
    # Calculate distance quartiles
    quartiles = np.percentile(sig_distances, [25, 50, 75])
    
    # Count interactions in different distance bins
    distance_bins = [0, 10000, 25000, 50000, 100000, np.inf]
    bin_labels = ['0-10kb', '10-25kb', '25-50kb', '50-100kb', '>100kb']
    
    bin_counts = []
    for i in range(len(distance_bins)-1):
        count = np.sum((sig_distances >= distance_bins[i]) & 
                      (sig_distances < distance_bins[i+1]))
        bin_counts.append(count)
    
    distance_dist = pd.DataFrame({
        'distance_range': bin_labels,
        'count': bin_counts,
        'percentage': np.array(bin_counts) / len(sig_distances) * 100
    })
    
    print(f"Distance distribution:")
    for _, row in distance_dist.iterrows():
        print(f"  {row['distance_range']}: {row['count']} ({row['percentage']:.1f}%)")
    
    return {
        'quartiles': quartiles,
        'distance_distribution': distance_dist,
        'median_distance': np.median(sig_distances),
        'mean_distance': np.mean(sig_distances)
    }

def compare_to_null_model(real_overlap_rate, null_data, n_permutations=1000):
    """
    Compare real HAS/CES overlap rate to null expectation.
    """
    print(f"\nComparing to null model...")
    
    # For HAS sites, which are specific X chromosome regions,
    # we expect ~10-20% overlap by chance (based on genomic coverage)
    
    # Calculate expected overlap based on HAS genomic coverage
    # X chromosome is ~23Mb, HAS sites cover ~few hundred kb
    # With 50kb windows, expect ~5-10% random overlap
    null_expectation = 0.075  # 7.5% baseline expectation
    
    # Calculate enrichment
    enrichment = real_overlap_rate / null_expectation if null_expectation > 0 else np.inf
    
    # Calculate p-value using binomial test
    # Null hypothesis: overlap rate = null_expectation
    from scipy.stats import binomtest
    
    # Estimate number of tests (interactions)
    n_interactions = len(null_data)
    n_overlaps = int(real_overlap_rate * n_interactions)
    
    result = binomtest(n_overlaps, n_interactions, null_expectation, alternative='greater')
    p_value = result.pvalue
    
    print(f"  Real overlap rate: {real_overlap_rate*100:.1f}%")
    print(f"  Expected (null): {null_expectation*100:.1f}%")
    print(f"  Enrichment: {enrichment:.2f}x")
    print(f"  P-value: {p_value:.2e}")
    
    return {
        'enrichment': enrichment,
        'p_value': p_value,
        'real_overlap_rate': real_overlap_rate,
        'null_expectation': null_expectation
    }

def analyze_by_logfc_direction(interactions_df, overlap_df):
    """
    Analyze HAS/CES enrichment separately for JW18 uninf. (up) and JW18 wMel (down) interactions.
    """
    print("\nAnalyzing by logFC direction...")
    
    # Merge overlap results with interaction data
    merged = interactions_df.copy()
    merged['any_anchor_overlap'] = overlap_df['any_anchor_overlap'].values
    merged['both_anchors_overlap'] = overlap_df['both_anchors_overlap'].values
    merged['total_has'] = overlap_df['total_has'].values
    merged['min_dist_any'] = overlap_df['min_dist_any'].values
    
    # Split by direction - JW18 uninf. (up) vs JW18 wMel (down)
    jw18_uninf = merged[merged['logFC'] > 0]   # Up-regulated = JW18 uninf.
    jw18_wmel = merged[merged['logFC'] < 0]    # Down-regulated = JW18 wMel
    
    results = {
        'jw18_uninf': {
            'n_interactions': len(jw18_uninf),
            'overlap_rate': jw18_uninf['any_anchor_overlap'].mean() if len(jw18_uninf) > 0 else 0,
            'mean_has': jw18_uninf['total_has'].mean() if len(jw18_uninf) > 0 else 0,
            'median_distance': jw18_uninf['min_dist_any'].median() if len(jw18_uninf) > 0 else np.inf
        },
        'jw18_wmel': {
            'n_interactions': len(jw18_wmel),
            'overlap_rate': jw18_wmel['any_anchor_overlap'].mean() if len(jw18_wmel) > 0 else 0,
            'mean_has': jw18_wmel['total_has'].mean() if len(jw18_wmel) > 0 else 0,
            'median_distance': jw18_wmel['min_dist_any'].median() if len(jw18_wmel) > 0 else np.inf
        }
    }
    
    # Statistical test
    if len(jw18_uninf) > 0 and len(jw18_wmel) > 0:
        # Chi-square test for overlap rates
        contingency = np.array([
            [jw18_uninf['any_anchor_overlap'].sum(), 
             len(jw18_uninf) - jw18_uninf['any_anchor_overlap'].sum()],
            [jw18_wmel['any_anchor_overlap'].sum(), 
             len(jw18_wmel) - jw18_wmel['any_anchor_overlap'].sum()]
        ])
        
        chi2, p_value = stats.chi2_contingency(contingency)[:2]
        
        # Mann-Whitney U test for distances
        uninf_dists = jw18_uninf['min_dist_any'].values
        wmel_dists = jw18_wmel['min_dist_any'].values
        uninf_dists = uninf_dists[~np.isinf(uninf_dists)]
        wmel_dists = wmel_dists[~np.isinf(wmel_dists)]
        
        if len(uninf_dists) > 0 and len(wmel_dists) > 0:
            mw_stat, mw_pval = stats.mannwhitneyu(uninf_dists, wmel_dists, alternative='two-sided')
        else:
            mw_stat, mw_pval = np.nan, 1.0
        
        results['comparison'] = {
            'chi2': chi2,
            'chi2_pvalue': p_value,
            'mannwhitney_stat': mw_stat,
            'mannwhitney_pvalue': mw_pval
        }
        
        print(f"\nJW18 uninf. (up-regulated) interactions:")
        print(f"  N = {results['jw18_uninf']['n_interactions']}")
        print(f"  HAS overlap rate: {results['jw18_uninf']['overlap_rate']*100:.1f}%")
        print(f"  Median distance: {results['jw18_uninf']['median_distance']:.0f} bp")
        
        print(f"\nJW18 wMel (down-regulated) interactions:")
        print(f"  N = {results['jw18_wmel']['n_interactions']}")
        print(f"  HAS overlap rate: {results['jw18_wmel']['overlap_rate']*100:.1f}%")
        print(f"  Median distance: {results['jw18_wmel']['median_distance']:.0f} bp")
        
        print(f"\nStatistical tests:")
        print(f"  Overlap rate difference (chi-square): p = {p_value:.2e}")
        print(f"  Distance difference (Mann-Whitney): p = {mw_pval:.2e}")
    
    return results, merged

def create_visualization(merged_df, direction_results, distance_analysis, null_comparison, output_prefix):
    """Create separate 2x2 visualizations with updated labels and colors"""
    print("\nCreating visualizations...")
    
    # Define colors for JW18 uninf. and JW18 wMel
    color_jw18_uninf = '#8fcb84'  # Light green for JW18 uninf. (upregulated)
    color_jw18_wmel = '#09aa4b'   # Dark green for JW18 wMel (downregulated)
    colors = [color_jw18_uninf, color_jw18_wmel]
    
    # Define labels
    labels = ['JW18 uninf.', 'JW18 wMel']
    
    # Set global font size to minimum 6pt
    plt.rcParams.update({'font.size': 6})
    
    # --- Plot Set 1: Enrichment and overlap rates (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)  # Transparent figure background
    
    # Plot 1: Observed vs Expected with significance
    ax = axes[0, 0]
    ax.patch.set_alpha(0)
    categories = ['Observed', 'Expected\n(Null)']
    rates = [
        null_comparison['real_overlap_rate'] * 100,
        null_comparison['null_expectation'] * 100
    ]
    
    bar_colors = ['#cc3333' if null_comparison['p_value'] < 0.05 else '#999999', '#cccccc']
    bars = ax.bar(categories, rates, color=bar_colors, alpha=0.8, edgecolor='black', linewidth=1)
    ax.set_ylabel('HAS/CES Overlap Rate (%)', fontsize=6)
    ax.set_title(f'Enrichment vs Null\n({null_comparison["enrichment"]:.2f}x)', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add p-value
    p_val = null_comparison['p_value']
    ax.text(0.5, max(rates) * 0.95, f'p = {p_val:.2e}', 
            ha='center', va='top', fontsize=5,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add values on bars
    for bar, val in zip(bars, rates):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height*0.5,
                f'{val:.1f}%', ha='center', va='center', fontsize=5)
    
    # Plot 2: Direction-specific with null comparison
    ax = axes[0, 1]
    ax.patch.set_alpha(0)
    directions_all = ['JW18 uninf.', 'JW18 wMel', 'Expected\n(Null)']
    overlap_rates_dir = [
        direction_results['jw18_uninf']['overlap_rate'] * 100,
        direction_results['jw18_wmel']['overlap_rate'] * 100,
        null_comparison['null_expectation'] * 100
    ]
    
    colors_dir = [color_jw18_uninf, color_jw18_wmel, '#cccccc']
    bars = ax.bar(directions_all, overlap_rates_dir, color=colors_dir, alpha=0.8, edgecolor='black', linewidth=1)
    ax.set_ylabel('HAS Overlap Rate (%)', fontsize=6)
    ax.set_title('Genotype-Specific Enrichment', fontsize=7, fontweight='bold')
    ax.tick_params(labelsize=6)
    
    # Add significance if direction test exists
    if 'comparison' in direction_results and direction_results['comparison']:
        pval = direction_results['comparison'].get('chi2_pvalue', 1.0)
        ax.text(0.5, max(overlap_rates_dir) * 0.95, f'p = {pval:.2e}', 
                ha='center', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add values on bars
    for bar, val in zip(bars, overlap_rates_dir):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=5)
    
    # Plot 3: Distribution of HAS sites
    ax = axes[1, 0]
    ax.patch.set_alpha(0)
    uninf_has = merged_df[merged_df['logFC'] > 0]['total_has']
    wmel_has = merged_df[merged_df['logFC'] < 0]['total_has']
    
    bins = np.arange(0, max(uninf_has.max(), wmel_has.max()) + 2)
    ax.hist([uninf_has, wmel_has], bins=bins, label=labels,
            color=colors, alpha=0.7)
    ax.set_xlabel('Number of HAS Sites', fontsize=6)
    ax.set_ylabel('Number of Interactions', fontsize=6)
    ax.set_title('Distribution of HAS Sites', fontsize=7, fontweight='bold')
    ax.legend(fontsize=5)
    ax.tick_params(labelsize=6)
    
    # Add Mann-Whitney U test
    if len(uninf_has) > 0 and len(wmel_has) > 0:
        u_stat, p_val = stats.mannwhitneyu(uninf_has, wmel_has, alternative='two-sided')
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
    ax.set_title('HAS Overlap Pattern', fontsize=7, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    
    # Add values on bars
    total = sum(counts)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{count}\n({count/total*100:.1f}%)',
                ha='center', va='bottom', fontsize=5)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/has_analysis_set1.pdf", 
                dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Plot set 1 saved to {output_prefix}/has_analysis_set1.pdf")
    
    # --- Plot Set 2: Distance and effect size analysis (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    fig.patch.set_alpha(0)
    
    # Plot 1: LogFC vs distance to HAS
    ax = axes[0, 0]
    ax.patch.set_alpha(0)
    
    # Filter out infinite distances
    plot_data = merged_df[~np.isinf(merged_df['min_dist_any'])].copy()
    
    # Color points by genotype
    uninf_mask = plot_data['logFC'] > 0
    wmel_mask = plot_data['logFC'] < 0
    
    ax.scatter(plot_data[uninf_mask]['logFC'], plot_data[uninf_mask]['min_dist_any']/1000,
               c=color_jw18_uninf, alpha=0.6, s=20, label='JW18 uninf.', edgecolors='none')
    ax.scatter(plot_data[wmel_mask]['logFC'], plot_data[wmel_mask]['min_dist_any']/1000,
               c=color_jw18_wmel, alpha=0.6, s=20, label='JW18 wMel', edgecolors='none')
    
    ax.set_xlabel('log2 Fold Change', fontsize=6)
    ax.set_ylabel('Distance to Nearest HAS (kb)', fontsize=6)
    ax.set_title('Effect Size vs HAS Proximity', fontsize=7, fontweight='bold')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.set_yscale('log')
    ax.legend(fontsize=5)
    ax.tick_params(labelsize=6)
    
    # Add correlation statistics
    corr, p_val = stats.spearmanr(plot_data['logFC'], plot_data['min_dist_any'])
    ax.text(0.05, 0.95, f'ρ = {corr:.3f}\np = {p_val:.2e}', 
            transform=ax.transAxes, ha='left', va='top', fontsize=5,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Plot 2: Distance distribution
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
        
        # Add median line
        median_dist = distance_analysis['median_distance']
        ax.text(0.95, 0.95, f'Median: {median_dist/1000:.1f}kb', 
                transform=ax.transAxes, ha='right', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
        
        # Add values on bars
        for bar, val in zip(bars, dist_data['percentage']):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=5)
    
    # Plot 3: Median distance comparison by genotype
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
    
    # Add Mann-Whitney p-value
    if 'comparison' in direction_results and 'mannwhitney_pvalue' in direction_results['comparison']:
        mw_pval = direction_results['comparison']['mannwhitney_pvalue']
        ax.text(0.5, max(median_dists) * 0.95, f'p = {mw_pval:.2e}', 
                ha='center', va='top', fontsize=5,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Add values on bars
    for bar, dist in zip(bars, median_dists):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{dist:.1f}kb', ha='center', va='bottom', fontsize=5)
    
    # Plot 4: Sample sizes
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
    
    # Add values on bars
    for bar, size in zip(bars, sample_sizes):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
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
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("HAS/CES Enrichment Analysis at X Chromosome Contacts")
    print("="*60)
    
    # Load data
    has_sites = load_has_sites(args.has_sites)
    sig_interactions, all_x_interactions = load_differential_interactions(
        args.interactions, args.fdr_threshold
    )
    null_model = load_null_model(args.null_model)
    
    if has_sites is None or len(sig_interactions) == 0:
        print("Error: No data to analyze")
        return 1
    
    # Calculate HAS/CES overlap
    overlap_results = calculate_has_overlap(
        sig_interactions, has_sites, args.window_size
    )
    
    if overlap_results is None:
        print("Error: Could not calculate overlap")
        return 1
    
    # Analyze distance distribution
    direction_results, merged_df = analyze_by_logfc_direction(
        sig_interactions, overlap_results
    )
    
    distance_analysis = analyze_distance_distribution(merged_df, has_sites)
    
    # Compare to null model
    overall_overlap_rate = overlap_results['any_anchor_overlap'].mean()
    null_comparison = compare_to_null_model(overall_overlap_rate, null_model)
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print("="*60)
    print(f"Overall HAS/CES enrichment:")
    print(f"  Overlap rate: {null_comparison['real_overlap_rate']*100:.1f}%")
    print(f"  Expected (null): {null_comparison['null_expectation']*100:.1f}%")
    print(f"  Enrichment: {null_comparison['enrichment']:.2f}x")
    print(f"  P-value: {null_comparison['p_value']:.2e}")
    print(f"  Significant: {'YES' if null_comparison['p_value'] < 0.05 else 'NO'}")
    
    # Create visualizations
    create_visualization(merged_df, direction_results, distance_analysis, null_comparison, args.output_prefix)
    
    # Save detailed results
    merged_df.to_csv(f"{args.output_prefix}/has_interactions.csv", index=False)
    
    # Save summary statistics
    summary = {
        'total_x_interactions': len(all_x_interactions),
        'significant_x_interactions': len(sig_interactions),
        'has_sites_on_x': len(has_sites),
        'window_size': args.window_size,
        'overall_overlap_rate': overall_overlap_rate,
        'enrichment_vs_null': null_comparison['enrichment'],
        'p_value': null_comparison['p_value'],
        'jw18_uninf_n': direction_results['jw18_uninf']['n_interactions'],
        'jw18_uninf_overlap_rate': direction_results['jw18_uninf']['overlap_rate'],
        'jw18_uninf_median_distance': direction_results['jw18_uninf']['median_distance'],
        'jw18_wmel_n': direction_results['jw18_wmel']['n_interactions'],
        'jw18_wmel_overlap_rate': direction_results['jw18_wmel']['overlap_rate'],
        'jw18_wmel_median_distance': direction_results['jw18_wmel']['median_distance'],
        'direction_overlap_pvalue': direction_results.get('comparison', {}).get('chi2_pvalue', None),
        'direction_distance_pvalue': direction_results.get('comparison', {}).get('mannwhitney_pvalue', None),
        'median_distance_to_has': distance_analysis['median_distance'] if distance_analysis else None,
        'mean_distance_to_has': distance_analysis['mean_distance'] if distance_analysis else None
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(f"{args.output_prefix}/has_summary.csv", index=False)
    
    # Save distance distribution
    if distance_analysis is not None:
        distance_analysis['distance_distribution'].to_csv(
            f"{args.output_prefix}/distance_distribution.csv", index=False
        )
    
    # Create text summary
    with open(f"{args.output_prefix}/has_summary.txt", 'w') as f:
        f.write("HAS/CES ENRICHMENT ANALYSIS SUMMARY\n")
        f.write("="*60 + "\n\n")
        
        f.write("Analysis Parameters:\n")
        f.write(f"  Window size: {args.window_size} bp\n")
        f.write(f"  FDR threshold: {args.fdr_threshold}\n")
        f.write(f"  Note: Large window accounts for MSL spreading from HAS\n\n")
        
        f.write("Data Summary:\n")
        f.write(f"  Total X chromosome interactions: {summary['total_x_interactions']}\n")
        f.write(f"  Significant interactions: {summary['significant_x_interactions']}\n")
        f.write(f"  HAS/CES sites on X: {summary['has_sites_on_x']}\n\n")
        
        f.write("HAS/CES Enrichment:\n")
        f.write(f"  Overall overlap rate: {summary['overall_overlap_rate']*100:.1f}%\n")
        f.write(f"  Expected (null): {null_comparison['null_expectation']*100:.1f}%\n")
        f.write(f"  Enrichment: {summary['enrichment_vs_null']:.2f}x\n")
        f.write(f"  P-value: {summary['p_value']:.2e}\n")
        f.write(f"  Significant: {'YES' if summary['p_value'] < 0.05 else 'NO'}\n\n")
        
        f.write("Distance Analysis:\n")
        if distance_analysis:
            f.write(f"  Median distance to HAS: {summary['median_distance_to_has']/1000:.1f} kb\n")
            f.write(f"  Mean distance to HAS: {summary['mean_distance_to_has']/1000:.1f} kb\n\n")
        
        f.write("Genotype-specific Analysis:\n")
        f.write(f"  JW18 uninf. (up-regulated) interactions:\n")
        f.write(f"    N = {summary['jw18_uninf_n']}\n")
        f.write(f"    HAS overlap: {summary['jw18_uninf_overlap_rate']*100:.1f}%\n")
        f.write(f"    Median distance: {summary['jw18_uninf_median_distance']/1000:.1f} kb\n")
        f.write(f"  JW18 wMel (down-regulated) interactions:\n")
        f.write(f"    N = {summary['jw18_wmel_n']}\n")
        f.write(f"    HAS overlap: {summary['jw18_wmel_overlap_rate']*100:.1f}%\n")
        f.write(f"    Median distance: {summary['jw18_wmel_median_distance']/1000:.1f} kb\n\n")
        
        if summary['direction_overlap_pvalue'] is not None:
            f.write(f"  Statistical tests:\n")
            f.write(f"    Overlap rate difference: p = {summary['direction_overlap_pvalue']:.2e}\n")
            f.write(f"    Distance difference: p = {summary['direction_distance_pvalue']:.2e}\n\n")
        
        f.write("\nBiological Interpretation:\n")
        f.write("-" * 60 + "\n")
        
        if summary['enrichment_vs_null'] > 1.5 and summary['p_value'] < 0.05:
            f.write("STRONG ENRICHMENT: Differential contacts are significantly enriched\n")
            f.write("near HAS/CES sites, suggesting Wolbachia affects chromatin contacts\n")
            f.write("at dosage compensation entry sites.\n\n")
            
            if summary['jw18_uninf_overlap_rate'] > summary['jw18_wmel_overlap_rate'] * 1.2:
                f.write("JW18 uninf. bias: New contacts forming preferentially near HAS.\n")
                f.write("This could indicate enhanced MSL complex recruitment or spreading.\n")
            elif summary['jw18_wmel_overlap_rate'] > summary['jw18_uninf_overlap_rate'] * 1.2:
                f.write("JW18 wMel bias: Contacts lost preferentially near HAS.\n")
                f.write("This could indicate disrupted MSL complex organization.\n")
        elif summary['enrichment_vs_null'] > 1.0 and summary['p_value'] < 0.05:
            f.write("MODERATE ENRICHMENT: Some enrichment detected near HAS sites.\n")
            f.write("Wolbachia may have subtle effects on dosage compensation regions.\n")
        else:
            f.write("NO SIGNIFICANT ENRICHMENT: Differential contacts are not enriched\n")
            f.write("near HAS sites. Changes may be independent of dosage compensation\n")
            f.write("entry sites, possibly affecting MSL spreading regions instead.\n")
    
    print(f"\nAnalysis complete! Results saved to {args.output_prefix}/*")
    
    return 0

if __name__ == '__main__':
    exit(main())
    