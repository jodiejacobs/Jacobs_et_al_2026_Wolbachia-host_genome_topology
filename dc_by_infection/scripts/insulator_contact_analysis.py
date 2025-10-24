#!/usr/bin/env python3
"""
Analyze insulator protein enrichment at differential chromatin contacts.
Uses permutation testing to assess statistical significance.

OPTIMIZED VERSION - Uses vectorized BedTools operations instead of per-interaction loops
FIXED VERSION - Properly sorts BED data before bedtools operations

Insulators tested:
- Class I (CTCF-dependent)
- Class II (CTCF-independent: Su(Hw), BEAF-32, etc.)

Based on:
- Schwartz et al. (2012) - Insulator classification
- Emberly et al. (2008) - CTCF insulator function
- Wood et al. (2011) - DCC and insulators
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
import warnings
warnings.filterwarnings('ignore')

def load_insulator_sites(class_i_file, class_ii_file):
    """
    Load Class I and Class II insulator sites.
    Class I: CTCF-dependent
    Class II: CTCF-independent (Su(Hw), BEAF-32, etc.)
    """
    print("Loading insulator sites...")
    
    insulators = {}
    
    # Load Class I (CTCF-dependent)
    try:
        class_i = []
        with open(class_i_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    class_i.append({
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'name': parts[3] if len(parts) > 3 else 'Class_I',
                        'class': 'I'
                    })
        
        class_i_df = pd.DataFrame(class_i)
        insulators['Class_I'] = class_i_df
        print(f"Class I (CTCF-dependent): {len(class_i_df)} sites")
        print(f"  Chromosomes: {sorted(class_i_df['chrom'].unique())}")
        
    except Exception as e:
        print(f"Error loading Class I insulators: {e}")
        insulators['Class_I'] = pd.DataFrame()
    
    # Load Class II (CTCF-independent)
    try:
        class_ii = []
        with open(class_ii_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    class_ii.append({
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'name': parts[3] if len(parts) > 3 else 'Class_II',
                        'class': 'II'
                    })
        
        class_ii_df = pd.DataFrame(class_ii)
        insulators['Class_II'] = class_ii_df
        print(f"Class II (CTCF-independent): {len(class_ii_df)} sites")
        print(f"  Chromosomes: {sorted(class_ii_df['chrom'].unique())}")
        
    except Exception as e:
        print(f"Error loading Class II insulators: {e}")
        insulators['Class_II'] = pd.DataFrame()
    
    # Combine for total count
    all_insulators = pd.concat([class_i_df, class_ii_df], ignore_index=True)
    insulators['All'] = all_insulators
    print(f"Total insulators: {len(all_insulators)} sites")
    
    return insulators

def load_differential_interactions(interactions_file, fdr_threshold=0.05):
    """Load significant differential interactions"""
    print(f"\nLoading differential interactions from {interactions_file}")
    
    interactions = pd.read_csv(interactions_file)
    
    print(f"Loaded {len(interactions)} total interactions")
    
    # Filter for significant
    sig_interactions = interactions[
        (interactions['FDR'] < fdr_threshold) & 
        (abs(interactions['logFC']) > 1)
    ].copy()
    
    print(f"Found {len(sig_interactions)} significant interactions")
    print(f"  FDR < {fdr_threshold}, |logFC| > 1")
    
    if len(sig_interactions) > 0:
        print(f"  Cis interactions: {sum(sig_interactions['chr1'] == sig_interactions['chr2'])}")
        print(f"  Trans interactions: {sum(sig_interactions['chr1'] != sig_interactions['chr2'])}")
        print(f"  Up-regulated: {sum(sig_interactions['logFC'] > 0)}")
        print(f"  Down-regulated: {sum(sig_interactions['logFC'] < 0)}")
    
    return sig_interactions, interactions

def load_null_model(null_file):
    """Load null model interactions"""
    print(f"\nLoading null model from {null_file}")
    null_data = pd.read_csv(null_file)
    print(f"Loaded {len(null_data)} null model interactions")
    return null_data

def calculate_insulator_overlap(interactions_df, insulator_sites, window_size=10000):
    """
    OPTIMIZED: Calculate overlap between interaction anchors and insulator sites.
    Uses vectorized BedTools operations instead of per-interaction loops.
    
    FIXED: Properly sorts all BED data before bedtools operations to avoid sorting errors.
    """
    print(f"\nAnalyzing insulator overlap (window size: {window_size}bp)")
    
    if len(insulator_sites) == 0:
        print("No insulator sites available")
        return None
    
    # Create extended insulator regions
    insulators_extended = insulator_sites.copy()
    insulators_extended['start'] = (insulators_extended['start'] - window_size).clip(lower=0)
    insulators_extended['end'] = insulators_extended['end'] + window_size
    
    # CRITICAL FIX: Sort before creating BedTool
    insulators_extended = insulators_extended.sort_values(['chrom', 'start', 'end'])
    
    insulator_bt = pybedtools.BedTool.from_dataframe(
        insulators_extended[['chrom', 'start', 'end']]
    )
    
    # OPTIMIZATION 1: Create BedTools for all anchors at once instead of in a loop
    anchor1_df = interactions_df[['chr1', 'start1', 'end1']].copy()
    anchor1_df.columns = ['chrom', 'start', 'end']
    anchor1_df['idx'] = interactions_df.index
    
    anchor2_df = interactions_df[['chr2', 'start2', 'end2']].copy()
    anchor2_df.columns = ['chrom', 'start', 'end']
    anchor2_df['idx'] = interactions_df.index
    
    # CRITICAL FIX: Sort anchor dataframes before creating BedTool objects
    anchor1_df = anchor1_df.sort_values(['chrom', 'start', 'end'])
    anchor2_df = anchor2_df.sort_values(['chrom', 'start', 'end'])
    
    anchor1_bt = pybedtools.BedTool.from_dataframe(anchor1_df)
    anchor2_bt = pybedtools.BedTool.from_dataframe(anchor2_df)
    
    # OPTIMIZATION 2: Use BedTools intersect operations in batch
    # Count overlaps for each anchor
    anchor1_overlaps = anchor1_bt.intersect(insulator_bt, c=True)
    anchor2_overlaps = anchor2_bt.intersect(insulator_bt, c=True)
    
    # Convert to DataFrames
    anchor1_counts = pd.DataFrame([
        (x.chrom, x.start, x.end, int(x.name), int(x.fields[3]))
        for x in anchor1_overlaps
    ], columns=['chrom', 'start', 'end', 'idx', 'n_ins_anchor1'])
    
    anchor2_counts = pd.DataFrame([
        (x.chrom, x.start, x.end, int(x.name), int(x.fields[3]))
        for x in anchor2_overlaps
    ], columns=['chrom', 'start', 'end', 'idx', 'n_ins_anchor2'])
    
    # OPTIMIZATION 3: Use BedTools closest operation for distances
    # This is MUCH faster than manual distance calculation
    # CRITICAL FIX: Sort before creating BedTool
    insulator_sites_sorted = insulator_sites.sort_values(['chrom', 'start', 'end'])
    insulator_base_bt = pybedtools.BedTool.from_dataframe(
        insulator_sites_sorted[['chrom', 'start', 'end']]
    )
    
    anchor1_closest = anchor1_bt.closest(insulator_base_bt, d=True)
    anchor2_closest = anchor2_bt.closest(insulator_base_bt, d=True)
    
    # Extract distances
    anchor1_dists = {}
    for x in anchor1_closest:
        idx = int(x.fields[3])
        dist = int(x.fields[-1])  # Last field is distance
        if idx not in anchor1_dists or dist < anchor1_dists[idx]:
            anchor1_dists[idx] = dist
    
    anchor2_dists = {}
    for x in anchor2_closest:
        idx = int(x.fields[3])
        dist = int(x.fields[-1])
        if idx not in anchor2_dists or dist < anchor2_dists[idx]:
            anchor2_dists[idx] = dist
    
    # Merge results
    results = []
    for idx in interactions_df.index:
        n_ins_1 = anchor1_counts[anchor1_counts['idx'] == idx]['n_ins_anchor1'].values
        n_ins_1 = n_ins_1[0] if len(n_ins_1) > 0 else 0
        
        n_ins_2 = anchor2_counts[anchor2_counts['idx'] == idx]['n_ins_anchor2'].values
        n_ins_2 = n_ins_2[0] if len(n_ins_2) > 0 else 0
        
        min_dist_1 = anchor1_dists.get(idx, np.inf)
        min_dist_2 = anchor2_dists.get(idx, np.inf)
        
        # Handle cases where distance is -1 (no overlap found on same chromosome)
        if min_dist_1 == -1:
            min_dist_1 = np.inf
        if min_dist_2 == -1:
            min_dist_2 = np.inf
        
        results.append({
            'interaction_idx': idx,
            'anchor1_overlap': n_ins_1 > 0,
            'anchor2_overlap': n_ins_2 > 0,
            'any_anchor_overlap': (n_ins_1 > 0) or (n_ins_2 > 0),
            'both_anchors_overlap': (n_ins_1 > 0) and (n_ins_2 > 0),
            'n_ins_anchor1': n_ins_1,
            'n_ins_anchor2': n_ins_2,
            'total_insulators': n_ins_1 + n_ins_2,
            'min_dist_anchor1': min_dist_1,
            'min_dist_anchor2': min_dist_2,
            'min_dist_any': min(min_dist_1, min_dist_2)
        })
    
    results_df = pd.DataFrame(results)
    
    print(f"\nInsulator Overlap Summary:")
    print(f"  Interactions with any anchor overlap: {results_df['any_anchor_overlap'].sum()} ({results_df['any_anchor_overlap'].mean()*100:.1f}%)")
    print(f"  Interactions with both anchors overlap: {results_df['both_anchors_overlap'].sum()} ({results_df['both_anchors_overlap'].mean()*100:.1f}%)")
    print(f"  Mean insulators per interaction: {results_df['total_insulators'].mean():.2f}")
    
    return results_df

def permutation_test_enrichment(interactions_df, insulator_sites, genome_file, 
                                window_size=10000, n_permutations=1000):
    """
    Perform permutation test to assess significance of insulator enrichment.
    
    Null hypothesis: Random genomic regions show similar insulator overlap.
    """
    print(f"\nPerforming permutation test ({n_permutations} permutations)...")
    
    # Calculate observed overlap
    overlap_results = calculate_insulator_overlap(interactions_df, insulator_sites, window_size)
    observed_rate = overlap_results['any_anchor_overlap'].mean()
    
    print(f"Observed overlap rate: {observed_rate*100:.1f}%")
    
    # Create BedTool for interactions
    interaction_regions = []
    for _, row in interactions_df.iterrows():
        interaction_regions.append([row['chr1'], row['start1'], row['end1']])
        interaction_regions.append([row['chr2'], row['start2'], row['end2']])
    
    interactions_bed_df = pd.DataFrame(interaction_regions, columns=['chrom', 'start', 'end'])
    interactions_bed_df = interactions_bed_df.drop_duplicates()
    
    # CRITICAL FIX: Sort before creating BedTool
    interactions_bed_df = interactions_bed_df.sort_values(['chrom', 'start', 'end'])
    
    interactions_bt = pybedtools.BedTool.from_dataframe(interactions_bed_df)
    
    # Extend insulator sites
    insulators_extended = insulator_sites.copy()
    insulators_extended['start'] = (insulators_extended['start'] - window_size).clip(lower=0)
    insulators_extended['end'] = insulators_extended['end'] + window_size
    
    # CRITICAL FIX: Sort before creating BedTool
    insulators_extended = insulators_extended.sort_values(['chrom', 'start', 'end'])
    
    insulator_bt = pybedtools.BedTool.from_dataframe(
        insulators_extended[['chrom', 'start', 'end']]
    )
    
    # Perform permutations
    null_overlap_rates = []
    
    for i in range(n_permutations):
        try:
            # Shuffle interaction regions
            shuffled = interactions_bt.shuffle(g=genome_file, chrom=True, seed=i)
            
            # Calculate overlap with insulators
            overlaps = shuffled.intersect(insulator_bt, u=True)
            overlap_count = len(list(overlaps))
            total_count = len(interactions_bed_df)
            
            null_rate = overlap_count / total_count if total_count > 0 else 0
            null_overlap_rates.append(null_rate)
            
            if (i + 1) % 100 == 0:
                print(f"  Completed {i + 1}/{n_permutations} permutations")
                
        except Exception as e:
            if i < 10:
                print(f"  Warning: Permutation {i} failed: {e}")
            continue
    
    null_overlap_rates = np.array(null_overlap_rates)
    
    if len(null_overlap_rates) == 0:
        print("Error: No successful permutations")
        return None
    
    # Calculate statistics
    expected_rate = np.median(null_overlap_rates)
    null_mean = np.mean(null_overlap_rates)
    null_std = np.std(null_overlap_rates)
    
    # Calculate enrichment
    enrichment = observed_rate / expected_rate if expected_rate > 0 else np.inf
    
    # Calculate p-value (one-tailed test for enrichment)
    p_value = (np.sum(null_overlap_rates >= observed_rate) + 1) / (len(null_overlap_rates) + 1)
    
    # Calculate z-score
    z_score = (observed_rate - null_mean) / null_std if null_std > 0 else 0
    
    print(f"\nPermutation Test Results:")
    print(f"  Observed: {observed_rate*100:.1f}%")
    print(f"  Expected (median null): {expected_rate*100:.1f}%")
    print(f"  Null mean ± SD: {null_mean*100:.1f}% ± {null_std*100:.1f}%")
    print(f"  Enrichment: {enrichment:.2f}x")
    print(f"  Z-score: {z_score:.2f}")
    print(f"  P-value: {p_value:.4f}")
    print(f"  Significant: {'YES' if p_value < 0.05 else 'NO'}")
    
    return {
        'observed_rate': observed_rate,
        'expected_rate': expected_rate,
        'null_mean': null_mean,
        'null_std': null_std,
        'null_distribution': null_overlap_rates,
        'enrichment': enrichment,
        'z_score': z_score,
        'p_value': p_value,
        'n_permutations': len(null_overlap_rates)
    }

def analyze_by_insulator_class(interactions_df, insulators_dict, genome_file, 
                               window_size=10000, n_permutations=1000):
    """
    Analyze enrichment separately for Class I and Class II insulators.
    """
    print("\n" + "="*60)
    print("ANALYZING BY INSULATOR CLASS")
    print("="*60)
    
    results = {}
    
    for ins_class, ins_sites in insulators_dict.items():
        if len(ins_sites) == 0:
            print(f"\nSkipping {ins_class}: no sites available")
            continue
        
        print(f"\n{'='*60}")
        print(f"Analyzing {ins_class} insulators")
        print(f"{'='*60}")
        
        # Calculate overlap
        overlap_results = calculate_insulator_overlap(interactions_df, ins_sites, window_size)
        
        if overlap_results is None:
            continue
        
        # Permutation test
        perm_results = permutation_test_enrichment(
            interactions_df, ins_sites, genome_file, window_size, n_permutations
        )
        
        results[ins_class] = {
            'overlap': overlap_results,
            'permutation': perm_results
        }
    
    return results

def analyze_by_logfc_direction(interactions_df, overlap_df):
    """Analyze insulator enrichment by up vs down-regulation"""
    print("\nAnalyzing by logFC direction...")
    
    merged = interactions_df.copy()
    merged['any_anchor_overlap'] = overlap_df['any_anchor_overlap'].values
    merged['total_insulators'] = overlap_df['total_insulators'].values
    merged['min_dist_any'] = overlap_df['min_dist_any'].values
    
    up_regulated = merged[merged['logFC'] > 0]
    down_regulated = merged[merged['logFC'] < 0]
    
    results = {
        'up_regulated': {
            'n_interactions': len(up_regulated),
            'overlap_rate': up_regulated['any_anchor_overlap'].mean() if len(up_regulated) > 0 else 0,
            'mean_insulators': up_regulated['total_insulators'].mean() if len(up_regulated) > 0 else 0,
            'median_distance': up_regulated['min_dist_any'].median() if len(up_regulated) > 0 else np.inf
        },
        'down_regulated': {
            'n_interactions': len(down_regulated),
            'overlap_rate': down_regulated['any_anchor_overlap'].mean() if len(down_regulated) > 0 else 0,
            'mean_insulators': down_regulated['total_insulators'].mean() if len(down_regulated) > 0 else 0,
            'median_distance': down_regulated['min_dist_any'].median() if len(down_regulated) > 0 else np.inf
        }
    }
    
    # Statistical test
    if len(up_regulated) > 0 and len(down_regulated) > 0:
        contingency = np.array([
            [up_regulated['any_anchor_overlap'].sum(), 
             len(up_regulated) - up_regulated['any_anchor_overlap'].sum()],
            [down_regulated['any_anchor_overlap'].sum(), 
             len(down_regulated) - down_regulated['any_anchor_overlap'].sum()]
        ])
        
        chi2, p_value = stats.chi2_contingency(contingency)[:2]
        results['comparison'] = {'chi2': chi2, 'p_value': p_value}
        
        print(f"\nUp-regulated: {results['up_regulated']['overlap_rate']*100:.1f}%")
        print(f"Down-regulated: {results['down_regulated']['overlap_rate']*100:.1f}%")
        print(f"Chi-square p-value: {p_value:.4f}")
    
    return results, merged

def check_specific_interaction(interactions_df, overlap_df, chr1, start1, end1, 
                               chr2, start2, end2):
    """
    Check if a specific interaction (like the NDF example) has insulator enrichment.
    Example: 2L:10000000-10008000_X:22432000-22440000
    """
    print(f"\n{'='*60}")
    print("CHECKING SPECIFIC INTERACTION")
    print(f"{'='*60}")
    print(f"Looking for: {chr1}:{start1}-{end1}_{chr2}:{start2}-{end2}")
    
    # Find matching interaction
    matches = interactions_df[
        (interactions_df['chr1'] == chr1) &
        (interactions_df['start1'] >= start1 - 1000) &
        (interactions_df['end1'] <= end1 + 1000) &
        (interactions_df['chr2'] == chr2) &
        (interactions_df['start2'] >= start2 - 1000) &
        (interactions_df['end2'] <= end2 + 1000)
    ]
    
    if len(matches) == 0:
        print("Interaction not found in dataset")
        return None
    
    print(f"\nFound {len(matches)} matching interaction(s)")
    
    for _, interaction in matches.iterrows():
        idx = interaction.name
        overlap_info = overlap_df[overlap_df['interaction_idx'] == idx]
        
        if len(overlap_info) == 0:
            continue
        
        overlap_info = overlap_info.iloc[0]
        
        print(f"\nInteraction Details:")
        print(f"  Coordinates: {interaction['chr1']}:{interaction['start1']}-{interaction['end1']}_"
              f"{interaction['chr2']}:{interaction['start2']}-{interaction['end2']}")
        print(f"  logFC: {interaction['logFC']:.2f}")
        print(f"  FDR: {interaction['FDR']:.2e}")
        
        print(f"\nInsulator Enrichment:")
        print(f"  Anchor 1 overlaps insulator: {'YES' if overlap_info['anchor1_overlap'] else 'NO'}")
        print(f"  Anchor 2 overlaps insulator: {'YES' if overlap_info['anchor2_overlap'] else 'NO'}")
        print(f"  Number of insulators (anchor 1): {overlap_info['n_ins_anchor1']}")
        print(f"  Number of insulators (anchor 2): {overlap_info['n_ins_anchor2']}")
        print(f"  Distance to nearest (anchor 1): {overlap_info['min_dist_anchor1']:.0f} bp")
        print(f"  Distance to nearest (anchor 2): {overlap_info['min_dist_anchor2']:.0f} bp")
        
        return {
            'interaction': interaction,
            'overlap': overlap_info
        }
    
    return None

def create_visualization(class_results, merged_df, perm_results, output_prefix):
    """Create comprehensive visualization"""
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Plot 1: Enrichment by insulator class
    ax = fig.add_subplot(gs[0, 0])
    classes = []
    enrichments = []
    p_values = []
    
    for ins_class, results in class_results.items():
        if results and results['permutation']:
            classes.append(ins_class)
            enrichments.append(results['permutation']['enrichment'])
            p_values.append(results['permutation']['p_value'])
    
    colors = ['#1f77b4' if p < 0.05 else '#d3d3d3' for p in p_values]
    bars = ax.bar(classes, enrichments, color=colors, alpha=0.7)
    ax.axhline(y=1, color='gray', linestyle='--', label='No enrichment')
    ax.set_ylabel('Enrichment vs Null')
    ax.set_title('Insulator Enrichment by Class')
    ax.legend()
    
    # Add significance stars
    for bar, pval in zip(bars, p_values):
        stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                stars, ha='center', va='bottom', fontsize=12)
    
    # Plot 2: Null distribution with observed value
    ax = fig.add_subplot(gs[0, 1])
    if perm_results and perm_results['null_distribution'] is not None:
        null_dist = perm_results['null_distribution'] * 100
        ax.hist(null_dist, bins=30, alpha=0.7, color='gray', edgecolor='black')
        ax.axvline(perm_results['observed_rate']*100, color='red', 
                  linestyle='--', linewidth=2, label='Observed')
        ax.axvline(perm_results['expected_rate']*100, color='blue',
                  linestyle='--', linewidth=2, label='Expected (median)')
        ax.set_xlabel('Overlap Rate (%)')
        ax.set_ylabel('Frequency')
        ax.set_title('Permutation Test Null Distribution')
        ax.legend()
    
    # Plot 3: Overlap rate by direction
    ax = fig.add_subplot(gs[0, 2])
    if 'logFC' in merged_df.columns:
        up_rate = merged_df[merged_df['logFC'] > 0]['any_anchor_overlap'].mean() * 100
        down_rate = merged_df[merged_df['logFC'] < 0]['any_anchor_overlap'].mean() * 100
        
        ax.bar(['Up-regulated', 'Down-regulated'], [up_rate, down_rate],
               color=['#d62728', '#2ca02c'], alpha=0.7)
        ax.set_ylabel('Insulator Overlap Rate (%)')
        ax.set_title('Direction-Specific Enrichment')
    
    # Plot 4: Distance distribution
    ax = fig.add_subplot(gs[1, :])
    plot_data = merged_df[~np.isinf(merged_df['min_dist_any'])].copy()
    
    if len(plot_data) > 0:
        bins = np.logspace(2, 6, 30)
        
        up_dists = plot_data[plot_data['logFC'] > 0]['min_dist_any']
        down_dists = plot_data[plot_data['logFC'] < 0]['min_dist_any']
        
        ax.hist([up_dists, down_dists], bins=bins, 
               label=['Up-regulated', 'Down-regulated'],
               color=['#d62728', '#2ca02c'], alpha=0.6)
        ax.set_xscale('log')
        ax.set_xlabel('Distance to Nearest Insulator (bp)')
        ax.set_ylabel('Number of Interactions')
        ax.set_title('Distance Distribution to Insulators')
        ax.legend()
        ax.axvline(10000, color='gray', linestyle='--', alpha=0.5, label='Window size')
    
    # Plot 5: LogFC vs insulators
    ax = fig.add_subplot(gs[2, 0])
    scatter = ax.scatter(merged_df['logFC'], merged_df['total_insulators'],
                        c=merged_df['any_anchor_overlap'], cmap='RdYlGn',
                        alpha=0.6)
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('Number of Insulators')
    ax.set_title('Effect Size vs Insulator Density')
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    plt.colorbar(scatter, ax=ax, label='Overlaps Insulator')
    
    # Plot 6: Overlap pattern
    ax = fig.add_subplot(gs[2, 1])
    overlap_types = ['No overlap', 'One anchor', 'Both anchors']
    counts = [
        (~merged_df['any_anchor_overlap']).sum(),
        (merged_df['any_anchor_overlap'] & ~merged_df['both_anchors_overlap']).sum(),
        merged_df['both_anchors_overlap'].sum()
    ]
    ax.bar(overlap_types, counts, color=['gray', 'orange', 'red'], alpha=0.7)
    ax.set_ylabel('Number of Interactions')
    ax.set_title('Insulator Overlap Pattern')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Plot 7: Cumulative distribution
    ax = fig.add_subplot(gs[2, 2])
    if len(plot_data) > 0:
        sorted_dists = np.sort(plot_data['min_dist_any'])
        cumsum = np.arange(1, len(sorted_dists) + 1) / len(sorted_dists)
        
        ax.plot(sorted_dists, cumsum * 100, linewidth=2)
        ax.set_xscale('log')
        ax.set_xlabel('Distance to Nearest Insulator (bp)')
        ax.set_ylabel('Cumulative Percentage')
        ax.set_title('Cumulative Distance Distribution')
        ax.grid(True, alpha=0.3)
        ax.axvline(10000, color='red', linestyle='--', alpha=0.5, label='10kb')
        ax.legend()
    
    plt.savefig(f"{output_prefix}/insulator_analysis.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nVisualization saved to {output_prefix}/insulator_analysis.pdf")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze insulator enrichment at differential contacts with permutation testing'
    )
    parser.add_argument('--interactions', required=True,
                       help='Significant differential interactions CSV file')
    parser.add_argument('--null_model', required=True,
                       help='Null model interactions CSV file')
    parser.add_argument('--class_i', required=True,
                       help='Class I insulator sites BED file (CTCF-dependent)')
    parser.add_argument('--class_ii', required=True,
                       help='Class II insulator sites BED file (CTCF-independent)')
    parser.add_argument('--genome', required=True,
                       help='Genome file for permutation testing')
    parser.add_argument('--window_size', type=int, default=10000,
                       help='Window size around insulators (bp)')
    parser.add_argument('--n_permutations', type=int, default=1000,
                       help='Number of permutations for significance testing')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significant interactions')
    parser.add_argument('--check_interaction', nargs=6, metavar=('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'),
                       help='Check specific interaction (e.g., 2L 10000000 10008000 X 22432000 22440000)')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("INSULATOR ENRICHMENT ANALYSIS")
    print("With Permutation Testing")
    print("OPTIMIZED VERSION")
    print("="*60)
    
    # Load data
    insulators = load_insulator_sites(args.class_i, args.class_ii)
    sig_interactions, all_interactions = load_differential_interactions(
        args.interactions, args.fdr_threshold
    )
    null_model = load_null_model(args.null_model)
    
    if len(sig_interactions) == 0:
        print("Error: No significant interactions to analyze")
        return 1
    
    # Analyze by insulator class
    class_results = analyze_by_insulator_class(
        sig_interactions, insulators, args.genome, 
        args.window_size, args.n_permutations
    )
    
    # Use 'All' insulators for general analysis
    if 'All' in class_results and class_results['All']:
        overlap_results = class_results['All']['overlap']
        perm_results = class_results['All']['permutation']
        
        # Direction analysis
        direction_results, merged_df = analyze_by_logfc_direction(
            sig_interactions, overlap_results
        )
        
        # Check specific interaction if requested
        if args.check_interaction:
            chr1, start1, end1, chr2, start2, end2 = args.check_interaction
            check_specific_interaction(
                sig_interactions, overlap_results,
                chr1, int(start1), int(end1),
                chr2, int(start2), int(end2)
            )
        
        # Create visualizations
        create_visualization(class_results, merged_df, perm_results, args.output_prefix)
        
        # Save results
        merged_df.to_csv(f"{args.output_prefix}/insulator_interactions.csv", index=False)
        
        # Save summary
        summary = {
            'total_interactions': len(all_interactions),
            'significant_interactions': len(sig_interactions),
            'window_size': args.window_size,
            'n_permutations': args.n_permutations
        }
        
        # Add class-specific results
        for ins_class, results in class_results.items():
            if results and results['permutation']:
                perm = results['permutation']
                summary[f'{ins_class}_n_sites'] = len(insulators[ins_class])
                summary[f'{ins_class}_observed_rate'] = perm['observed_rate']
                summary[f'{ins_class}_expected_rate'] = perm['expected_rate']
                summary[f'{ins_class}_enrichment'] = perm['enrichment']
                summary[f'{ins_class}_z_score'] = perm['z_score']
                summary[f'{ins_class}_p_value'] = perm['p_value']
        
        summary_df = pd.DataFrame([summary])
        summary_df.to_csv(f"{args.output_prefix}/insulator_summary.csv", index=False)
        
        # Create text summary
        with open(f"{args.output_prefix}/insulator_summary.txt", 'w') as f:
            f.write("INSULATOR ENRICHMENT ANALYSIS SUMMARY\n")
            f.write("="*60 + "\n\n")
            
            f.write("Analysis Parameters:\n")
            f.write(f"  Window size: {args.window_size} bp\n")
            f.write(f"  Permutations: {args.n_permutations}\n")
            f.write(f"  FDR threshold: {args.fdr_threshold}\n\n")
            
            f.write("Data Summary:\n")
            f.write(f"  Total interactions: {summary['total_interactions']}\n")
            f.write(f"  Significant interactions: {summary['significant_interactions']}\n\n")
            
            for ins_class in ['Class_I', 'Class_II', 'All']:
                if f'{ins_class}_p_value' in summary:
                    f.write(f"\n{ins_class} Insulators:\n")
                    f.write(f"  Number of sites: {summary[f'{ins_class}_n_sites']}\n")
                    f.write(f"  Observed overlap: {summary[f'{ins_class}_observed_rate']*100:.1f}%\n")
                    f.write(f"  Expected overlap: {summary[f'{ins_class}_expected_rate']*100:.1f}%\n")
                    f.write(f"  Enrichment: {summary[f'{ins_class}_enrichment']:.2f}x\n")
                    f.write(f"  Z-score: {summary[f'{ins_class}_z_score']:.2f}\n")
                    f.write(f"  P-value: {summary[f'{ins_class}_p_value']:.4f}\n")
                    f.write(f"  Significant: {'YES' if summary[f'{ins_class}_p_value'] < 0.05 else 'NO'}\n")
        
        print(f"\nAnalysis complete! Results saved to {args.output_prefix}/*")
        
    return 0

if __name__ == '__main__':
    exit(main())