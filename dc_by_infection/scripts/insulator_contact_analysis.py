#!/usr/bin/env python3
"""
Analyze insulator protein enrichment at differential chromatin contacts.
Uses permutation testing to assess statistical significance.

FIXED VERSION - Correctly counts overlaps per anchor instead of total insulators

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
    FIXED: Calculate overlap between interaction anchors and insulator sites.
    Now correctly counts overlaps per anchor instead of total insulators.
    """
    print(f"\nAnalyzing insulator overlap (window size: {window_size}bp)")
    
    if len(insulator_sites) == 0:
        print("No insulator sites available")
        return None
    
    # Create extended insulator regions
    insulators_extended = insulator_sites.copy()
    insulators_extended['start'] = (insulators_extended['start'] - window_size).clip(lower=0)
    insulators_extended['end'] = insulators_extended['end'] + window_size
    
    # Sort before creating BedTool
    insulators_extended = insulators_extended.sort_values(['chrom', 'start', 'end'])
    
    insulator_bt = pybedtools.BedTool.from_dataframe(
        insulators_extended[['chrom', 'start', 'end']]
    )
    
    # Create anchor dataframes with indices
    anchor1_df = interactions_df[['chr1', 'start1', 'end1']].copy()
    anchor1_df.columns = ['chrom', 'start', 'end']
    anchor1_df['idx'] = interactions_df.index
    
    anchor2_df = interactions_df[['chr2', 'start2', 'end2']].copy()
    anchor2_df.columns = ['chrom', 'start', 'end']
    anchor2_df['idx'] = interactions_df.index
    
    # Sort anchor dataframes
    anchor1_df = anchor1_df.sort_values(['chrom', 'start', 'end'])
    anchor2_df = anchor2_df.sort_values(['chrom', 'start', 'end'])
    
    anchor1_bt = pybedtools.BedTool.from_dataframe(anchor1_df)
    anchor2_bt = pybedtools.BedTool.from_dataframe(anchor2_df)
    
    # Use bedtools intersect with -c flag to count overlaps
    print("  Counting overlaps for anchor 1...")
    anchor1_overlaps = anchor1_bt.intersect(insulator_bt, c=True)
    
    print("  Counting overlaps for anchor 2...")
    anchor2_overlaps = anchor2_bt.intersect(insulator_bt, c=True)
    
    # FIXED: Parse bedtools output correctly
    # bedtools intersect -c adds the count as the last column
    anchor1_counts = {}
    for feature in anchor1_overlaps:
        # Convert to string and split to get all fields
        fields = str(feature).strip().split('\t')
        idx = int(fields[3])  # Index is in 4th column (0-indexed as 3)
        count = int(fields[4])  # Count is in 5th column (0-indexed as 4)
        anchor1_counts[idx] = count
    
    anchor2_counts = {}
    for feature in anchor2_overlaps:
        fields = str(feature).strip().split('\t')
        idx = int(fields[3])
        count = int(fields[4])
        anchor2_counts[idx] = count
    
    # DEBUG: Print some examples
    print(f"\n  DEBUG - Sample counts:")
    sample_indices = list(interactions_df.index[:3])
    for idx in sample_indices:
        a1_count = anchor1_counts.get(idx, 0)
        a2_count = anchor2_counts.get(idx, 0)
        print(f"    Interaction {idx}: anchor1={a1_count}, anchor2={a2_count}")
    
    # Calculate distances using bedtools closest
    print("  Calculating distances to nearest insulators...")
    insulator_sites_sorted = insulator_sites.sort_values(['chrom', 'start', 'end'])
    insulator_base_bt = pybedtools.BedTool.from_dataframe(
        insulator_sites_sorted[['chrom', 'start', 'end']]
    )
    
    anchor1_closest = anchor1_bt.closest(insulator_base_bt, d=True)
    anchor2_closest = anchor2_bt.closest(insulator_base_bt, d=True)
    
    # Extract distances
    anchor1_dists = {}
    for feature in anchor1_closest:
        fields = str(feature).strip().split('\t')
        idx = int(fields[3])
        dist = int(fields[-1])  # Distance is the last field
        if idx not in anchor1_dists or dist < anchor1_dists[idx]:
            anchor1_dists[idx] = dist
    
    anchor2_dists = {}
    for feature in anchor2_closest:
        fields = str(feature).strip().split('\t')
        idx = int(fields[3])
        dist = int(fields[-1])
        if idx not in anchor2_dists or dist < anchor2_dists[idx]:
            anchor2_dists[idx] = dist
    
    # Compile results
    results = []
    for idx in interactions_df.index:
        n_ins_1 = anchor1_counts.get(idx, 0)
        n_ins_2 = anchor2_counts.get(idx, 0)
        
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
    print(f"  Mean insulators per anchor1: {results_df['n_ins_anchor1'].mean():.2f}")
    print(f"  Mean insulators per anchor2: {results_df['n_ins_anchor2'].mean():.2f}")
    
    return results_df

def calculate_genomic_background_model(interactions_df, insulator_sites, genome_file, 
                                       overlap_df, window_size=10000):
    """
    Calculate expected overlap probability based on genome-wide insulator distribution.
    
    This answers: "Given the distribution of insulators across the genome, 
    what is the likelihood of an interaction overlapping insulators?"
    
    Uses binomial probability based on:
    - Total genome size
    - Total insulator coverage (with windows)
    - Interaction anchor sizes
    """
    print(f"\n{'='*60}")
    print("GENOMIC BACKGROUND MODEL")
    print("Calculating expected overlap based on insulator density")
    print(f"{'='*60}")
    
    # Load genome sizes
    genome_sizes = {}
    total_genome_size = 0
    
    with open(genome_file, 'r') as f:
        for line in f:
            chrom, size = line.strip().split()
            genome_sizes[chrom] = int(size)
            total_genome_size += int(size)
    
    print(f"\nGenome Information:")
    print(f"  Total genome size: {total_genome_size:,} bp")
    print(f"  Number of chromosomes: {len(genome_sizes)}")
    
    # Calculate total insulator coverage (with windows)
    insulators_extended = insulator_sites.copy()
    insulators_extended['start'] = (insulators_extended['start'] - window_size).clip(lower=0)
    insulators_extended['end'] = insulators_extended['end'] + window_size
    
    # Use bedtools to merge overlapping insulator regions
    insulators_extended = insulators_extended.sort_values(['chrom', 'start', 'end'])
    insulator_bt = pybedtools.BedTool.from_dataframe(
        insulators_extended[['chrom', 'start', 'end']]
    ).merge()
    
    # Calculate total insulator coverage
    total_insulator_coverage = 0
    insulator_coverage_by_chr = {}
    
    for feature in insulator_bt:
        chrom = feature.chrom
        coverage = feature.end - feature.start
        total_insulator_coverage += coverage
        
        if chrom not in insulator_coverage_by_chr:
            insulator_coverage_by_chr[chrom] = 0
        insulator_coverage_by_chr[chrom] += coverage
    
    # Calculate genome-wide insulator density
    genome_insulator_density = total_insulator_coverage / total_genome_size
    
    print(f"\nInsulator Coverage:")
    print(f"  Total insulator coverage: {total_insulator_coverage:,} bp")
    print(f"  Genome-wide density: {genome_insulator_density*100:.2f}%")
    print(f"\nPer-chromosome density:")
    for chrom in sorted(insulator_coverage_by_chr.keys()):
        if chrom in genome_sizes:
            chr_density = insulator_coverage_by_chr[chrom] / genome_sizes[chrom]
            print(f"    {chrom}: {chr_density*100:.2f}%")
    
    # Calculate expected probabilities for each interaction
    results = []
    
    for idx, row in interactions_df.iterrows():
        # Get chromosome-specific densities
        chr1 = row['chr1']
        chr2 = row['chr2']
        
        # Use chromosome-specific density if available, otherwise genome-wide
        density1 = insulator_coverage_by_chr.get(chr1, 0) / genome_sizes.get(chr1, 1) if chr1 in genome_sizes else genome_insulator_density
        density2 = insulator_coverage_by_chr.get(chr2, 0) / genome_sizes.get(chr2, 1) if chr2 in genome_sizes else genome_insulator_density
        
        # Calculate anchor sizes
        anchor1_size = row['end1'] - row['start1']
        anchor2_size = row['end2'] - row['start2']
        
        # Expected probability that each anchor overlaps an insulator
        # P(overlap) ≈ insulator_density (for small regions)
        p_anchor1_overlap = density1
        p_anchor2_overlap = density2
        
        # Expected probability that ANY anchor overlaps
        p_any_overlap = p_anchor1_overlap + p_anchor2_overlap - (p_anchor1_overlap * p_anchor2_overlap)
        
        # Expected probability that BOTH anchors overlap
        p_both_overlap = p_anchor1_overlap * p_anchor2_overlap
        
        # Get observed overlap from overlap_df
        overlap_info = overlap_df[overlap_df['interaction_idx'] == idx]
        if len(overlap_info) > 0:
            overlap_info = overlap_info.iloc[0]
            observed_any = overlap_info['any_anchor_overlap']
            observed_both = overlap_info['both_anchors_overlap']
        else:
            observed_any = False
            observed_both = False
        
        results.append({
            'interaction_idx': idx,
            'chr1': chr1,
            'chr2': chr2,
            'density1': density1,
            'density2': density2,
            'p_anchor1_overlap': p_anchor1_overlap,
            'p_anchor2_overlap': p_anchor2_overlap,
            'p_any_overlap': p_any_overlap,
            'p_both_overlap': p_both_overlap,
            'observed_any': observed_any,
            'observed_both': observed_both,
            'logFC': row['logFC']
        })
    
    results_df = pd.DataFrame(results)
    
    # Calculate overall statistics
    expected_any = results_df['p_any_overlap'].sum()
    observed_any = results_df['observed_any'].sum()
    
    expected_both = results_df['p_both_overlap'].sum()
    observed_both = results_df['observed_both'].sum()
    
    # Binomial test for significance
    from scipy.stats import binom_test
    
    n_interactions = len(results_df)
    p_any_mean = results_df['p_any_overlap'].mean()
    p_both_mean = results_df['p_both_overlap'].mean()
    
    p_value_any = binom_test(observed_any, n_interactions, p_any_mean, alternative='greater')
    p_value_both = binom_test(observed_both, n_interactions, p_both_mean, alternative='greater')
    
    # Calculate enrichment
    enrichment_any = observed_any / expected_any if expected_any > 0 else np.inf
    enrichment_both = observed_both / expected_both if expected_both > 0 else np.inf
    
    print(f"\n{'='*60}")
    print("BACKGROUND MODEL RESULTS")
    print(f"{'='*60}")
    print(f"\nAny Anchor Overlap:")
    print(f"  Expected: {expected_any:.1f} interactions ({p_any_mean*100:.1f}%)")
    print(f"  Observed: {observed_any} interactions ({observed_any/n_interactions*100:.1f}%)")
    print(f"  Enrichment: {enrichment_any:.2f}x")
    print(f"  Binomial p-value: {p_value_any:.4e}")
    print(f"  Significant: {'YES' if p_value_any < 0.05 else 'NO'}")
    
    print(f"\nBoth Anchors Overlap:")
    print(f"  Expected: {expected_both:.1f} interactions ({p_both_mean*100:.2f}%)")
    print(f"  Observed: {observed_both} interactions ({observed_both/n_interactions*100:.1f}%)")
    print(f"  Enrichment: {enrichment_both:.2f}x")
    print(f"  Binomial p-value: {p_value_both:.4e}")
    print(f"  Significant: {'YES' if p_value_both < 0.05 else 'NO'}")
    
    summary = {
        'genome_size': total_genome_size,
        'insulator_coverage': total_insulator_coverage,
        'insulator_density': genome_insulator_density,
        'n_interactions': n_interactions,
        'expected_any': expected_any,
        'observed_any': observed_any,
        'enrichment_any': enrichment_any,
        'p_value_any': p_value_any,
        'expected_both': expected_both,
        'observed_both': observed_both,
        'enrichment_both': enrichment_both,
        'p_value_both': p_value_both
    }
    
    return results_df, summary

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
    
    # Sort before creating BedTool
    interactions_bed_df = interactions_bed_df.sort_values(['chrom', 'start', 'end'])
    
    interactions_bt = pybedtools.BedTool.from_dataframe(interactions_bed_df)
    
    # Extend insulator sites
    insulators_extended = insulator_sites.copy()
    insulators_extended['start'] = (insulators_extended['start'] - window_size).clip(lower=0)
    insulators_extended['end'] = insulators_extended['end'] + window_size
    
    # Sort before creating BedTool
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
    merged['both_anchors_overlap'] = overlap_df['both_anchors_overlap'].values
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

def check_specific_interaction(interactions_df, overlap_df, background_results_df, 
                               chr1, start1, end1, chr2, start2, end2):
    """
    Check if a specific interaction (like the NDF example) has insulator enrichment.
    Example: 2L:10000000-10008000_X:22432000-22440000
    
    Now includes p-value for whether insulator overlap is significant for THIS specific interaction.
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
        
        # Get background model results for this interaction
        background_info = background_results_df[background_results_df['interaction_idx'] == idx]
        
        print(f"\nInteraction Details:")
        print(f"  Coordinates: {interaction['chr1']}:{interaction['start1']}-{interaction['end1']}_"
              f"{interaction['chr2']}:{interaction['start2']}-{interaction['end2']}")
        print(f"  logFC: {interaction['logFC']:.2f}")
        print(f"  FDR (differential): {interaction['FDR']:.2e}")
        
        print(f"\nInsulator Enrichment:")
        print(f"  Anchor 1 overlaps insulator: {'YES' if overlap_info['anchor1_overlap'] else 'NO'}")
        print(f"  Anchor 2 overlaps insulator: {'YES' if overlap_info['anchor2_overlap'] else 'NO'}")
        print(f"  Number of insulators (anchor 1): {overlap_info['n_ins_anchor1']}")
        print(f"  Number of insulators (anchor 2): {overlap_info['n_ins_anchor2']}")
        print(f"  Distance to nearest (anchor 1): {overlap_info['min_dist_anchor1']:.0f} bp")
        print(f"  Distance to nearest (anchor 2): {overlap_info['min_dist_anchor2']:.0f} bp")
        
        # Add background model statistics for this specific interaction
        if len(background_info) > 0:
            background_info = background_info.iloc[0]
            
            print(f"\nBackground Model for This Interaction:")
            print(f"  Insulator density on {interaction['chr1']}: {background_info['density1']*100:.2f}%")
            print(f"  Insulator density on {interaction['chr2']}: {background_info['density2']*100:.2f}%")
            print(f"  Expected P(any anchor overlaps): {background_info['p_any_overlap']*100:.2f}%")
            print(f"  Expected P(both anchors overlap): {background_info['p_both_overlap']*100:.2f}%")
            
            # Calculate p-value for this specific interaction
            # Using binomial test: is observing overlap at this one interaction surprising?
            from scipy.stats import binom_test
            
            observed_any = int(overlap_info['any_anchor_overlap'])
            expected_prob_any = background_info['p_any_overlap']
            
            # Binomial test for single trial
            if observed_any:
                p_value_any = expected_prob_any if expected_prob_any > 0 else 1.0
                # For a single Bernoulli trial, p-value = probability of success
                enrichment_interpretation = f"{1/expected_prob_any:.2f}x less likely by chance" if expected_prob_any > 0 else "N/A"
            else:
                p_value_any = 1 - expected_prob_any
                enrichment_interpretation = "No overlap (expected)"
            
            print(f"\nStatistical Significance:")
            print(f"  Observed overlap: {'YES' if observed_any else 'NO'}")
            print(f"  Likelihood by chance: {expected_prob_any*100:.2f}%")
            print(f"  P-value: {min(expected_prob_any, 1-expected_prob_any):.4f}")
            
            if expected_prob_any < 0.05 and observed_any:
                print(f"  Interpretation: Overlap is RARE by chance (expected in <5% of interactions)")
                print(f"                  This overlap is likely BIOLOGICALLY MEANINGFUL")
            elif expected_prob_any < 0.20 and observed_any:
                print(f"  Interpretation: Overlap is UNCOMMON by chance (expected in {expected_prob_any*100:.1f}% of interactions)")
                print(f"                  This overlap may be biologically relevant")
            elif observed_any:
                print(f"  Interpretation: Overlap is COMMON by chance (expected in {expected_prob_any*100:.1f}% of interactions)")
                print(f"                  This overlap may be coincidental")
            else:
                print(f"  Interpretation: No overlap observed (expected given low insulator density)")
        
        return {
            'interaction': interaction,
            'overlap': overlap_info,
            'background': background_info if len(background_info) > 0 else None
        }
    
    return None

def create_visualization(class_results, merged_df, perm_results, output_prefix):
    """Create individual plot files for each visualization"""
    
    print("\nGenerating individual plot files...")
    
    # Plot 1: Enrichment by insulator class
    fig, ax = plt.subplots(figsize=(8, 6))
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
    ax.set_ylabel('Enrichment vs Null', fontsize=12)
    ax.set_title('Insulator Enrichment by Class', fontsize=14)
    ax.legend()
    
    # Add significance stars
    for bar, pval in zip(bars, p_values):
        stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                stars, ha='center', va='bottom', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/1_enrichment_by_class.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: 1_enrichment_by_class.pdf")
    
    # Plot 2: Null distribution with observed value
    if perm_results and perm_results['null_distribution'] is not None:
        fig, ax = plt.subplots(figsize=(8, 6))
        null_dist = perm_results['null_distribution'] * 100
        ax.hist(null_dist, bins=30, alpha=0.7, color='gray', edgecolor='black')
        ax.axvline(perm_results['observed_rate']*100, color='red', 
                  linestyle='--', linewidth=2, label='Observed')
        ax.axvline(perm_results['expected_rate']*100, color='blue',
                  linestyle='--', linewidth=2, label='Expected (median)')
        ax.set_xlabel('Overlap Rate (%)', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title('Permutation Test Null Distribution', fontsize=14)
        ax.legend()
        plt.tight_layout()
        plt.savefig(f"{output_prefix}/2_null_distribution.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: 2_null_distribution.pdf")
    
    # Plot 3: Overlap rate by direction
    if 'logFC' in merged_df.columns:
        fig, ax = plt.subplots(figsize=(6, 6))
        up_rate = merged_df[merged_df['logFC'] > 0]['any_anchor_overlap'].mean() * 100
        down_rate = merged_df[merged_df['logFC'] < 0]['any_anchor_overlap'].mean() * 100

        ax.bar(['Uninf.', 'wMel'], [up_rate, down_rate],
               color=['#8fcb84', '#09aa4b'], alpha=0.7)
        ax.set_ylabel('Insulator Overlap Rate (%)', fontsize=12)
        ax.set_title('Direction-Specific Enrichment', fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}/3_direction_specific.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: 3_direction_specific.pdf")
    
    # Plot 4: Distance distribution
    plot_data = merged_df[~np.isinf(merged_df['min_dist_any'])].copy()
    
    if len(plot_data) > 0:
        fig, ax = plt.subplots(figsize=(10, 6))
        bins = np.logspace(2, 6, 30)
        
        up_dists = plot_data[plot_data['logFC'] > 0]['min_dist_any']
        down_dists = plot_data[plot_data['logFC'] < 0]['min_dist_any']
        
        ax.hist([up_dists, down_dists], bins=bins, 
               label=['Uninf.', 'wMel'],
               color=['#8fcb84', '#09aa4b'], alpha=0.6)
        ax.set_xscale('log')
        ax.set_xlabel('Distance to Nearest Insulator (bp)', fontsize=12)
        ax.set_ylabel('Number of Interactions', fontsize=12)
        ax.set_title('Distance Distribution to Insulators', fontsize=14)
        ax.legend()
        ax.axvline(1000, color='gray', linestyle='--', alpha=0.5, label='1kb')
        plt.tight_layout()
        plt.savefig(f"{output_prefix}/4_distance_distribution.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: 4_distance_distribution.pdf")
    
    # Plot 5: LogFC vs insulators
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(merged_df['logFC'], merged_df['total_insulators'],
                        c=merged_df['any_anchor_overlap'], cmap='RdYlGn',
                        alpha=0.6, s=50)
    ax.set_xlabel('log2 Fold Change', fontsize=12)
    ax.set_ylabel('Number of Insulators', fontsize=12)
    ax.set_title('Effect Size vs Insulator Density', fontsize=14)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Overlaps Insulator', fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/5_logfc_vs_insulators.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: 5_logfc_vs_insulators.pdf")
    
    # Plot 6: Overlap pattern
    fig, ax = plt.subplots(figsize=(8, 6))
    overlap_types = ['No overlap', 'One anchor', 'Both anchors']
    counts = [
        (~merged_df['any_anchor_overlap']).sum(),
        (merged_df['any_anchor_overlap'] & ~merged_df['both_anchors_overlap']).sum(),
        merged_df['both_anchors_overlap'].sum()
    ]
    ax.bar(overlap_types, counts, color=['gray', 'orange', 'red'], alpha=0.7)
    ax.set_ylabel('Number of Interactions', fontsize=12)
    ax.set_title('Insulator Overlap Pattern', fontsize=14)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/6_overlap_pattern.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: 6_overlap_pattern.pdf")
    
    # Plot 7: Cumulative distribution
    if len(plot_data) > 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        sorted_dists = np.sort(plot_data['min_dist_any'])
        cumsum = np.arange(1, len(sorted_dists) + 1) / len(sorted_dists)
        
        ax.plot(sorted_dists, cumsum * 100, linewidth=2)
        ax.set_xscale('log')
        ax.set_xlabel('Distance to Nearest Insulator (bp)', fontsize=12)
        ax.set_ylabel('Cumulative Percentage', fontsize=12)
        ax.set_title('Cumulative Distance Distribution', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.axvline(10000, color='red', linestyle='--', alpha=0.5, label='10kb')
        ax.axvline(1000, color='orange', linestyle='--', alpha=0.5, label='1kb')
        ax.legend()
        plt.tight_layout()
        plt.savefig(f"{output_prefix}/7_cumulative_distance.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: 7_cumulative_distance.pdf")
    
    print(f"\nAll plots saved to {output_prefix}/")
    print("Plot files:")
    print("  1_enrichment_by_class.pdf")
    print("  2_null_distribution.pdf")
    print("  3_direction_specific.pdf")
    print("  4_distance_distribution.pdf")
    print("  5_logfc_vs_insulators.pdf")
    print("  6_overlap_pattern.pdf")
    print("  7_cumulative_distance.pdf")

def create_background_model_plots(background_results_df, background_summary, output_prefix):
    """Create plots for genomic background model analysis"""
    
    print("\nGenerating background model plots...")
    
    # Plot 1: Expected vs Observed overlap rates
    fig, ax = plt.subplots(figsize=(8, 6))
    
    categories = ['Any Anchor\nOverlap', 'Both Anchors\nOverlap']
    expected = [
        background_summary['expected_any'] / background_summary['n_interactions'] * 100,
        background_summary['expected_both'] / background_summary['n_interactions'] * 100
    ]
    observed = [
        background_summary['observed_any'] / background_summary['n_interactions'] * 100,
        background_summary['observed_both'] / background_summary['n_interactions'] * 100
    ]
    
    x = np.arange(len(categories))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, expected, width, label='Expected (background)', 
                   color='lightblue', alpha=0.7)
    bars2 = ax.bar(x + width/2, observed, width, label='Observed', 
                   color='darkblue', alpha=0.7)
    
    ax.set_ylabel('Percentage of Interactions', fontsize=12)
    ax.set_title('Observed vs Expected Insulator Overlap\n(Based on Genome-wide Insulator Density)', 
                 fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend()
    
    # Add p-value annotations
    p_values = [background_summary['p_value_any'], background_summary['p_value_both']]
    for i, (cat, pval) in enumerate(zip(categories, p_values)):
        stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        y_pos = max(expected[i], observed[i]) + 2
        ax.text(i, y_pos, stars, ha='center', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/8_background_model_comparison.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: 8_background_model_comparison.pdf")
    
    # Plot 2: Enrichment factors
    fig, ax = plt.subplots(figsize=(8, 6))
    
    enrichments = [
        background_summary['enrichment_any'],
        background_summary['enrichment_both']
    ]
    
    colors_sig = ['darkgreen' if p < 0.05 else 'gray' for p in p_values]
    bars = ax.bar(categories, enrichments, color=colors_sig, alpha=0.7)
    ax.axhline(y=1, color='red', linestyle='--', linewidth=2, label='No enrichment')
    ax.set_ylabel('Enrichment Factor', fontsize=12)
    ax.set_title('Enrichment Over Genomic Background', fontsize=14)
    ax.legend()
    
    # Add enrichment values on bars
    for bar, enrich, pval in zip(bars, enrichments, p_values):
        height = bar.get_height()
        stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{enrich:.2f}x\n{stars}', ha='center', va='bottom', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/9_enrichment_over_background.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: 9_enrichment_over_background.pdf")
    
    # Plot 3: Distribution of expected probabilities
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(background_results_df['p_any_overlap'] * 100, bins=30, 
           alpha=0.6, label='Any anchor', color='blue', edgecolor='black')
    ax.hist(background_results_df['p_both_overlap'] * 100, bins=30, 
           alpha=0.6, label='Both anchors', color='green', edgecolor='black')
    
    ax.set_xlabel('Expected Overlap Probability (%)', fontsize=12)
    ax.set_ylabel('Number of Interactions', fontsize=12)
    ax.set_title('Distribution of Expected Overlap Probabilities', fontsize=14)
    ax.legend()
    
    # Add mean lines
    ax.axvline(background_results_df['p_any_overlap'].mean() * 100, 
              color='blue', linestyle='--', linewidth=2, alpha=0.7)
    ax.axvline(background_results_df['p_both_overlap'].mean() * 100, 
              color='green', linestyle='--', linewidth=2, alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/10_expected_probability_distribution.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: 10_expected_probability_distribution.pdf")
    
    # Plot 4: Chromosome-specific densities
    fig, ax = plt.subplots(figsize=(10, 6))
    
    chr_densities = background_results_df.groupby('chr1').agg({
        'density1': 'mean'
    }).reset_index()
    chr_densities = chr_densities.sort_values('density1', ascending=False)
    
    ax.bar(chr_densities['chr1'], chr_densities['density1'] * 100, 
           color='steelblue', alpha=0.7)
    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('Insulator Density (%)', fontsize=12)
    ax.set_title('Chromosome-Specific Insulator Density', fontsize=14)
    ax.axhline(y=background_summary['insulator_density'] * 100, 
              color='red', linestyle='--', linewidth=2, label='Genome-wide average')
    ax.legend()
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/11_chromosome_densities.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: 11_chromosome_densities.pdf")
    
    # Plot 5: Observed vs Expected per interaction (scatter)
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Create binary observed values
    observed_any_binary = background_results_df['observed_any'].astype(int)
    
    # Jitter for visualization
    jitter = 0.1
    x_jittered = background_results_df['p_any_overlap'] * 100 + np.random.uniform(-jitter, jitter, len(background_results_df))
    y_jittered = observed_any_binary + np.random.uniform(-jitter, jitter, len(background_results_df))
    
    colors = ['green' if obs else 'gray' for obs in background_results_df['observed_any']]
    ax.scatter(x_jittered, y_jittered, alpha=0.5, c=colors, s=30)
    
    ax.set_xlabel('Expected Overlap Probability (%)', fontsize=12)
    ax.set_ylabel('Observed Overlap (0=No, 1=Yes)', fontsize=12)
    ax.set_title('Expected Probability vs Observed Overlap\n(Each point = one interaction)', 
                 fontsize=14)
    ax.set_ylim(-0.3, 1.3)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['No', 'Yes'])
    
    # Add reference line at mean expected probability
    mean_prob = background_results_df['p_any_overlap'].mean() * 100
    ax.axvline(mean_prob, color='blue', linestyle='--', linewidth=2, 
              alpha=0.7, label=f'Mean expected ({mean_prob:.1f}%)')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}/12_expected_vs_observed_scatter.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: 12_expected_vs_observed_scatter.pdf")
    
    print(f"\nBackground model plots saved to {output_prefix}/")
    print("Additional plot files:")
    print("  8_background_model_comparison.pdf")
    print("  9_enrichment_over_background.pdf")
    print("  10_expected_probability_distribution.pdf")
    print("  11_chromosome_densities.pdf")
    print("  12_expected_vs_observed_scatter.pdf")

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
    output_dir = Path(args.output_prefix)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("INSULATOR ENRICHMENT ANALYSIS")
    print("With Permutation Testing")
    print("FIXED VERSION - Correct overlap counting")
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
        
        # Calculate genomic background model
        print("\n" + "="*60)
        print("CALCULATING GENOMIC BACKGROUND MODEL")
        print("="*60)
        background_results_df, background_summary = calculate_genomic_background_model(
            sig_interactions, insulators['All'], args.genome, 
            overlap_results, args.window_size
        )
        
        # Direction analysis
        direction_results, merged_df = analyze_by_logfc_direction(
            sig_interactions, overlap_results
        )
        
        # Check specific interaction if requested
        if args.check_interaction:
            chr1, start1, end1, chr2, start2, end2 = args.check_interaction
            check_specific_interaction(
                sig_interactions, overlap_results, background_results_df,
                chr1, int(start1), int(end1),
                chr2, int(start2), int(end2)
            )
        
        # Create visualizations
        create_visualization(class_results, merged_df, perm_results, args.output_prefix)
        
        # Create background model visualizations
        create_background_model_plots(background_results_df, background_summary, args.output_prefix)
        
        # Save results
        merged_df.to_csv(f"{args.output_prefix}/insulator_interactions.csv", index=False)
        background_results_df.to_csv(f"{args.output_prefix}/background_model_results.csv", index=False)
        
        # Save summary
        summary = {
            'total_interactions': len(all_interactions),
            'significant_interactions': len(sig_interactions),
            'window_size': args.window_size,
            'n_permutations': args.n_permutations
        }
        
        # Add background model results to summary
        summary.update({
            'genome_size': background_summary['genome_size'],
            'insulator_coverage': background_summary['insulator_coverage'],
            'insulator_density_percent': background_summary['insulator_density'] * 100,
            'bg_expected_any': background_summary['expected_any'],
            'bg_observed_any': background_summary['observed_any'],
            'bg_enrichment_any': background_summary['enrichment_any'],
            'bg_pvalue_any': background_summary['p_value_any'],
            'bg_expected_both': background_summary['expected_both'],
            'bg_observed_both': background_summary['observed_both'],
            'bg_enrichment_both': background_summary['enrichment_both'],
            'bg_pvalue_both': background_summary['p_value_both']
        })
        
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
            
            f.write("="*60 + "\n")
            f.write("GENOMIC BACKGROUND MODEL\n")
            f.write("="*60 + "\n\n")
            f.write("This analysis answers: Given the genome-wide distribution of\n")
            f.write("insulators, what is the expected probability of overlap?\n\n")
            
            f.write(f"Genome Information:\n")
            f.write(f"  Total genome size: {summary['genome_size']:,} bp\n")
            f.write(f"  Insulator coverage: {summary['insulator_coverage']:,} bp\n")
            f.write(f"  Genome-wide insulator density: {summary['insulator_density_percent']:.2f}%\n\n")
            
            f.write(f"Any Anchor Overlap:\n")
            f.write(f"  Expected: {summary['bg_expected_any']:.1f} interactions\n")
            f.write(f"  Observed: {summary['bg_observed_any']} interactions\n")
            f.write(f"  Enrichment: {summary['bg_enrichment_any']:.2f}x\n")
            f.write(f"  P-value: {summary['bg_pvalue_any']:.4e}\n")
            f.write(f"  Significant: {'YES' if summary['bg_pvalue_any'] < 0.05 else 'NO'}\n\n")
            
            f.write(f"Both Anchors Overlap:\n")
            f.write(f"  Expected: {summary['bg_expected_both']:.1f} interactions\n")
            f.write(f"  Observed: {summary['bg_observed_both']} interactions\n")
            f.write(f"  Enrichment: {summary['bg_enrichment_both']:.2f}x\n")
            f.write(f"  P-value: {summary['bg_pvalue_both']:.4e}\n")
            f.write(f"  Significant: {'YES' if summary['bg_pvalue_both'] < 0.05 else 'NO'}\n\n")
            
            f.write("="*60 + "\n")
            f.write("PERMUTATION TEST RESULTS\n")
            f.write("="*60 + "\n\n")
            
            for ins_class in ['Class_I', 'Class_II', 'All']:
                if f'{ins_class}_p_value' in summary:
                    f.write(f"\n{ins_class} Insulators:\n")
                    f.write(f"  Number of sites: {summary[f'{ins_class}_n_sites']}\n")
                    f.write(f"  Observed overlap: {summary[f'{ins_class}_observed_rate']*100:.1f}%\n")
                    f.write(f"  Expected overlap (permutation): {summary[f'{ins_class}_expected_rate']*100:.1f}%\n")
                    f.write(f"  Enrichment: {summary[f'{ins_class}_enrichment']:.2f}x\n")
                    f.write(f"  Z-score: {summary[f'{ins_class}_z_score']:.2f}\n")
                    f.write(f"  P-value: {summary[f'{ins_class}_p_value']:.4f}\n")
                    f.write(f"  Significant: {'YES' if summary[f'{ins_class}_p_value'] < 0.05 else 'NO'}\n")
            
            f.write("\n" + "="*60 + "\n")
            f.write("INTERPRETATION\n")
            f.write("="*60 + "\n\n")
            f.write("Background Model: Tests if interactions overlap insulators more than\n")
            f.write("expected based on genome-wide insulator density.\n\n")
            f.write("Permutation Test: Tests if interactions overlap insulators more than\n")
            f.write("random genomic regions of the same size.\n\n")
            
            if summary['bg_pvalue_any'] < 0.05 and summary['All_p_value'] < 0.05:
                f.write("Both tests are significant, indicating strong evidence that\n")
                f.write("differential interactions are enriched at insulator sites.\n")
            elif summary['bg_pvalue_any'] < 0.05:
                f.write("Background model is significant: interactions overlap insulators\n")
                f.write("more than expected from genome-wide insulator density.\n")
            elif summary['All_p_value'] < 0.05:
                f.write("Permutation test is significant: interactions overlap insulators\n")
                f.write("more than random genomic regions.\n")
            else:
                f.write("Neither test is significant: insulator enrichment is not\n")
                f.write("statistically significant.\n")
        
        print(f"\nAnalysis complete! Results saved to {args.output_prefix}/*")
        
    return 0

if __name__ == '__main__':
    exit(main())
