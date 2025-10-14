#!/usr/bin/env python3
"""
This is a working version as of 6/16/25
Compare chromatin structure across conditions using existing diffHic results.
Uses null_model_results.csv for null model and all_results_combined.csv for interaction data.

MODIFIED: Now accepts variable number of conditions (not hardcoded to 4)
"""
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import cooler
import cooltools
from cooltools import expected_cis
from cooltools import dots
from cooltools import eigs_cis
from cooltools import insulation

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool
import bioframe

def get_closest_resolution(mcool_file, target_resolution):
    """Get the closest available resolution in the mcool file."""
    try:
        f = cooler.fileops.list_coolers(mcool_file)
        available_resolutions = [int(path.split('/')[-1]) for path in f]
        available_resolutions.sort()
        
        # Find the closest resolution
        closest = min(available_resolutions, key=lambda x: abs(x - target_resolution))
        print(f"Target resolution: {target_resolution}, using: {closest}")
        return closest
    except Exception as e:
        print(f"Error getting resolutions from {mcool_file}: {e}")
        return None

def load_diffhic_results(results_file):
    """Load the combined diffHic results."""
    print(f"Loading diffHic results from {results_file}...")
    results_df = pd.read_csv(results_file)
    
    # Filter for significant interactions
    significant_results = results_df[
        (results_df['FDR'] < 0.05) & 
        (abs(results_df['logFC']) > 1)
    ].copy()
    
    print(f"Loaded {len(results_df)} total interactions, {len(significant_results)} significant")
    return results_df, significant_results

def load_null_model(null_file):
    """Load null model results from diffHic analysis"""
    print(f"Loading null model from {null_file}...")
    null_df = pd.read_csv(null_file)
    print(f"Loaded null model with {len(null_df)} entries")
    return null_df

def identify_tad_boundaries_from_interactions(significant_interactions, conditions, window_size=50000):
    """
    Identify TAD boundaries based on significant differential interactions.
    TAD boundaries are regions with enriched differential interactions.
    """
    print("Identifying TAD boundaries from differential interactions...")
    
    all_boundaries = {}
    
    # Check if there's an infection/condition column
    condition_col = None
    for possible_col in ['infection', 'condition', 'comparison', 'contrast', 'group']:
        if possible_col in significant_interactions.columns:
            condition_col = possible_col
            print(f"  Found condition column: {condition_col}")
            break
    
    # If we have a condition column, group by it
    if condition_col is not None:
        groups_to_process = significant_interactions[condition_col].unique()
        print(f"  Processing {len(groups_to_process)} groups: {groups_to_process}")
        
        for group_name in groups_to_process:
            if pd.isna(group_name):
                continue
                
            print(f"  Processing {group_name}...")
            group_data = significant_interactions[
                significant_interactions[condition_col] == group_name
            ].copy()
            
            # Process this group
            boundaries = process_interactions_for_boundaries(group_data, window_size)
            
            if boundaries:
                boundaries_df = pd.DataFrame(boundaries)
                # Clean up the condition name
                clean_name = str(group_name).replace('infectionJW18', '').replace('infection', '')
                boundaries_df['condition'] = clean_name
                all_boundaries[clean_name] = boundaries_df
                print(f"  Found {len(boundaries_df)} boundaries for {clean_name}")
    else:
        # No condition column - process all interactions together for each specified condition
        # This is a fallback - we'll create boundaries based on all significant interactions
        print(f"  No condition column found. Processing all interactions together.")
        print(f"  Will assign boundaries to conditions based on logFC direction")
        
        # Separate by logFC direction to approximate different conditions
        for condition in conditions:
            print(f"  Creating boundaries for {condition}...")
            if condition in ['DOX', 'uninfected']:
                # Use all interactions for reference
                condition_data = significant_interactions.copy()
            else:
                # Use positive logFC interactions for infected conditions
                condition_data = significant_interactions[
                    significant_interactions['logFC'] > 0
                ].copy()
            
            boundaries = process_interactions_for_boundaries(condition_data, window_size)
            
            if boundaries:
                boundaries_df = pd.DataFrame(boundaries)
                boundaries_df['condition'] = condition
                all_boundaries[condition] = boundaries_df
                print(f"  Found {len(boundaries_df)} boundaries for {condition}")
    
    return all_boundaries

def process_interactions_for_boundaries(interaction_data, window_size=50000):
    """
    Helper function to process interactions and identify boundaries.
    """
    boundaries = []
    
    # Process each chromosome
    for chrom in interaction_data['chr1'].unique():
        chrom_interactions = interaction_data[
            (interaction_data['chr1'] == chrom) & 
            (interaction_data['chr2'] == chrom)  # cis interactions only
        ].copy()
        
        if len(chrom_interactions) == 0:
            continue
        
        # Create bins across the chromosome
        max_pos = max(
            chrom_interactions['end1'].max(),
            chrom_interactions['end2'].max()
        )
        bins = np.arange(0, max_pos + window_size, window_size)
        
        # Count interactions in each bin
        bin_counts = np.zeros(len(bins) - 1)
        
        for _, interaction in chrom_interactions.iterrows():
            # Find which bins this interaction overlaps
            start_bin = np.searchsorted(bins, interaction['start1'], side='right') - 1
            end_bin = np.searchsorted(bins, interaction['end2'], side='right') - 1
            
            # Add to bins
            for b in range(max(0, start_bin), min(len(bin_counts), end_bin + 1)):
                bin_counts[b] += abs(interaction['logFC'])
        
        # Find peaks in interaction density (potential boundaries)
        from scipy.signal import find_peaks
        peaks, properties = find_peaks(
            bin_counts, 
            height=np.percentile(bin_counts[bin_counts > 0], 75) if np.any(bin_counts > 0) else 0,
            distance=2  # Minimum 2 bins apart
        )
        
        # Convert peaks to genomic coordinates
        for peak in peaks:
            if peak < len(bins) - 1:
                boundary_start = bins[peak]
                boundary_end = bins[peak + 1]
                
                boundaries.append({
                    'chrom': chrom,
                    'start': boundary_start,
                    'end': boundary_end,
                    'boundary_strength': bin_counts[peak],
                    'n_interactions': np.sum(
                        (chrom_interactions['start1'] >= boundary_start) &
                        (chrom_interactions['end2'] <= boundary_end)
                    )
                })
    
    return boundaries

def compare_tad_boundaries(tad_boundaries, null_model, fdr_threshold=0.05):
    """
    Compare TAD boundaries between conditions using the null model.
    """
    print("Comparing TAD boundaries between conditions...")
    
    if not tad_boundaries:
        print("No TAD boundaries to compare")
        return pd.DataFrame()
    
    # Get reference condition (uninfected)
    ref_condition = None
    for condition in ['DOX', 'uninfected']:
        if condition in tad_boundaries:
            ref_condition = condition
            break
    
    if ref_condition is None:
        print("Warning: No reference condition found")
        return pd.DataFrame()
    
    ref_boundaries = tad_boundaries[ref_condition]
    comparison_results = []
    
    # Compare each other condition to reference
    for infected_condition in tad_boundaries.keys():
        if infected_condition == ref_condition:
            continue
        
        infected_boundaries = tad_boundaries[infected_condition]
        
        print(f"  Comparing {ref_condition} ({len(ref_boundaries)}) vs {infected_condition} ({len(infected_boundaries)})")
        
        # Find overlapping boundaries
        overlaps = []
        
        for _, ref_boundary in ref_boundaries.iterrows():
            # Find overlapping boundaries in infected condition
            chrom_boundaries = infected_boundaries[
                infected_boundaries['chrom'] == ref_boundary['chrom']
            ]
            
            for _, inf_boundary in chrom_boundaries.iterrows():
                # Check for overlap (50% reciprocal overlap)
                overlap_start = max(ref_boundary['start'], inf_boundary['start'])
                overlap_end = min(ref_boundary['end'], inf_boundary['end'])
                
                if overlap_start < overlap_end:
                    overlap_length = overlap_end - overlap_start
                    ref_length = ref_boundary['end'] - ref_boundary['start']
                    inf_length = inf_boundary['end'] - inf_boundary['start']
                    
                    if (overlap_length / ref_length >= 0.5 and 
                        overlap_length / inf_length >= 0.5):
                        
                        # Calculate change in boundary strength
                        strength_change = (inf_boundary['boundary_strength'] - 
                                         ref_boundary['boundary_strength'])
                        
                        overlaps.append({
                            'chrom': ref_boundary['chrom'],
                            'ref_start': ref_boundary['start'],
                            'ref_end': ref_boundary['end'],
                            'inf_start': inf_boundary['start'],
                            'inf_end': inf_boundary['end'],
                            'ref_strength': ref_boundary['boundary_strength'],
                            'inf_strength': inf_boundary['boundary_strength'],
                            'strength_change': strength_change,
                            'log2_fold_change': np.log2(
                                (inf_boundary['boundary_strength'] + 1) / 
                                (ref_boundary['boundary_strength'] + 1)
                            )
                        })
        
        if not overlaps:
            print(f"    No overlapping boundaries found")
            continue
        
        overlaps_df = pd.DataFrame(overlaps)
        
        # Statistical testing using permutation
        n_permutations = 1000
        null_changes = []
        
        for _ in range(n_permutations):
            # Randomly shuffle the strength values
            perm_strengths = np.random.permutation(overlaps_df['inf_strength'])
            null_change = np.mean(perm_strengths - overlaps_df['ref_strength'])
            null_changes.append(null_change)
        
        null_changes = np.array(null_changes)
        observed_change = overlaps_df['strength_change'].mean()
        
        # Calculate p-value
        p_value = np.sum(np.abs(null_changes) >= np.abs(observed_change)) / n_permutations
        
        # Per-boundary p-values
        boundary_p_values = []
        for _, row in overlaps_df.iterrows():
            # Use null model to test individual boundary changes
            null_dist = np.random.normal(0, null_changes.std(), 1000)
            p_val = np.sum(np.abs(null_dist) >= np.abs(row['strength_change'])) / 1000
            boundary_p_values.append(p_val)
        
        overlaps_df['p_value'] = boundary_p_values
        
        # FDR correction
        if len(boundary_p_values) > 0:
            _, overlaps_df['fdr'], _, _ = multipletests(boundary_p_values, method='fdr_bh')
            overlaps_df['significant'] = overlaps_df['fdr'] < fdr_threshold
        else:
            overlaps_df['fdr'] = []
            overlaps_df['significant'] = []
        
        overlaps_df['comparison'] = f"{ref_condition}_vs_{infected_condition}"
        overlaps_df['global_p_value'] = p_value
        
        comparison_results.append(overlaps_df)
        
        print(f"    Found {len(overlaps_df)} overlapping boundaries, {overlaps_df['significant'].sum() if len(overlaps_df) > 0 else 0} significant")
    
    if comparison_results:
        return pd.concat(comparison_results, ignore_index=True)
    else:
        return pd.DataFrame()

def calculate_compartments_all_conditions(mcool_files, conditions, chromosomes, 
                                         resolution=64000, n_eigs=3):
    """
    Calculate A/B compartments for all conditions using cooltools.
    """
    all_compartments = {}
    
    for condition, mcool_file in zip(conditions, mcool_files):
        print(f"\nCalculating compartments for {condition}...")
        
        # Get the closest available resolution
        actual_resolution = get_closest_resolution(mcool_file, resolution)
        if actual_resolution is None:
            print(f"Error: Could not determine resolution for {mcool_file}")
            continue
        
        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_resolution}")
        
        # Create a proper bioframe viewframe
        from cooltools.lib.common import make_cooler_view
        view_df = make_cooler_view(clr)
        
        # Filter to only include the chromosomes we want (exclude Y and Mt)
        filtered_chroms = [c for c in chromosomes if c not in ['Y', 'Mt']]
        view_df = view_df[view_df['chrom'].isin(filtered_chroms)]
        
        if view_df.empty:
            print(f"Warning: No valid chromosomes found for {condition}")
            continue
        
        try:
            # Calculate eigenvectors using eigs_cis
            print(f"  Calculating eigenvectors for {condition}...")
            eigvals, eigvecs = eigs_cis(
                clr, 
                view_df=view_df,
                n_eigs=n_eigs,
                clr_weight_name='weight',
                ignore_diags=2,
                clip_percentile=99.9
            )
            
            # eigvecs is a DataFrame with bins and eigenvector columns
            if not eigvecs.empty:
                # Add condition and resolution info
                eigvecs['condition'] = condition
                eigvecs['resolution'] = actual_resolution
                
                # Add compartment labels based on first eigenvector (E1)
                if 'E1' in eigvecs.columns:
                    eigvecs['compartment'] = np.where(eigvecs['E1'] > 0, 'A', 'B')
                    all_compartments[condition] = eigvecs
                else:
                    print(f"Warning: No E1 column found in eigenvector results for {condition}")
                    
        except Exception as e:
            print(f"Warning: Could not calculate compartments for {condition}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    return all_compartments

def compare_compartments_to_null(compartment_data, null_model, conditions, fdr_threshold=0.05):
    """
    Compare compartment changes between conditions with FDR correction.
    """
    print("\nComparing compartment changes...")
    
    # Find reference condition
    ref_condition = None
    for cond in ['DOX', 'uninfected']:
        if cond in compartment_data:
            ref_condition = cond
            break
    
    if ref_condition is None:
        print("Warning: No reference condition found in compartment data")
        return pd.DataFrame()
    
    uninfected = compartment_data[ref_condition]
    comparison_results = []
    
    for infected_condition in conditions:
        if infected_condition == ref_condition or infected_condition not in compartment_data:
            continue
            
        infected = compartment_data[infected_condition]
        
        # Merge on genomic bins
        merged = pd.merge(
            uninfected[['chrom', 'start', 'end', 'E1', 'compartment']],
            infected[['chrom', 'start', 'end', 'E1', 'compartment']],
            on=['chrom', 'start', 'end'],
            suffixes=('_uninf', '_inf')
        )
        
        if len(merged) == 0:
            print(f"Warning: No matching bins found between {ref_condition} and {infected_condition}")
            continue
        
        # Identify compartment switches
        merged['switch'] = merged['compartment_uninf'] != merged['compartment_inf']
        merged['E1_diff'] = merged['E1_inf'] - merged['E1_uninf']
        
        # Statistical test for E1 differences using null model
        null_mean = null_model['logFC'].mean() if 'logFC' in null_model.columns else 0
        null_std = null_model['logFC'].std() if 'logFC' in null_model.columns else 1
        
        # Per-bin p-values based on E1 differences
        bin_p_values = []
        for _, row in merged.iterrows():
            # Test if E1 change is significant compared to null
            z_score = (row['E1_diff'] - null_mean) / null_std if null_std > 0 else 0
            p_val = 2 * (1 - stats.norm.cdf(abs(z_score)))  # Two-tailed test
            bin_p_values.append(p_val)
        
        merged['p_value'] = bin_p_values
        
        # FDR correction
        _, merged['fdr'], _, _ = multipletests(merged['p_value'], method='fdr_bh')
        merged['significant'] = merged['fdr'] < fdr_threshold
        merged['comparison'] = f"{ref_condition}_vs_{infected_condition}"
        
        # Overall switching statistics
        observed_switches = merged['switch'].sum()
        total_bins = len(merged)
        switch_rate = observed_switches / total_bins
        
        # Test switch rate against null expectation (25% random switching)
        expected_switches = total_bins * 0.25
        switch_p_value = stats.binom_test(observed_switches, total_bins, 0.25)
        
        merged['switch_rate'] = switch_rate
        merged['switch_p_value'] = switch_p_value
        
        comparison_results.append(merged)
        
        print(f"  {infected_condition}: {observed_switches}/{total_bins} switches ({switch_rate:.2%}), p={switch_p_value:.3f}")
    
    if comparison_results:
        return pd.concat(comparison_results, ignore_index=True)
    else:
        return pd.DataFrame()

def call_loops_all_conditions(mcool_files, conditions, chromosomes, 
                             resolution=8000, fdr_threshold=0.1):
    """
    Call loops for all conditions using cooltools dots.
    """
    all_loops = {}
    
    for condition, mcool_file in zip(conditions, mcool_files):
        print(f"\nCalling loops for {condition}...")
        
        # Get the closest available resolution
        actual_resolution = get_closest_resolution(mcool_file, resolution)
        if actual_resolution is None:
            print(f"Error: Could not determine resolution for {mcool_file}")
            continue
        
        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_resolution}")
        
        try:
            # Calculate expected
            print(f"Calculating expected for {condition}...")
            expected_df = expected_cis(clr, nproc=4)
            
            # Call dots with correct parameters
            print(f"Calling dots for {condition}...")
            dot_calls = dots(
                clr,
                expected_df,
                expected_value_col="balanced.avg",
                clr_weight_name="weight",
                max_loci_separation=2000000,
                lambda_bin_fdr=fdr_threshold,
                clustering_radius=20000,
                nproc=4
            )
            
            dot_calls['condition'] = condition
            dot_calls['resolution'] = actual_resolution
            all_loops[condition] = dot_calls
            
        except Exception as e:
            print(f"Warning: Could not call loops for {condition}: {e}")
            # Try a simpler approach without expected
            try:
                print(f"  Trying simplified loop calling for {condition}...")
                dot_calls = dots(
                    clr,
                    expected_df,  # Still try with expected first
                    max_loci_separation=2000000,
                    lambda_bin_fdr=fdr_threshold,
                    clustering_radius=20000,
                    nproc=1  # Use single core for fallback
                )
                dot_calls['condition'] = condition
                dot_calls['resolution'] = actual_resolution
                all_loops[condition] = dot_calls
            except Exception as e2:
                print(f"  Simplified approach also failed: {e2}")
                continue
    
    return all_loops

def compare_loops_to_null(loop_data, null_model, conditions, fdr_threshold=0.05):
    """
    Compare loop calls between conditions and test for significance.
    """
    print("\nComparing loops between conditions...")
    
    # Find reference condition
    ref_condition = None
    for cond in ['DOX', 'uninfected']:
        if cond in loop_data:
            ref_condition = cond
            break
    
    if ref_condition is None:
        print("Warning: No reference condition found in loop data")
        print(f"Available conditions: {list(loop_data.keys())}")
        return pd.DataFrame()
    
    uninfected_loops = loop_data[ref_condition]
    comparison_results = []
    
    for infected_condition in conditions:
        if infected_condition == ref_condition or infected_condition not in loop_data:
            continue
            
        infected_loops = loop_data[infected_condition]
        
        # Find differential loops
        def create_loop_set(df):
            """Create a set of loop coordinates with some tolerance"""
            loop_set = set()
            for _, row in df.iterrows():
                # Round coordinates to nearest 10kb for matching
                chrom1 = row['chrom1']
                start1 = round(row['start1'] / 10000) * 10000
                chrom2 = row['chrom2'] 
                start2 = round(row['start2'] / 10000) * 10000
                loop_set.add((chrom1, start1, chrom2, start2))
            return loop_set
        
        print(f"  Comparing {len(uninfected_loops)} uninfected loops vs {len(infected_loops)} infected loops")
        
        uninf_set = create_loop_set(uninfected_loops)
        inf_set = create_loop_set(infected_loops)
        
        # Find gained and lost loops
        gained = inf_set - uninf_set
        lost = uninf_set - inf_set
        
        print(f"  Found {len(gained)} gained loops and {len(lost)} lost loops")
        
        # Statistical test using null model
        # Use the null model's logFC distribution to set significance thresholds
        if 'logFC' in null_model.columns:
            null_mean = null_model['logFC'].mean()
            null_std = null_model['logFC'].std()
        else:
            null_mean = 0
            null_std = 1
        
        # Simple permutation test for loop differences
        n_permutations = 1000
        null_gained = []
        null_lost = []
        
        # Combine all loops for permutation
        all_loops_list = list(uninf_set | inf_set)
        n_uninf = len(uninf_set)
        n_total = len(all_loops_list)
        
        if n_total == 0:
            print(f"  No loops found for comparison in {infected_condition}")
            continue
        
        for _ in range(n_permutations):
            # Randomly assign loops to conditions
            perm_indices = np.random.permutation(n_total)
            perm_uninf_indices = perm_indices[:min(n_uninf, n_total)]
            
            perm_uninf = set(all_loops_list[i] for i in perm_uninf_indices)
            perm_inf = set(all_loops_list) - perm_uninf
            
            null_gained.append(len(perm_inf - perm_uninf))
            null_lost.append(len(perm_uninf - perm_inf))
        
        # Calculate p-values
        if len(null_gained) > 0:
            p_gained = np.sum(np.array(null_gained) >= len(gained)) / n_permutations
            p_lost = np.sum(np.array(null_lost) >= len(lost)) / n_permutations
        else:
            p_gained = 1.0
            p_lost = 1.0
        
        # Apply multiple testing correction
        n_comparisons = len(conditions) - 1
        fdr_gained = min(p_gained * n_comparisons, 1.0)
        fdr_lost = min(p_lost * n_comparisons, 1.0)
        
        comparison_results.append({
            'comparison': f"{ref_condition}_vs_{infected_condition}",
            'gained_loops': len(gained),
            'lost_loops': len(lost),
            'total_uninf_loops': len(uninf_set),
            'total_inf_loops': len(inf_set),
            'p_value_gained': p_gained,
            'p_value_lost': p_lost,
            'fdr_gained': fdr_gained,
            'fdr_lost': fdr_lost,
            'significant_gained': fdr_gained < fdr_threshold,
            'significant_lost': fdr_lost < fdr_threshold
        })
        
        print(f"  {infected_condition}: gained={len(gained)} (p={p_gained:.3f}), lost={len(lost)} (p={p_lost:.3f})")
    
    return pd.DataFrame(comparison_results)

def create_summary_plots(tad_comparison, compartment_comp, loop_comp, output_prefix):
    """
    Create summary visualizations for all comparisons.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot 1: TAD boundary changes
    ax = axes[0, 0]
    if not tad_comparison.empty and 'comparison' in tad_comparison.columns:
        try:
            boundary_stats = tad_comparison.groupby('comparison')['significant'].sum()
            if len(boundary_stats) > 0:
                ax.bar(range(len(boundary_stats)), boundary_stats.values)
                ax.set_xticks(range(len(boundary_stats)))
                ax.set_xticklabels([x.replace('_vs_', ' vs ') for x in boundary_stats.index], rotation=45)
        except Exception as e:
            print(f"Warning: Could not plot TAD boundary changes: {e}")
    ax.set_ylabel('Number of Significant Changes')
    ax.set_title('TAD Boundary Changes')
    
    # Plot 2: TAD boundary strength changes
    ax = axes[0, 1]
    if not tad_comparison.empty and 'log2_fold_change' in tad_comparison.columns:
        try:
            for comparison in tad_comparison['comparison'].unique():
                data = tad_comparison[tad_comparison['comparison'] == comparison]
                sig_data = data[data['significant']]
                if len(sig_data) > 0:
                    ax.hist(sig_data['log2_fold_change'], alpha=0.5, 
                           label=comparison.replace('_vs_', ' vs '), bins=20)
            ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
            ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot boundary strength changes: {e}")
    ax.set_xlabel('Log2 Fold Change in Boundary Strength')
    ax.set_ylabel('Count')
    ax.set_title('TAD Boundary Strength Changes')
    
    # Plot 3: Compartment switches
    ax = axes[0, 2]
    if not compartment_comp.empty and 'comparison' in compartment_comp.columns:
        try:
            comparisons = compartment_comp['comparison'].unique()
            switch_rates = []
            sig_switch_rates = []
            
            for comparison in comparisons:
                data = compartment_comp[compartment_comp['comparison'] == comparison]
                switch_rate = data['switch'].sum() / len(data) * 100 if len(data) > 0 else 0
                sig_switch_rate = data[data['significant']]['switch'].sum() / len(data) * 100 if len(data) > 0 else 0
                switch_rates.append(switch_rate)
                sig_switch_rates.append(sig_switch_rate)
            
            if len(comparisons) > 0:
                x = np.arange(len(comparisons))
                width = 0.35
                ax.bar(x - width/2, switch_rates, width, label='All', alpha=0.7)
                ax.bar(x + width/2, sig_switch_rates, width, label='Significant', alpha=0.7)
                ax.set_xticks(x)
                ax.set_xticklabels([x.replace('_vs_', ' vs ') for x in comparisons], rotation=45)
                ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot compartment switches: {e}")
    ax.set_ylabel('Switch Rate (%)')
    ax.set_title('Compartment Switch Rates')
    
    # Plot 4: E1 value changes
    ax = axes[1, 0]
    if not compartment_comp.empty and 'E1_uninf' in compartment_comp.columns:
        try:
            for comparison in compartment_comp['comparison'].unique():
                data = compartment_comp[compartment_comp['comparison'] == comparison]
                sig_data = data[data['significant']]
                if len(sig_data) > 0:
                    valid_data = sig_data.dropna(subset=['E1_uninf', 'E1_inf'])
                    if len(valid_data) > 0:
                        ax.scatter(valid_data['E1_uninf'], valid_data['E1_inf'], 
                                  alpha=0.5, s=10, label=comparison.replace('_vs_', ' vs '))
            ax.plot([-1, 1], [-1, 1], 'k--', alpha=0.5)
            ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot E1 changes: {e}")
    ax.set_xlabel('E1 Uninfected')
    ax.set_ylabel('E1 Infected')
    ax.set_title('Compartment Strength Changes')
    
    # Plot 5: Loop changes
    ax = axes[1, 1]
    if not loop_comp.empty and 'gained_loops' in loop_comp.columns:
        try:
            x = np.arange(len(loop_comp))
            width = 0.35
            ax.bar(x - width/2, loop_comp['gained_loops'], width, label='Gained')
            ax.bar(x + width/2, loop_comp['lost_loops'], width, label='Lost')
            ax.set_xticks(x)
            ax.set_xticklabels([x.replace('_vs_', ' vs ') for x in loop_comp['comparison']], rotation=45)
            ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot loop changes: {e}")
    ax.set_ylabel('Number of Loops')
    ax.set_title('Loop Changes')
    
    # Plot 6: Summary statistics
    ax = axes[1, 2]
    summary_data = []
    
    # Get all unique comparisons
    all_comparisons = set()
    if not tad_comparison.empty and 'comparison' in tad_comparison.columns:
        all_comparisons.update(tad_comparison['comparison'].unique())
    if not compartment_comp.empty and 'comparison' in compartment_comp.columns:
        all_comparisons.update(compartment_comp['comparison'].unique())
    if not loop_comp.empty and 'comparison' in loop_comp.columns:
        all_comparisons.update(loop_comp['comparison'].unique())
    
    for comp in all_comparisons:
        tad_sig = 0
        comp_sig = 0
        loop_count = 0
        
        try:
            if not tad_comparison.empty and 'comparison' in tad_comparison.columns:
                tad_sig = tad_comparison[tad_comparison['comparison'] == comp]['significant'].sum()
            
            if not compartment_comp.empty and 'comparison' in compartment_comp.columns:
                comp_sig = compartment_comp[compartment_comp['comparison'] == comp]['significant'].sum()
            
            if not loop_comp.empty and 'comparison' in loop_comp.columns:
                loop_data = loop_comp[loop_comp['comparison'] == comp]
                if len(loop_data) > 0:
                    loop_count = loop_data['gained_loops'].values[0] + loop_data['lost_loops'].values[0]
        except Exception as e:
            print(f"Warning: Could not calculate summary for {comp}: {e}")
        
        # Extract comparison name
        comp_name = comp.replace('DOX_vs_', '').replace('uninfected_vs_', '')
        
        summary_data.append({
            'Comparison': comp_name,
            'TAD boundaries': tad_sig,
            'Compartments': comp_sig,
            'Loops': loop_count
        })
    
    summary_df = pd.DataFrame(summary_data)
    if not summary_df.empty:
        try:
            summary_df.set_index('Comparison').plot(kind='bar', ax=ax)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        except Exception as e:
            print(f"Warning: Could not plot summary: {e}")
    ax.set_ylabel('Count')
    ax.set_title('Summary of Significant Changes')
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_chromatin_structure_comparison.pdf", dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Compare chromatin structure using existing diffHic results')
    parser.add_argument('--mcool_files', nargs='+', required=True,
                       help='Micro-C mcool files (must match number of conditions)')
    parser.add_argument('--conditions', nargs='+', required=True,
                       help='Condition names (must match number of mcool files)')
    parser.add_argument('--diffhic_results', required=True,
                       help='Combined diffHic results file (all_results_combined.csv)')
    parser.add_argument('--null_model', required=True,
                       help='Null model results file (null_model_results.csv)')
    parser.add_argument('--chromosomes', nargs='+', 
                       default=['2L', '2R', '3L', '3R', '4', 'X'],
                       help='Chromosomes to analyze')
    parser.add_argument('--resolution_compartment', type=int, default=50000,
                       help='Resolution for compartment analysis')
    parser.add_argument('--resolution_loop', type=int, default=5000,
                       help='Resolution for loop calling')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significance')
    parser.add_argument('--tad_window_size', type=int, default=50000,
                       help='Window size for TAD boundary identification')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    # Validation
    if len(args.mcool_files) != len(args.conditions):
        raise ValueError(f"Number of mcool files ({len(args.mcool_files)}) must match number of conditions ({len(args.conditions)})")
    
    print(f"Analyzing {len(args.conditions)} conditions: {', '.join(args.conditions)}")
    print(f"Using cooltools version: {cooltools.__version__}")
    
    # Load diffHic results and null model
    all_results, significant_results = load_diffhic_results(args.diffhic_results)
    null_model = load_null_model(args.null_model)
    
    # Identify TAD boundaries from differential interactions
    print("\n" + "="*60)
    print("IDENTIFYING TAD BOUNDARIES FROM DIFFERENTIAL INTERACTIONS")
    print("="*60)
    
    tad_boundaries = identify_tad_boundaries_from_interactions(
        significant_results, args.conditions, window_size=args.tad_window_size
    )
    
    # Compare TAD boundaries between conditions
    tad_comparison = compare_tad_boundaries(
        tad_boundaries, null_model, args.fdr_threshold
    )
    
    # Calculate compartments using cooltools (still useful for this analysis)
    print("\n" + "="*60)
    print("CALCULATING COMPARTMENTS")
    print("="*60)
    
    compartment_data = calculate_compartments_all_conditions(
        args.mcool_files, args.conditions, args.chromosomes,
        resolution=args.resolution_compartment
    )
    
    # Compare compartments
    compartment_comparison = compare_compartments_to_null(
        compartment_data, null_model, args.conditions, args.fdr_threshold
    )
    
    # Call loops using cooltools
    print("\n" + "="*60)
    print("CALLING LOOPS")
    print("="*60)
    
    loop_data = call_loops_all_conditions(
        args.mcool_files, args.conditions, args.chromosomes,
        resolution=args.resolution_loop
    )
    
    # Compare loops
    loop_comparison = compare_loops_to_null(loop_data, null_model, args.conditions, args.fdr_threshold)
    
    # Create visualizations
    print("\n" + "="*60)
    print("CREATING VISUALIZATIONS")
    print("="*60)
    
    create_summary_plots(
        tad_comparison, compartment_comparison, loop_comparison,
        args.output_prefix
    )
    
    # Save results
    print("\n" + "="*60)
    print("SAVING RESULTS")
    print("="*60)
    
    if not tad_comparison.empty:
        tad_comparison.to_csv(
            f"{args.output_prefix}_tad_boundary_comparison.tsv", sep='\t', index=False
        )
        print(f"Saved TAD boundary comparison: {args.output_prefix}_tad_boundary_comparison.tsv")
    
    if not compartment_comparison.empty:
        compartment_comparison.to_csv(
            f"{args.output_prefix}_compartment_comparison.tsv", sep='\t', index=False
        )
        print(f"Saved compartment comparison: {args.output_prefix}_compartment_comparison.tsv")
    
    if not loop_comparison.empty:
        loop_comparison.to_csv(
            f"{args.output_prefix}_loop_comparison.tsv", sep='\t', index=False
        )
        print(f"Saved loop comparison: {args.output_prefix}_loop_comparison.tsv")
    
    # Save TAD boundaries for each condition as BED files
    for condition, boundaries in tad_boundaries.items():
        bed_file = f"{args.output_prefix}_{condition}_tad_boundaries.bed"
        boundaries[['chrom', 'start', 'end']].to_csv(
            bed_file, sep='\t', index=False, header=False
        )
        print(f"Saved {condition} TAD boundaries: {bed_file}")
    
    # Create summary statistics
    summary_stats = {
        'total_interactions_analyzed': len(all_results),
        'significant_interactions': len(significant_results),
        'tad_boundaries_identified': sum(len(df) for df in tad_boundaries.values()),
        'significant_tad_changes': tad_comparison['significant'].sum() if not tad_comparison.empty else 0,
        'significant_compartment_changes': compartment_comparison['significant'].sum() if not compartment_comparison.empty else 0,
        'total_loop_changes': loop_comparison[['gained_loops', 'lost_loops']].sum().sum() if not loop_comparison.empty else 0
    }
    
    # Save summary
    with open(f"{args.output_prefix}_analysis_summary.txt", 'w') as f:
        f.write("CHROMATIN STRUCTURE ANALYSIS SUMMARY\n")
        f.write("="*50 + "\n\n")
        f.write(f"Conditions analyzed: {', '.join(args.conditions)}\n\n")
        f.write(f"Total interactions analyzed: {summary_stats['total_interactions_analyzed']}\n")
        f.write(f"Significant differential interactions: {summary_stats['significant_interactions']}\n")
        f.write(f"TAD boundaries identified: {summary_stats['tad_boundaries_identified']}\n")
        f.write(f"Significant TAD boundary changes: {summary_stats['significant_tad_changes']}\n")
        f.write(f"Significant compartment changes: {summary_stats['significant_compartment_changes']}\n")
        f.write(f"Total loop changes: {summary_stats['total_loop_changes']}\n\n")
        
        # Add per-condition breakdown
        f.write("TAD BOUNDARIES PER CONDITION:\n")
        for condition, boundaries in tad_boundaries.items():
            f.write(f"  {condition}: {len(boundaries)} boundaries\n")
        
        f.write("\nSIGNIFICANT CHANGES PER COMPARISON:\n")
        if not tad_comparison.empty:
            for comparison in tad_comparison['comparison'].unique():
                sig_count = tad_comparison[
                    tad_comparison['comparison'] == comparison
                ]['significant'].sum()
                f.write(f"  TAD boundaries {comparison}: {sig_count}\n")
        
        if not compartment_comparison.empty:
            for comparison in compartment_comparison['comparison'].unique():
                sig_count = compartment_comparison[
                    compartment_comparison['comparison'] == comparison
                ]['significant'].sum()
                switch_rate = compartment_comparison[
                    compartment_comparison['comparison'] == comparison
                ]['switch_rate'].iloc[0] if len(compartment_comparison[
                    compartment_comparison['comparison'] == comparison
                ]) > 0 else 0
                f.write(f"  Compartments {comparison}: {sig_count} significant changes, {switch_rate:.2%} switch rate\n")
        
        if not loop_comparison.empty:
            for _, row in loop_comparison.iterrows():
                f.write(f"  Loops {row['comparison']}: +{row['gained_loops']} gained, -{row['lost_loops']} lost\n")
    
    print(f"\nAnalysis complete! Results saved to {args.output_prefix}_*")
    print(f"Summary saved to: {args.output_prefix}_analysis_summary.txt")

if __name__ == '__main__':
    main()
