#!/usr/bin/env python3
"""
Call TADs, loops, and compartments from Micro-C data using cooltools.
Compare across uninfected and three Wolbachia-infected conditions with null model testing.
"""
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import cooler
import cooltools
# For cooltools 0.7.1, import functions directly
from cooltools import expected_cis
from cooltools import saddle
from cooltools import dots
from cooltools import eigs_cis
from cooltools import insulation  # Import the insulation function directly

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

def load_null_model(null_file):
    """Load null model results from diffHic analysis"""
    print("Loading null model...")
    null_df = pd.read_csv(null_file)
    return null_df


def calculate_insulation_all_conditions(mcool_files, conditions, chromosomes, 
                                       resolution=16000, window_sizes=[160000, 320000, 640000]):
    """
    Calculate insulation scores for all conditions using cooltools.
    """
    all_insulation = {}
    
    for condition, mcool_file in zip(conditions, mcool_files):
        print(f"\nCalculating insulation for {condition}...")
        
        # Get the closest available resolution
        actual_resolution = get_closest_resolution(mcool_file, resolution)
        if actual_resolution is None:
            print(f"Error: Could not determine resolution for {mcool_file}")
            continue
            
        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_resolution}")
        
        # Create a proper bioframe viewframe
        from cooltools.lib.common import make_cooler_view
        view_df = make_cooler_view(clr)
        
        # Filter to only include the chromosomes we want
        view_df = view_df[view_df['chrom'].isin(chromosomes)]
        
        if view_df.empty:
            print(f"Warning: No valid chromosomes found for {condition}")
            continue
        
        try:
            print(f"  Calculating insulation for all window sizes: {window_sizes}")
            
            # Calculate insulation for all window sizes at once
            # The insulation function can handle multiple windows in one call
            ins_table = insulation(
                clr, 
                window_bp=window_sizes,  # Pass the entire list
                view_df=view_df,
                ignore_diags=2,
                min_frac_valid_pixels=0.66,
                clr_weight_name='weight'
            )
            
            ins_table['condition'] = condition
            all_insulation[condition] = ins_table
            
        except Exception as e:
            print(f"Warning: Could not calculate insulation for {condition}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    return all_insulation

def compare_insulation_to_null(insulation_data, null_model, fdr_threshold=0.05):
    """
    Compare insulation scores between conditions and test against null model.
    """
    print("\nComparing insulation scores to null model...")
    
    # Check if we have data
    if not insulation_data:
        print("No insulation data available for comparison")
        return pd.DataFrame()
    
    # Prepare data for comparison
    conditions = list(insulation_data.keys())
    
    if 'DOX' not in conditions:
        print("Error: DOX condition not found in insulation data")
        print(f"Available conditions: {conditions}")
        return pd.DataFrame()
    
    uninfected_data = insulation_data['DOX']
    
    comparison_results = []
    
    # Check what column names are available
    print("Available columns in insulation data:", uninfected_data.columns.tolist())
    
    # Find all insulation score columns (they follow the pattern log2_insulation_score_WINDOW)
    insulation_cols = [col for col in uninfected_data.columns if col.startswith('log2_insulation_score_')]
    
    if not insulation_cols:
        print("Error: Could not find any insulation score columns")
        print("Available columns:", uninfected_data.columns.tolist())
        return pd.DataFrame()
    
    print(f"Found insulation columns: {insulation_cols}")
    
    for infected_condition in ['wMel', 'wRi', 'wWil']:
        if infected_condition not in insulation_data:
            continue
            
        infected_data = insulation_data[infected_condition]
        
        # Process each window size separately
        for insulation_col in insulation_cols:
            # Extract window size from column name
            window_size = insulation_col.split('_')[-1]
            
            # Match bins between conditions
            merged = pd.merge(
                uninfected_data[['chrom', 'start', 'end', insulation_col]],
                infected_data[['chrom', 'start', 'end', insulation_col]],
                on=['chrom', 'start', 'end'],
                suffixes=('_uninf', '_inf')
            )
            
            if len(merged) == 0:
                print(f"Warning: No matching bins found between DOX and {infected_condition} for {insulation_col}")
                continue
            
            # Calculate differences
            merged['insulation_diff'] = merged[f'{insulation_col}_inf'] - merged[f'{insulation_col}_uninf']
            
            # Generate null distribution by permutation
            n_permutations = 1000
            null_diffs = []
            
            for _ in range(n_permutations):
                # Randomly shuffle the labels
                shuffled = merged.copy()
                shuffle_mask = np.random.random(len(shuffled)) < 0.5
                shuffled.loc[shuffle_mask, [f'{insulation_col}_uninf', f'{insulation_col}_inf']] = \
                    shuffled.loc[shuffle_mask, [f'{insulation_col}_inf', f'{insulation_col}_uninf']].values
                
                null_diff = shuffled[f'{insulation_col}_inf'] - shuffled[f'{insulation_col}_uninf']
                null_diffs.append(null_diff.values)
            
            null_diffs = np.array(null_diffs)
            
            # Calculate p-values
            p_values = []
            for i, obs_diff in enumerate(merged['insulation_diff']):
                null_dist = null_diffs[:, i]
                p_val = np.sum(np.abs(null_dist) >= np.abs(obs_diff)) / n_permutations
                p_values.append(p_val)
            
            merged['p_value'] = p_values
            
            # FDR correction
            _, merged['fdr'], _, _ = multipletests(merged['p_value'], method='fdr_bh')
            
            # Mark significant changes
            merged['significant'] = merged['fdr'] < fdr_threshold
            merged['comparison'] = f"DOX_vs_{infected_condition}"
            merged['window_size'] = window_size
            
            comparison_results.append(merged)
    
    if comparison_results:
        return pd.concat(comparison_results, ignore_index=True)
    else:
        return pd.DataFrame()

def call_tads_all_conditions(insulation_data, min_boundary_strength=0.1):
    """
    Call TAD boundaries from insulation scores for all conditions.
    """
    all_boundaries = {}
    
    for condition, ins_data in insulation_data.items():
        print(f"\nCalling TAD boundaries for {condition}...")
        
        boundaries = []
        
        # Find all boundary strength columns (they follow the pattern boundary_strength_WINDOW)
        boundary_cols = [col for col in ins_data.columns if col.startswith('boundary_strength_')]
        insulation_cols = [col for col in ins_data.columns if col.startswith('log2_insulation_score_')]
        
        if not boundary_cols or not insulation_cols:
            print(f"Warning: Could not find boundary or insulation columns for {condition}")
            continue
        
        # Process each window size
        for boundary_col in boundary_cols:
            # Extract window size from column name
            window_size = boundary_col.split('_')[-1]
            insulation_col = f'log2_insulation_score_{window_size}'
            
            if insulation_col not in ins_data.columns:
                print(f"Warning: No matching insulation column for {boundary_col}")
                continue
            
            # Group by chromosome for boundary calling
            for chrom in ins_data['chrom'].unique():
                try:
                    chrom_data = ins_data[ins_data['chrom'] == chrom].copy()
                    
                    # Filter out NaN values
                    valid_data = chrom_data.dropna(subset=[insulation_col, boundary_col])
                    
                    if len(valid_data) < 3:
                        continue
                    
                    # Find boundaries using the boundary_strength column
                    # Boundaries are typically defined where boundary_strength > threshold
                    is_boundary_col = f'is_boundary_{window_size}'
                    if is_boundary_col in valid_data.columns:
                        # Use the pre-computed boundary calls
                        boundaries_mask = valid_data[is_boundary_col] == True
                    else:
                        # Fallback: use boundary strength threshold
                        boundaries_mask = valid_data[boundary_col] > min_boundary_strength
                    
                    if np.any(boundaries_mask):
                        boundary_data = valid_data[boundaries_mask].copy()
                        boundary_data['condition'] = condition
                        boundary_data['window_size'] = window_size
                        boundaries.append(boundary_data)
                        
                except Exception as e:
                    print(f"Warning: Could not find boundaries for {condition} {chrom} window {window_size}: {e}")
                    continue
        
        if boundaries:
            all_boundaries[condition] = pd.concat(boundaries, ignore_index=True)
    
    return all_boundaries

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

def compare_compartments_to_null(compartment_data, null_model, fdr_threshold=0.05):
    """
    Compare compartment changes between conditions with FDR correction.
    """
    print("\nComparing compartment changes...")
    
    if 'DOX' not in compartment_data:
        print("Warning: DOX condition not found in compartment data")
        return pd.DataFrame()
    
    uninfected = compartment_data['DOX']
    comparison_results = []
    
    for infected_condition in ['wMel', 'wRi', 'wWil']:
        if infected_condition not in compartment_data:
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
            print(f"Warning: No matching bins found between DOX and {infected_condition}")
            continue
        
        # Identify compartment switches
        merged['switch'] = merged['compartment_uninf'] != merged['compartment_inf']
        merged['E1_diff'] = merged['E1_inf'] - merged['E1_uninf']
        
        # Statistical test for E1 differences
        # Use permutation test
        n_permutations = 1000
        null_switches = []
        
        for _ in range(n_permutations):
            # Randomly assign compartments
            perm_data = merged.copy()
            shuffle_idx = np.random.permutation(len(perm_data))
            perm_data['compartment_inf'] = perm_data['compartment_inf'].iloc[shuffle_idx].values
            perm_switches = (perm_data['compartment_uninf'] != perm_data['compartment_inf']).sum()
            null_switches.append(perm_switches)
        
        # Calculate p-value for overall switching rate
        observed_switches = merged['switch'].sum()
        p_value = np.sum(null_switches >= observed_switches) / n_permutations
        
        # Per-bin p-values based on E1 differences
        bin_p_values = []
        for _, row in merged.iterrows():
            # Test if E1 change is significant
            null_dist = np.random.normal(0, np.std(merged['E1_diff']), n_permutations)
            p_val = np.sum(np.abs(null_dist) >= np.abs(row['E1_diff'])) / n_permutations
            bin_p_values.append(p_val)
        
        merged['p_value'] = bin_p_values
        
        # FDR correction
        _, merged['fdr'], _, _ = multipletests(merged['p_value'], method='fdr_bh')
        merged['significant'] = merged['fdr'] < fdr_threshold
        merged['comparison'] = f"DOX_vs_{infected_condition}"
        merged['global_switch_pvalue'] = p_value
        
        comparison_results.append(merged)
    
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


def create_summary_plots(insulation_comp, compartment_comp, loop_comp, output_prefix):
    """
    Create summary visualizations for all comparisons.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot 1: Significant insulation changes
    ax = axes[0, 0]
    if not insulation_comp.empty and 'comparison' in insulation_comp.columns and 'insulation_diff' in insulation_comp.columns:
        for comparison in insulation_comp['comparison'].unique():
            data = insulation_comp[insulation_comp['comparison'] == comparison]
            sig_data = data[data['significant']]
            if len(sig_data) > 0 and not sig_data['insulation_diff'].isna().all():
                # Remove NaN values before plotting
                valid_diff = sig_data['insulation_diff'].dropna()
                if len(valid_diff) > 0:
                    ax.hist(valid_diff, alpha=0.5, label=comparison, bins=30)
    ax.set_xlabel('Insulation Score Change')
    ax.set_ylabel('Count')
    ax.set_title('Significant Insulation Changes')
    ax.legend()
    
    # Plot 2: TAD boundary changes
    ax = axes[0, 1]
    if not insulation_comp.empty and 'comparison' in insulation_comp.columns:
        try:
            boundary_stats = insulation_comp.groupby('comparison')['significant'].sum()
            if len(boundary_stats) > 0:
                ax.bar(range(len(boundary_stats)), boundary_stats.values)
                ax.set_xticks(range(len(boundary_stats)))
                ax.set_xticklabels(boundary_stats.index, rotation=45)
        except Exception as e:
            print(f"Warning: Could not plot boundary changes: {e}")
    ax.set_ylabel('Number of Significant Changes')
    ax.set_title('TAD Boundary Changes')
    
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
                ax.set_xticklabels(comparisons, rotation=45)
                ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot compartment switches: {e}")
    ax.set_ylabel('Switch Rate (%)')
    ax.set_title('Compartment Switch Rates')
    
    # Plot 4: E1 value changes
    ax = axes[1, 0]
    if not compartment_comp.empty and 'comparison' in compartment_comp.columns:
        try:
            for comparison in compartment_comp['comparison'].unique():
                data = compartment_comp[compartment_comp['comparison'] == comparison]
                sig_data = data[data['significant']]
                if len(sig_data) > 0 and 'E1_uninf' in sig_data.columns and 'E1_inf' in sig_data.columns:
                    # Remove NaN values
                    valid_data = sig_data.dropna(subset=['E1_uninf', 'E1_inf'])
                    if len(valid_data) > 0:
                        ax.scatter(valid_data['E1_uninf'], valid_data['E1_inf'], 
                                  alpha=0.5, s=10, label=comparison)
            ax.plot([-1, 1], [-1, 1], 'k--', alpha=0.5)
            ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot E1 changes: {e}")
    ax.set_xlabel('E1 Uninfected')
    ax.set_ylabel('E1 Infected')
    ax.set_title('Compartment Strength Changes')
    
    # Plot 5: Loop changes
    ax = axes[1, 1]
    if not loop_comp.empty and 'gained_loops' in loop_comp.columns and 'lost_loops' in loop_comp.columns:
        try:
            x = np.arange(len(loop_comp))
            width = 0.35
            ax.bar(x - width/2, loop_comp['gained_loops'], width, label='Gained')
            ax.bar(x + width/2, loop_comp['lost_loops'], width, label='Lost')
            ax.set_xticks(x)
            ax.set_xticklabels(loop_comp['comparison'], rotation=45)
            ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot loop changes: {e}")
    ax.set_ylabel('Number of Loops')
    ax.set_title('Loop Changes')
    
    # Plot 6: Summary statistics
    ax = axes[1, 2]
    summary_data = []
    for comp in ['DOX_vs_wMel', 'DOX_vs_wRi', 'DOX_vs_wWil']:
        ins_sig = 0
        comp_sig = 0
        loop_count = 0
        
        try:
            if not insulation_comp.empty and 'comparison' in insulation_comp.columns:
                ins_sig = insulation_comp[insulation_comp['comparison'] == comp]['significant'].sum()
            
            if not compartment_comp.empty and 'comparison' in compartment_comp.columns:
                comp_sig = compartment_comp[compartment_comp['comparison'] == comp]['significant'].sum()
            
            if not loop_comp.empty and 'comparison' in loop_comp.columns:
                loop_data = loop_comp[loop_comp['comparison'] == comp]
                if len(loop_data) > 0:
                    loop_count = loop_data['gained_loops'].values[0] + loop_data['lost_loops'].values[0]
        except Exception as e:
            print(f"Warning: Could not calculate summary for {comp}: {e}")
        
        summary_data.append({
            'Comparison': comp.replace('DOX_vs_', ''),
            'TAD boundaries': ins_sig,
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

def compare_loops_to_null(loop_data, null_model, fdr_threshold=0.05):
    """
    Compare loop calls between conditions and test for significance.
    """
    print("\nComparing loops between conditions...")
    
    if 'DOX' not in loop_data:
        print("Warning: DOX condition not found in loop data")
        print(f"Available conditions: {list(loop_data.keys())}")
        return pd.DataFrame()
    
    uninfected_loops = loop_data['DOX']
    comparison_results = []
    
    for infected_condition in ['wMel', 'wRi', 'wWil']:
        if infected_condition not in loop_data:
            print(f"Warning: {infected_condition} not found in loop data")
            continue
            
        infected_loops = loop_data[infected_condition]
        
        # Find differential loops
        # Create sets of loop coordinates for comparison
        # Use a tolerance for coordinate matching since resolutions might differ slightly
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
        
        # Statistical test using permutation
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
            # Use numpy random indices instead of choice on complex objects
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
        
        # Apply multiple testing correction (Bonferroni for 3 comparisons)
        fdr_gained = min(p_gained * 3, 1.0)
        fdr_lost = min(p_lost * 3, 1.0)
        
        comparison_results.append({
            'comparison': f"DOX_vs_{infected_condition}",
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

def main():
    parser = argparse.ArgumentParser(description='Compare chromatin structure across conditions')
    parser.add_argument('--mcool_files', nargs=4, required=True,
                       help='Micro-C mcool files for DOX, wMel, wRi, wWil')
    parser.add_argument('--conditions', nargs=4, default=['DOX', 'wMel', 'wRi', 'wWil'],
                       help='Condition names')
    parser.add_argument('--diffhic_results', required=True,
                       help='DiffHic results file')
    parser.add_argument('--chromosomes', nargs='+', 
                       default=['2L', '2R', '3L', '3R', '4', 'X'],
                       help='Chromosomes to analyze')
    parser.add_argument('--resolution_tad', type=int, default=10000,
                       help='Resolution for TAD calling')
    parser.add_argument('--resolution_compartment', type=int, default=50000,
                       help='Resolution for compartment analysis')
    parser.add_argument('--resolution_loop', type=int, default=5000,
                       help='Resolution for loop calling')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significance')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    print(f"Using cooltools version: {cooltools.__version__}")
    
    # Load null model
    null_model = load_null_model(args.diffhic_results)
    
    # Calculate insulation scores for all conditions
    insulation_data = calculate_insulation_all_conditions(
        args.mcool_files, args.conditions, args.chromosomes,
        resolution=args.resolution_tad
    )
    
    # Compare insulation scores
    insulation_comparison = compare_insulation_to_null(
        insulation_data, null_model, args.fdr_threshold
    )
    
    # Call TAD boundaries
    tad_boundaries = call_tads_all_conditions(insulation_data)
    
    # Calculate compartments
    compartment_data = calculate_compartments_all_conditions(
        args.mcool_files, args.conditions, args.chromosomes,
        resolution=args.resolution_compartment
    )
    
    # Compare compartments
    compartment_comparison = compare_compartments_to_null(
        compartment_data, null_model, args.fdr_threshold
    )
    
    # Call loops
    loop_data = call_loops_all_conditions(
        args.mcool_files, args.conditions, args.chromosomes,
        resolution=args.resolution_loop
    )
    
    # Compare loops
    loop_comparison = compare_loops_to_null(loop_data, null_model, args.fdr_threshold)
    
    # Create visualizations
    create_summary_plots(
        insulation_comparison, compartment_comparison, loop_comparison,
        args.output_prefix
    )
    
    # Save results
    if not insulation_comparison.empty:
        insulation_comparison.to_csv(
            f"{args.output_prefix}_insulation_comparison.tsv", sep='\t', index=False
        )
    
    if not compartment_comparison.empty:
        compartment_comparison.to_csv(
            f"{args.output_prefix}_compartment_comparison.tsv", sep='\t', index=False
        )
    
    if not loop_comparison.empty:
        loop_comparison.to_csv(
            f"{args.output_prefix}_loop_comparison.tsv", sep='\t', index=False
        )
    
    # Save TAD boundaries for each condition
    for condition, boundaries in tad_boundaries.items():
        boundaries.to_csv(
            f"{args.output_prefix}_{condition}_tad_boundaries.bed",
            sep='\t', index=False
        )
    
    print(f"\nAnalysis complete. Results saved to {args.output_prefix}_*")

if __name__ == '__main__':
    main()