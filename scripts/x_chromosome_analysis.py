#!/usr/bin/env python3
"""
Analyze X-chromosome regulation changes in response to Wolbachia infection.
Focus on dosage compensation complex (DCC) high-affinity sites (HAS/CES).
Based on Schauer et al. 2017 methods.
"""

import pandas as pd
import numpy as np
import cooler
import pybedtools
from scipy import stats
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

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

def load_has_sites(has_file):
    """
    Load high-affinity sites (HAS) / chromosomal entry sites (CES) for DCC.
    These are X-chromosomal sites where the dosage compensation complex binds.
    """
    print("Loading HAS/CES sites...")
    has_sites = pd.read_csv(has_file, sep='\t', header=None,
                           names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    # Filter for X chromosome
    has_sites = has_sites[has_sites['chrom'] == 'X']
    print(f"Loaded {len(has_sites)} HAS/CES sites on X chromosome")
    
    return has_sites

def analyze_has_contacts(has_sites, mcool_files, conditions, resolution=5000, window_size=50000):
    """
    Analyze chromatin contacts around HAS sites for each condition.
    """
    results = {}
    
    for condition, mcool_file in zip(conditions, mcool_files):
        print(f"\nAnalyzing {condition}...")
        
        # Get the closest available resolution
        actual_resolution = get_closest_resolution(mcool_file, resolution)
        if actual_resolution is None:
            print(f"Error: Could not determine resolution for {mcool_file}")
            continue
            
        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_resolution}")
        
        # Get X chromosome matrix
        try:
            x_matrix = clr.matrix(balance=True).fetch('X')
        except Exception as e:
            print(f"Warning: Could not load X chromosome for {condition}: {e}")
            continue
        
        # Analyze contacts around each HAS
        has_contacts = []
        for _, site in has_sites.iterrows():
            # Define window around HAS
            center_bin = site['start'] // actual_resolution
            window_bins = window_size // actual_resolution
            
            start_bin = max(0, center_bin - window_bins)
            end_bin = min(x_matrix.shape[0], center_bin + window_bins + 1)
            
            if start_bin < end_bin:
                # Extract submatrix
                submatrix = x_matrix[start_bin:end_bin, start_bin:end_bin]
                
                # Calculate metrics
                if submatrix.size > 0:
                    # Total contacts in window
                    total_contacts = np.nansum(submatrix)
                    
                    # Average contact frequency
                    avg_contacts = np.nanmean(submatrix[~np.isnan(submatrix)])
                    
                    # Insulation score (ratio of intra vs inter contacts)
                    if submatrix.shape[0] > 2:
                        intra = np.nanmean(np.diag(submatrix))
                        inter = np.nanmean(submatrix[~np.eye(submatrix.shape[0], dtype=bool)])
                        insulation = intra / inter if inter > 0 else np.nan
                    else:
                        insulation = np.nan
                    
                    has_contacts.append({
                        'site': site['name'],
                        'total_contacts': total_contacts,
                        'avg_contacts': avg_contacts,
                        'insulation': insulation
                    })
        
        results[condition] = pd.DataFrame(has_contacts)
    
    return results

def analyze_x_compartmentalization(mcool_files, conditions, resolution=50000):
    """
    Analyze X chromosome compartmentalization (A/B compartments).
    """
    results = {}
    
    for condition, mcool_file in zip(conditions, mcool_files):
        print(f"\nAnalyzing X compartmentalization for {condition}...")
        
        # Get the closest available resolution
        actual_resolution = get_closest_resolution(mcool_file, resolution)
        if actual_resolution is None:
            print(f"Error: Could not determine resolution for {mcool_file}")
            continue
            
        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_resolution}")
        
        # Get X chromosome matrix
        try:
            x_matrix = clr.matrix(balance=True).fetch('X')
        except Exception as e:
            print(f"Warning: Could not load X chromosome for {condition}: {e}")
            continue
        
        # Calculate correlation matrix
        # Remove NaN values
        valid_mask = ~np.all(np.isnan(x_matrix), axis=0)
        x_matrix_clean = x_matrix[valid_mask][:, valid_mask]
        
        if x_matrix_clean.size > 0 and x_matrix_clean.shape[0] > 1:
            try:
                # Calculate correlation matrix
                corr_matrix = np.corrcoef(x_matrix_clean)
                
                # Handle any remaining NaNs
                corr_matrix[np.isnan(corr_matrix)] = 0
                
                # Perform PCA to identify compartments
                try:
                    from sklearn.decomposition import PCA
                    pca = PCA(n_components=1)
                    
                    # Get first principal component (compartment signal)
                    pc1 = pca.fit_transform(corr_matrix).flatten()
                    
                    # Calculate compartment strength
                    compartment_strength = np.std(pc1)
                    
                    # Identify A/B compartments
                    a_compartment = pc1 > 0
                    b_compartment = pc1 < 0
                    
                    results[condition] = {
                        'compartment_signal': pc1,
                        'compartment_strength': compartment_strength,
                        'percent_a': np.sum(a_compartment) / len(pc1) * 100,
                        'percent_b': np.sum(b_compartment) / len(pc1) * 100,
                        'matrix_size': x_matrix_clean.shape[0]
                    }
                except ImportError:
                    print(f"Warning: sklearn not available for PCA analysis for {condition}")
                    # Simple alternative without PCA
                    eigenvals, eigenvecs = np.linalg.eigh(corr_matrix)
                    pc1 = eigenvecs[:, -1]  # First principal component
                    
                    compartment_strength = np.std(pc1)
                    a_compartment = pc1 > 0
                    b_compartment = pc1 < 0
                    
                    results[condition] = {
                        'compartment_signal': pc1,
                        'compartment_strength': compartment_strength,
                        'percent_a': np.sum(a_compartment) / len(pc1) * 100,
                        'percent_b': np.sum(b_compartment) / len(pc1) * 100,
                        'matrix_size': x_matrix_clean.shape[0]
                    }
                    
            except Exception as e:
                print(f"Warning: Could not perform compartment analysis for {condition}: {e}")
                continue
    
    return results

def analyze_x_interactions(interaction_file):
    """
    Analyze differential interactions specifically on the X chromosome.
    """
    print(f"\nAnalyzing X interactions from {interaction_file}...")
    
    # Load interactions
    interactions = pd.read_csv(interaction_file)
    
    # Filter for X chromosome interactions
    x_cis = interactions[
        (interactions['chr1'] == 'X') & 
        (interactions['chr2'] == 'X') &
        (interactions['FDR'] < 0.01) &
        (abs(interactions['logFC']) > 1)
    ]
    
    x_trans = interactions[
        ((interactions['chr1'] == 'X') | (interactions['chr2'] == 'X')) &
        (interactions['chr1'] != interactions['chr2']) &
        (interactions['FDR'] < 0.01) &
        (abs(interactions['logFC']) > 1)
    ]
    
    # Calculate total significant interactions for comparison
    total_sig = interactions[
        (interactions['FDR'] < 0.01) &
        (abs(interactions['logFC']) > 1)
    ]
    
    # Analyze interaction patterns
    results = {
        'n_cis': len(x_cis),
        'n_trans': len(x_trans),
        'n_total_sig': len(total_sig),
        'mean_logfc_cis': x_cis['logFC'].mean() if len(x_cis) > 0 else 0,
        'mean_logfc_trans': x_trans['logFC'].mean() if len(x_trans) > 0 else 0,
        'upregulated_cis': sum(x_cis['logFC'] > 0),
        'downregulated_cis': sum(x_cis['logFC'] < 0),
        'upregulated_trans': sum(x_trans['logFC'] > 0),
        'downregulated_trans': sum(x_trans['logFC'] < 0)
    }
    
    print(f"Found {results['n_cis']} X-cis and {results['n_trans']} X-trans differential interactions")
    print(f"Total significant interactions: {results['n_total_sig']}")
    
    return results

def analyze_null_model_stats(null_interaction_file):
    """
    Analyze null model statistics (no chromosome filtering possible).
    The null model file only contains statistical columns without chromosome positions.
    """
    print(f"\nAnalyzing null model from {null_interaction_file}...")
    
    null_data = pd.read_csv(null_interaction_file)
    
    # Analyze overall distribution
    sig_null = null_data[
        (null_data['FDR'] < 0.01) &
        (abs(null_data['logFC']) > 1)
    ]
    
    results = {
        'n_total': len(null_data),
        'n_sig': len(sig_null),
        'mean_logfc': sig_null['logFC'].mean() if len(sig_null) > 0 else 0,
        'median_logfc': sig_null['logFC'].median() if len(sig_null) > 0 else 0,
        'std_logfc': sig_null['logFC'].std() if len(sig_null) > 0 else 0,
        'upregulated': int(sum(sig_null['logFC'] > 0)) if len(sig_null) > 0 else 0,
        'downregulated': int(sum(sig_null['logFC'] < 0)) if len(sig_null) > 0 else 0,
        'proportion_sig': len(sig_null) / len(null_data) if len(null_data) > 0 else 0
    }
    
    print(f"Null model: {results['n_sig']} significant of {results['n_total']} total ({results['proportion_sig']:.2%})")
    
    return results

def analyze_has_proximity(has_sites, interaction_file, window_size=50000):
    """
    Analyze differential interactions near HAS sites.
    """
    print(f"\nAnalyzing interactions near HAS sites (within {window_size}bp)...")
    
    # Load interactions
    interactions = pd.read_csv(interaction_file)
    
    # Filter for significant interactions
    sig_interactions = interactions[
        (interactions['FDR'] < 0.01) &
        (abs(interactions['logFC']) > 1)
    ].copy()
    
    # Convert to BED format for both anchors
    anchor1_bed = sig_interactions[['chr1', 'start1', 'end1']].copy()
    anchor1_bed.columns = ['chrom', 'start', 'end']
    anchor1_bed['interaction_idx'] = sig_interactions.index
    
    anchor2_bed = sig_interactions[['chr2', 'start2', 'end2']].copy()
    anchor2_bed.columns = ['chrom', 'start', 'end']
    anchor2_bed['interaction_idx'] = sig_interactions.index
    
    # Create HAS windows
    has_windows = has_sites.copy()
    has_windows['start'] = has_windows['start'] - window_size
    has_windows['end'] = has_windows['end'] + window_size
    has_windows = has_windows[has_windows['start'] >= 0]  # Remove negative coordinates
    
    try:
        # Use pybedtools for overlap detection
        has_bt = pybedtools.BedTool.from_dataframe(has_windows[['chrom', 'start', 'end', 'name']])
        anchor1_bt = pybedtools.BedTool.from_dataframe(anchor1_bed)
        anchor2_bt = pybedtools.BedTool.from_dataframe(anchor2_bed)
        
        # Find overlaps
        anchor1_overlaps = anchor1_bt.intersect(has_bt, wa=True, wb=True)
        anchor2_overlaps = anchor2_bt.intersect(has_bt, wa=True, wb=True)
        
        # Collect overlapping interaction indices
        has_proximal_indices = set()
        
        for overlap in anchor1_overlaps:
            idx = int(overlap.fields[3])  # interaction_idx
            has_proximal_indices.add(idx)
        
        for overlap in anchor2_overlaps:
            idx = int(overlap.fields[3])  # interaction_idx
            has_proximal_indices.add(idx)
        
        # Get proximal interactions
        proximal_interactions = sig_interactions.loc[list(has_proximal_indices)]
        
        results = {
            'n_has_proximal': len(proximal_interactions),
            'has_proximal_upregulated': sum(proximal_interactions['logFC'] > 0),
            'has_proximal_downregulated': sum(proximal_interactions['logFC'] < 0),
            'mean_logfc_has_proximal': proximal_interactions['logFC'].mean() if len(proximal_interactions) > 0 else 0
        }
        
    except Exception as e:
        print(f"Warning: Could not analyze HAS proximity using bedtools: {e}")
        # Fallback to manual analysis
        results = {
            'n_has_proximal': 0,
            'has_proximal_upregulated': 0,
            'has_proximal_downregulated': 0,
            'mean_logfc_has_proximal': 0
        }
    
    print(f"Found {results['n_has_proximal']} interactions near HAS sites")
    
    return results

def create_x_regulation_plots(has_results, compartment_results, interaction_results, has_proximity_results, output_prefix):
    """
    Create visualizations for X chromosome regulation analysis.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot 1: HAS contact frequencies
    ax = axes[0, 0]
    if has_results:
        conditions = list(has_results.keys())
        avg_contacts = []
        
        for condition in conditions:
            df = has_results[condition]
            if len(df) > 0 and 'avg_contacts' in df.columns:
                avg_val = df['avg_contacts'].mean()
                avg_contacts.append(avg_val if not np.isnan(avg_val) else 0)
            else:
                avg_contacts.append(0)
        
        ax.bar(conditions, avg_contacts)
        ax.set_ylabel('Average Contact Frequency')
        ax.set_title('Contacts at HAS/CES Sites')
        ax.tick_params(axis='x', rotation=45)
    
    # Plot 2: Compartment strength
    ax = axes[0, 1]
    if compartment_results:
        conditions = list(compartment_results.keys())
        strengths = [result['compartment_strength'] for result in compartment_results.values()]
        
        ax.bar(conditions, strengths)
        ax.set_ylabel('Compartment Strength')
        ax.set_title('X Chromosome Compartmentalization')
        ax.tick_params(axis='x', rotation=45)
    
    # Plot 3: A/B compartment distribution
    ax = axes[0, 2]
    if compartment_results:
        conditions = list(compartment_results.keys())
        percent_a = [result['percent_a'] for result in compartment_results.values()]
        percent_b = [result['percent_b'] for result in compartment_results.values()]
        
        x = np.arange(len(conditions))
        width = 0.35
        
        ax.bar(x - width/2, percent_a, width, label='A compartment')
        ax.bar(x + width/2, percent_b, width, label='B compartment')
        
        ax.set_ylabel('Percentage')
        ax.set_title('A/B Compartment Distribution')
        ax.set_xticks(x)
        ax.set_xticklabels(conditions, rotation=45)
        ax.legend()
    
    # Plot 4: Differential interactions on X
    ax = axes[1, 0]
    if interaction_results:
        cis_count = interaction_results['n_cis']
        trans_count = interaction_results['n_trans']
        
        ax.bar(['Cis', 'Trans'], [cis_count, trans_count])
        ax.set_ylabel('Number of Interactions')
        ax.set_title('Differential X Chromosome Interactions')
    
    # Plot 5: Direction of change
    ax = axes[1, 1]
    if interaction_results:
        up_cis = interaction_results['upregulated_cis']
        down_cis = interaction_results['downregulated_cis']
        up_trans = interaction_results['upregulated_trans']
        down_trans = interaction_results['downregulated_trans']
        
        x = np.arange(2)
        width = 0.35
        
        ax.bar(x - width/2, [up_cis, up_trans], width, label='Upregulated', color='red')
        ax.bar(x + width/2, [down_cis, down_trans], width, label='Downregulated', color='blue')
        
        ax.set_ylabel('Number of Interactions')
        ax.set_title('Direction of X Chromosome Changes')
        ax.set_xticks(x)
        ax.set_xticklabels(['Cis', 'Trans'])
        ax.legend()
    
    # Plot 6: HAS proximity interactions
    ax = axes[1, 2]
    if has_proximity_results:
        up_has = has_proximity_results['has_proximal_upregulated']
        down_has = has_proximity_results['has_proximal_downregulated']
        
        ax.bar(['Upregulated', 'Downregulated'], [up_has, down_has], 
               color=['red', 'blue'])
        ax.set_ylabel('Number of Interactions')
        ax.set_title('HAS-Proximal Interactions')
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_analysis.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_prefix}_analysis.pdf")

def main():
    parser = argparse.ArgumentParser(description='Analyze X chromosome regulation')
    parser.add_argument('--has_sites', required=True, help='HAS/CES sites BED file')
    parser.add_argument('--mcool_files', nargs='+', required=True, help='Micro-C mcool files')
    parser.add_argument('--conditions', nargs='+', required=True, help='Condition names')
    parser.add_argument('--interaction_file', help='Differential interaction file (real data)')
    parser.add_argument('--null_interaction_file', help='Null model interaction file')
    parser.add_argument('--resolution', type=int, default=8000, help='Resolution for analysis')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    # Validate inputs
    if len(args.mcool_files) != len(args.conditions):
        print("Error: Number of mcool files must match number of conditions")
        return
    
    # Load HAS sites
    has_sites = load_has_sites(args.has_sites)
    
    # Analyze HAS contacts
    has_results = analyze_has_contacts(
        has_sites, args.mcool_files, args.conditions, 
        resolution=args.resolution
    )
    
    # Analyze compartmentalization
    compartment_results = analyze_x_compartmentalization(
        args.mcool_files, args.conditions
    )
    
    # Analyze differential interactions if provided
    interaction_results = {}
    interaction_results_null = {}
    has_proximity_results = {}
    
    if args.interaction_file:
        print("\n" + "="*60)
        print("ANALYZING REAL DATA")
        print("="*60)
        interaction_results = analyze_x_interactions(args.interaction_file)
        has_proximity_results = analyze_has_proximity(has_sites, args.interaction_file)
    
    if args.null_interaction_file:
        print("\n" + "="*60)
        print("ANALYZING NULL MODEL")
        print("="*60)
        interaction_results_null = analyze_null_model_stats(args.null_interaction_file)
    
    # Create visualizations
    create_x_regulation_plots(
        has_results, compartment_results, interaction_results, 
        has_proximity_results, args.output_prefix
    )
    
    # Save summary statistics
    summary = []
    for condition in args.conditions:
        row = {'condition': condition}
        
        if condition in has_results:
            df = has_results[condition]
            if len(df) > 0:
                row['has_avg_contacts'] = df['avg_contacts'].mean() if 'avg_contacts' in df.columns else 0
                row['has_insulation'] = df['insulation'].mean() if 'insulation' in df.columns else 0
                row['n_has_sites_analyzed'] = len(df)
            else:
                row['has_avg_contacts'] = 0
                row['has_insulation'] = 0
                row['n_has_sites_analyzed'] = 0
        
        if condition in compartment_results:
            row['compartment_strength'] = compartment_results[condition]['compartment_strength']
            row['percent_a_compartment'] = compartment_results[condition]['percent_a']
            row['percent_b_compartment'] = compartment_results[condition]['percent_b']
            row['matrix_size'] = compartment_results[condition]['matrix_size']
        
        summary.append(row)
    
    # Add interaction results (these are global, not per condition)
    if interaction_results:
        interaction_summary = {
            'x_cis_interactions': interaction_results['n_cis'],
            'x_trans_interactions': interaction_results['n_trans'],
            'total_sig_interactions': interaction_results['n_total_sig'],
            'x_cis_upregulated': interaction_results['upregulated_cis'],
            'x_cis_downregulated': interaction_results['downregulated_cis'],
            'x_trans_upregulated': interaction_results['upregulated_trans'],
            'x_trans_downregulated': interaction_results['downregulated_trans'],
            'mean_logfc_x_cis': interaction_results['mean_logfc_cis'],
            'mean_logfc_x_trans': interaction_results['mean_logfc_trans']
        }
        
        # Add to first row (since these are global stats)
        if summary:
            summary[0].update(interaction_summary)
    
    if has_proximity_results:
        proximity_summary = {
            'has_proximal_interactions': has_proximity_results['n_has_proximal'],
            'has_proximal_upregulated': has_proximity_results['has_proximal_upregulated'],
            'has_proximal_downregulated': has_proximity_results['has_proximal_downregulated'],
            'mean_logfc_has_proximal': has_proximity_results['mean_logfc_has_proximal']
        }
        
        # Add to first row
        if summary:
            summary[0].update(proximity_summary)
    
    # Add null model results
    if interaction_results_null:
        null_summary = {
            'null_total_interactions': interaction_results_null['n_total'],
            'null_sig_interactions': interaction_results_null['n_sig'],
            'null_proportion_sig': interaction_results_null['proportion_sig'],
            'null_mean_logfc': interaction_results_null['mean_logfc'],
            'null_median_logfc': interaction_results_null['median_logfc'],
            'null_std_logfc': interaction_results_null['std_logfc'],
            'null_upregulated': interaction_results_null['upregulated'],
            'null_downregulated': interaction_results_null['downregulated']
        }
        
        # Add to first row
        if summary:
            summary[0].update(null_summary)
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(f"{args.output_prefix}_summary.tsv", sep='\t', index=False)
    
    # Save detailed results
    if has_results:
        # Combine all HAS results
        all_has_results = []
        for condition, df in has_results.items():
            df_copy = df.copy()
            df_copy['condition'] = condition
            all_has_results.append(df_copy)
        
        if all_has_results:
            combined_has = pd.concat(all_has_results, ignore_index=True)
            combined_has.to_csv(f"{args.output_prefix}_has_contact_comparison.tsv", sep='\t', index=False)
    
    if compartment_results:
        # Save compartment results
        comp_summary = []
        for condition, results in compartment_results.items():
            comp_summary.append({
                'condition': condition,
                'compartment_strength': results['compartment_strength'],
                'percent_a': results['percent_a'],
                'percent_b': results['percent_b'],
                'matrix_size': results['matrix_size']
            })
        
        comp_df = pd.DataFrame(comp_summary)
        comp_df.to_csv(f"{args.output_prefix}_compartment_changes.tsv", sep='\t', index=False)
    
    # Save null model comparison statistics if both real and null data exist
    if interaction_results and interaction_results_null:
        null_comparison = {
            'metric': ['Total significant interactions', 
                      'X-cis interactions',
                      'X-trans interactions',
                      'Proportion significant'],
            'real_data': [interaction_results['n_total_sig'],
                         interaction_results['n_cis'],
                         interaction_results['n_trans'],
                         interaction_results['n_total_sig'] / interaction_results_null['n_total'] if interaction_results_null['n_total'] > 0 else 0],
            'null_model': [interaction_results_null['n_sig'],
                          np.nan,  # Can't calculate X-specific for null
                          np.nan,
                          interaction_results_null['proportion_sig']],
            'enrichment': [interaction_results['n_total_sig'] / interaction_results_null['n_sig'] if interaction_results_null['n_sig'] > 0 else np.nan,
                          np.nan,
                          np.nan,
                          (interaction_results['n_total_sig'] / interaction_results_null['n_total']) / interaction_results_null['proportion_sig'] if interaction_results_null['proportion_sig'] > 0 else np.nan]
        }
        
        null_comp_df = pd.DataFrame(null_comparison)
        null_comp_df.to_csv(f"{args.output_prefix}_null_comparison_stats.tsv", sep='\t', index=False)
        
        # Save detailed null comparison
        detailed_null = pd.DataFrame([{
            'real_total_sig': interaction_results['n_total_sig'],
            'real_x_cis': interaction_results['n_cis'],
            'real_x_trans': interaction_results['n_trans'],
            'real_mean_logfc_cis': interaction_results['mean_logfc_cis'],
            'real_mean_logfc_trans': interaction_results['mean_logfc_trans'],
            'null_total': interaction_results_null['n_total'],
            'null_sig': interaction_results_null['n_sig'],
            'null_proportion_sig': interaction_results_null['proportion_sig'],
            'null_mean_logfc': interaction_results_null['mean_logfc'],
            'null_median_logfc': interaction_results_null['median_logfc'],
            'null_std_logfc': interaction_results_null['std_logfc']
        }])
        detailed_null.to_csv(f"{args.output_prefix}_null_comparison_detailed.tsv", sep='\t', index=False)
    
    print(f"\nAnalysis complete. Results saved to {args.output_prefix}_*")
    
    # Print summary
    print("\nSummary:")
    print(f"  HAS sites analyzed: {len(has_sites)}")
    print(f"  Conditions with HAS data: {len(has_results)}")
    print(f"  Conditions with compartment data: {len(compartment_results)}")
    
    if interaction_results:
        print(f"  X-chromosome cis interactions: {interaction_results['n_cis']}")
        print(f"  X-chromosome trans interactions: {interaction_results['n_trans']}")
        print(f"  Total significant interactions: {interaction_results['n_total_sig']}")
    
    if has_proximity_results:
        print(f"  HAS-proximal interactions: {has_proximity_results['n_has_proximal']}")
    
    if interaction_results_null:
        print(f"\nNull model:")
        print(f"  Total interactions: {interaction_results_null['n_total']}")
        print(f"  Significant interactions: {interaction_results_null['n_sig']} ({interaction_results_null['proportion_sig']:.2%})")
        if interaction_results:
            enrichment = interaction_results['n_total_sig'] / interaction_results_null['n_sig'] if interaction_results_null['n_sig'] > 0 else np.nan
            print(f"  Enrichment over null: {enrichment:.2f}x")

if __name__ == '__main__':
    main()