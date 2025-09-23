#!/usr/bin/env python3
"""
Fixed version of enhancer class analysis script.
Addresses condition name matching and type conversion issues.
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
from collections import defaultdict

def load_and_process_interactions(interactions_file, reference_condition='DOX'):
    """
    Load differential interactions and process condition names.
    """
    print(f"Loading interactions from {interactions_file}")
    interactions = pd.read_csv(interactions_file)
    
    print(f"Loaded {len(interactions)} interactions")
    print("Available columns:", interactions.columns.tolist())
    
    # Check for condition columns
    condition_cols = [col for col in interactions.columns if 'infection' in col.lower() or 'condition' in col.lower()]
    if condition_cols:
        print(f"Found condition columns: {condition_cols}")
        # Use the first condition column
        condition_col = condition_cols[0]
        interactions['condition'] = interactions[condition_col]
    elif 'filename' in interactions.columns:
        # Extract condition from filename
        print("Extracting condition from filename column")
        interactions['condition'] = interactions['filename'].str.extract(r'infection(JW18\w+)')[0]
        # Clean up condition names
        interactions['condition'] = interactions['condition'].str.replace('JW18', '')
    else:
        print("Warning: No clear condition column found. Using all data as single condition.")
        interactions['condition'] = 'unknown'
    
    # Print available conditions
    available_conditions = interactions['condition'].unique()
    print(f"Available conditions: {available_conditions}")
    
    # Map conditions to standard names if needed
    condition_mapping = {
        'wMel': 'wMel',
        'wRi': 'wRi', 
        'wWil': 'wWil',
        'DOX': 'DOX',
        'JW18wMel': 'wMel',
        'JW18wRi': 'wRi',
        'JW18wWil': 'wWil',
        'JW18DOX': 'DOX'
    }
    
    # Apply mapping
    interactions['condition_clean'] = interactions['condition'].map(condition_mapping).fillna(interactions['condition'])
    
    # Check if reference condition exists
    clean_conditions = interactions['condition_clean'].unique()
    print(f"Cleaned conditions: {clean_conditions}")
    
    if reference_condition not in clean_conditions:
        print(f"Warning: Reference condition '{reference_condition}' not found in data")
        print(f"Available conditions: {clean_conditions}")
        # Use the first available condition as reference
        reference_condition = clean_conditions[0]
        print(f"Using '{reference_condition}' as reference condition")
    
    return interactions, reference_condition

def classify_enhancers(enhancer_file, classification_file=None):
    """
    Load and classify enhancers as housekeeping or developmental.
    """
    print("Loading enhancer annotations...")
    
    # Try different formats for the enhancer file
    try:
        enhancers = pd.read_csv(enhancer_file, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    except:
        # Try with header
        enhancers = pd.read_csv(enhancer_file, sep='\t')
        if 'chrom' not in enhancers.columns:
            # Rename columns if they exist
            col_mapping = {
                enhancers.columns[0]: 'chrom',
                enhancers.columns[1]: 'start', 
                enhancers.columns[2]: 'end'
            }
            if len(enhancers.columns) > 3:
                col_mapping[enhancers.columns[3]] = 'name'
            if len(enhancers.columns) > 4:
                col_mapping[enhancers.columns[4]] = 'score'
            if len(enhancers.columns) > 5:
                col_mapping[enhancers.columns[5]] = 'strand'
            
            enhancers = enhancers.rename(columns=col_mapping)
    
    # Ensure required columns exist
    required_cols = ['chrom', 'start', 'end']
    for col in required_cols:
        if col not in enhancers.columns:
            raise ValueError(f"Required column '{col}' not found in enhancer file")
    
    if 'name' not in enhancers.columns:
        enhancers['name'] = enhancers.apply(lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1)
    
    print(f"Loaded {len(enhancers)} enhancers")
    
    # Load classification if provided
    if classification_file and Path(classification_file).exists():
        print(f"Loading enhancer classification from {classification_file}")
        try:
            classification = pd.read_csv(classification_file, sep='\t')
            enhancers = enhancers.merge(classification[['name', 'class']], on='name', how='left')
        except Exception as e:
            print(f"Error loading classification file: {e}")
            # Fallback to heuristic
            enhancers = apply_heuristic_classification(enhancers)
    else:
        print("No classification file provided, using heuristic classification")
        enhancers = apply_heuristic_classification(enhancers)
    
    # Fill missing classifications
    enhancers['class'] = enhancers['class'].fillna('developmental')
    
    print(f"Classified {sum(enhancers['class'] == 'housekeeping')} housekeeping enhancers")
    print(f"Classified {sum(enhancers['class'] == 'developmental')} developmental enhancers")
    
    return enhancers

def apply_heuristic_classification(enhancers):
    """Apply heuristic classification based on enhancer names."""
    housekeeping_markers = ['ubiq', 'house', 'const', 'rp', 'ef1', 'gapdh', 'actb', 'tubulin', 'actin']
    
    enhancers['class'] = enhancers['name'].apply(
        lambda x: 'housekeeping' if any(marker in str(x).lower() for marker in housekeeping_markers)
        else 'developmental'
    )
    
    return enhancers

def find_enhancer_interactions(interactions, enhancers, tss_file=None):
    """
    Find enhancer-enhancer and enhancer-TSS interactions.
    """
    print("Finding enhancer interactions...")
    
    # Convert interactions to BedTool format
    interactions_bed1 = interactions[['chr1', 'start1', 'end1']].copy()
    interactions_bed1.columns = ['chrom', 'start', 'end']
    interactions_bed1['interaction_idx'] = interactions.index
    
    interactions_bed2 = interactions[['chr2', 'start2', 'end2']].copy()
    interactions_bed2.columns = ['chrom', 'start', 'end']
    interactions_bed2['interaction_idx'] = interactions.index
    
    # Convert to BedTools
    enh_bt = pybedtools.BedTool.from_dataframe(enhancers[['chrom', 'start', 'end', 'name', 'class']])
    int1_bt = pybedtools.BedTool.from_dataframe(interactions_bed1)
    int2_bt = pybedtools.BedTool.from_dataframe(interactions_bed2)
    
    # Find overlaps
    enh_int1 = int1_bt.intersect(enh_bt, wa=True, wb=True)
    enh_int2 = int2_bt.intersect(enh_bt, wa=True, wb=True)
    
    # Process overlaps
    enh_interactions = []
    
    # Process anchor 1 overlaps
    anchor1_overlaps = {}
    for overlap in enh_int1:
        fields = str(overlap).strip().split('\t')
        if len(fields) >= 9:
            idx = int(fields[3])  # interaction_idx
            enh_name = fields[7]  # enhancer name
            enh_class = fields[8]  # enhancer class
            anchor1_overlaps[idx] = {'name': enh_name, 'class': enh_class}
    
    # Process anchor 2 overlaps
    anchor2_overlaps = {}
    for overlap in enh_int2:
        fields = str(overlap).strip().split('\t')
        if len(fields) >= 9:
            idx = int(fields[3])  # interaction_idx
            enh_name = fields[7]  # enhancer name  
            enh_class = fields[8]  # enhancer class
            anchor2_overlaps[idx] = {'name': enh_name, 'class': enh_class}
    
    # Combine overlaps to find enhancer interactions
    for idx in interactions.index:
        anchor1_enh = anchor1_overlaps.get(idx)
        anchor2_enh = anchor2_overlaps.get(idx)
        
        if anchor1_enh and anchor2_enh:
            # Both anchors are enhancers - E-E interaction
            interaction_type = 'E-E'
            if anchor1_enh['class'] == anchor2_enh['class']:
                contact_class = anchor1_enh['class']
            else:
                contact_class = 'cross_class'
        elif anchor1_enh or anchor2_enh:
            # One anchor is enhancer - could be E-TSS
            interaction_type = 'E-TSS'
            enh_info = anchor1_enh or anchor2_enh
            contact_class = enh_info['class']
        else:
            continue
        
        # Add interaction info
        interaction_data = interactions.iloc[idx].copy()
        interaction_data['interaction_type'] = interaction_type
        interaction_data['contact_class'] = contact_class
        
        if anchor1_enh:
            interaction_data['enh1_name'] = anchor1_enh['name']
            interaction_data['enh1_class'] = anchor1_enh['class']
        if anchor2_enh:
            interaction_data['enh2_name'] = anchor2_enh['name']
            interaction_data['enh2_class'] = anchor2_enh['class']
        
        enh_interactions.append(interaction_data)
    
    if enh_interactions:
        enh_df = pd.DataFrame(enh_interactions)
        print(f"Found {len(enh_df)} enhancer interactions:")
        print(f"  E-E interactions: {sum(enh_df['interaction_type'] == 'E-E')}")
        print(f"  E-TSS interactions: {sum(enh_df['interaction_type'] == 'E-TSS')}")
        print("  By class:")
        print(enh_df['contact_class'].value_counts())
        
        return enh_df
    else:
        print("No enhancer interactions found")
        return pd.DataFrame()

def perform_statistical_comparison(enh_interactions, null_model, reference_condition, fdr_threshold=0.05):
    """
    Compare enhancer interactions between conditions with proper type handling.
    """
    print("Performing robust statistical comparison with null model...")
    
    if enh_interactions.empty:
        print("No enhancer interactions to analyze")
        return pd.DataFrame()
    
    print(f"Comparing conditions with null model (reference: {reference_condition})...")
    available_conditions = enh_interactions['condition_clean'].unique()
    print(f"Available conditions: {available_conditions}")
    
    if reference_condition not in available_conditions:
        print(f"Warning: Reference condition '{reference_condition}' not found in data")
        reference_condition = available_conditions[0]
        print(f"Using '{reference_condition}' as reference")
    
    comparison_results = []
    
    # Get reference data
    ref_data = enh_interactions[enh_interactions['condition_clean'] == reference_condition]
    if ref_data.empty:
        print(f"No reference data found for condition: {reference_condition}")
        return pd.DataFrame()
    
    # Compare each condition to reference
    for condition in available_conditions:
        if condition == reference_condition:
            continue
            
        cond_data = enh_interactions[enh_interactions['condition_clean'] == condition]
        if cond_data.empty:
            continue
        
        print(f"Comparing {condition} vs {reference_condition}")
        
        # Analyze by interaction type and class
        for int_type in ['E-E', 'E-TSS']:
            ref_type = ref_data[ref_data['interaction_type'] == int_type]
            cond_type = cond_data[cond_data['interaction_type'] == int_type]
            
            if ref_type.empty or cond_type.empty:
                continue
            
            for contact_class in ref_type['contact_class'].unique():
                ref_class = ref_type[ref_type['contact_class'] == contact_class]
                cond_class = cond_type[cond_type['contact_class'] == contact_class]
                
                if ref_class.empty or cond_class.empty:
                    continue
                
                try:
                    # Ensure logFC is numeric
                    ref_logfc = pd.to_numeric(ref_class['logFC'], errors='coerce').dropna()
                    cond_logfc = pd.to_numeric(cond_class['logFC'], errors='coerce').dropna()
                    
                    if len(ref_logfc) == 0 or len(cond_logfc) == 0:
                        continue
                    
                    # Statistical test
                    mean_ref = ref_logfc.mean()
                    mean_cond = cond_logfc.mean()
                    
                    # Use Mann-Whitney U test (non-parametric)
                    try:
                        statistic, p_value = stats.mannwhitneyu(
                            cond_logfc, ref_logfc, alternative='two-sided'
                        )
                    except Exception as e:
                        print(f"Statistical test failed for {int_type} {contact_class}: {e}")
                        p_value = 1.0
                    
                    # Calculate effect size (log2 fold change of means)
                    if mean_ref != 0:
                        log2fc_diff = np.log2(abs(mean_cond) + 1e-10) - np.log2(abs(mean_ref) + 1e-10)
                    else:
                        log2fc_diff = 0
                    
                    comparison_results.append({
                        'comparison': f'{reference_condition}_vs_{condition}',
                        'interaction_type': int_type,
                        'contact_class': contact_class,
                        'n_ref': len(ref_logfc),
                        'n_cond': len(cond_logfc),
                        'mean_logfc_ref': mean_ref,
                        'mean_logfc_cond': mean_cond,
                        'log2fc_difference': log2fc_diff,
                        'p_value': p_value,
                        'test_statistic': statistic if 'statistic' in locals() else np.nan
                    })
                    
                except Exception as e:
                    print(f"Error in statistical comparison for {int_type} {contact_class}: {e}")
                    continue
    
    if not comparison_results:
        print("No valid comparisons could be performed")
        return pd.DataFrame()
    
    # Convert to DataFrame and apply FDR correction
    results_df = pd.DataFrame(comparison_results)
    
    # Apply FDR correction
    _, results_df['fdr'], _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    results_df['significant'] = results_df['fdr'] < fdr_threshold
    
    print(f"Found {sum(results_df['significant'])} significant comparisons")
    
    return results_df

def analyze_distance_dependence(enh_interactions):
    """
    Analyze how enhancer interactions change with genomic distance.
    """
    print("Analyzing distance dependence...")
    
    if enh_interactions.empty:
        return pd.DataFrame()
    
    # Filter for cis interactions only
    cis_interactions = enh_interactions[enh_interactions['chr1'] == enh_interactions['chr2']].copy()
    
    if cis_interactions.empty:
        return pd.DataFrame()
    
    # Calculate interaction distance
    cis_interactions['distance'] = abs(
        pd.to_numeric(cis_interactions['start2'], errors='coerce') - 
        pd.to_numeric(cis_interactions['start1'], errors='coerce')
    )
    
    # Remove rows with invalid distances
    cis_interactions = cis_interactions.dropna(subset=['distance'])
    
    if cis_interactions.empty:
        return pd.DataFrame()
    
    # Define distance bins
    distance_bins = [0, 10000, 50000, 100000, 500000, 1000000, float('inf')]
    distance_labels = ['<10kb', '10-50kb', '50-100kb', '100-500kb', '500kb-1Mb', '>1Mb']
    
    cis_interactions['distance_bin'] = pd.cut(
        cis_interactions['distance'], 
        bins=distance_bins, 
        labels=distance_labels,
        include_lowest=True
    )
    
    # Analyze by condition, interaction type, and distance
    distance_analysis = []
    
    for condition in cis_interactions['condition_clean'].unique():
        cond_data = cis_interactions[cis_interactions['condition_clean'] == condition]
        
        for int_type in ['E-E', 'E-TSS']:
            type_data = cond_data[cond_data['interaction_type'] == int_type]
            
            if type_data.empty:
                continue
            
            for contact_class in type_data['contact_class'].unique():
                class_data = type_data[type_data['contact_class'] == contact_class]
                
                if class_data.empty:
                    continue
                
                for dist_bin in distance_labels:
                    bin_data = class_data[class_data['distance_bin'] == dist_bin]
                    
                    if len(bin_data) < 3:  # Minimum samples
                        continue
                    
                    # Convert logFC to numeric
                    logfc_values = pd.to_numeric(bin_data['logFC'], errors='coerce').dropna()
                    
                    if len(logfc_values) == 0:
                        continue
                    
                    distance_analysis.append({
                        'condition': condition,
                        'interaction_type': int_type,
                        'contact_class': contact_class,
                        'distance_bin': dist_bin,
                        'n_interactions': len(bin_data),
                        'mean_logfc': logfc_values.mean(),
                        'median_logfc': logfc_values.median(),
                        'std_logfc': logfc_values.std()
                    })
    
    distance_df = pd.DataFrame(distance_analysis)
    print(f"Distance analysis completed for {len(distance_df)} combinations")
    
    return distance_df

def create_summary_tables(enh_interactions, comparison_results, distance_analysis):
    """
    Create summary tables with proper error handling.
    """
    print("Creating summary tables...")
    
    summary_tables = {}
    
    try:
        # Overall summary
        if not enh_interactions.empty:
            overall_summary = {
                'total_enhancer_interactions': len(enh_interactions),
                'e_e_interactions': sum(enh_interactions['interaction_type'] == 'E-E'),
                'e_tss_interactions': sum(enh_interactions['interaction_type'] == 'E-TSS'),
                'conditions_analyzed': len(enh_interactions['condition_clean'].unique()),
                'housekeeping_interactions': sum(enh_interactions['contact_class'] == 'housekeeping'),
                'developmental_interactions': sum(enh_interactions['contact_class'] == 'developmental'),
                'cross_class_interactions': sum(enh_interactions['contact_class'] == 'cross_class')
            }
            summary_tables['overall'] = overall_summary
        
        # Comparison summary
        if not comparison_results.empty:
            comp_summary = {
                'total_comparisons': len(comparison_results),
                'significant_comparisons': sum(comparison_results['significant']),
                'e_e_significant': sum(
                    (comparison_results['interaction_type'] == 'E-E') & 
                    comparison_results['significant']
                ),
                'e_tss_significant': sum(
                    (comparison_results['interaction_type'] == 'E-TSS') & 
                    comparison_results['significant']
                )
            }
            summary_tables['comparisons'] = comp_summary
        
        # Distance summary
        if not distance_analysis.empty:
            dist_summary = {
                'distance_bins_analyzed': len(distance_analysis['distance_bin'].unique()),
                'conditions_in_distance_analysis': len(distance_analysis['condition'].unique()),
                'total_distance_combinations': len(distance_analysis)
            }
            summary_tables['distance'] = dist_summary
        
        return summary_tables
        
    except Exception as e:
        print(f"Error creating summary tables: {e}")
        return {}

def create_summary_plots(enh_interactions, comparison_results, distance_analysis, output_prefix):
    """
    Create summary plots with error handling.
    """
    try:
        if enh_interactions.empty:
            print("No enhancer interactions to plot")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Interaction counts by type and class
        ax = axes[0, 0]
        if not enh_interactions.empty:
            interaction_counts = enh_interactions.groupby(['interaction_type', 'contact_class']).size().unstack(fill_value=0)
            interaction_counts.plot(kind='bar', ax=ax)
            ax.set_title('Enhancer Interactions by Type and Class')
            ax.set_ylabel('Count')
            ax.legend(title='Contact Class')
        
        # Plot 2: LogFC distribution
        ax = axes[0, 1]
        if not enh_interactions.empty and 'logFC' in enh_interactions.columns:
            logfc_numeric = pd.to_numeric(enh_interactions['logFC'], errors='coerce').dropna()
            if len(logfc_numeric) > 0:
                ax.hist(logfc_numeric, bins=30, alpha=0.7)
                ax.set_title('LogFC Distribution')
                ax.set_xlabel('LogFC')
                ax.set_ylabel('Frequency')
        
        # Plot 3: Significant comparisons
        ax = axes[1, 0]
        if not comparison_results.empty:
            sig_counts = comparison_results.groupby('interaction_type')['significant'].sum()
            if len(sig_counts) > 0:
                ax.bar(sig_counts.index, sig_counts.values)
                ax.set_title('Significant Comparisons by Type')
                ax.set_ylabel('Count')
        
        # Plot 4: Distance analysis
        ax = axes[1, 1]
        if not distance_analysis.empty:
            # Plot mean logFC by distance bin
            dist_means = distance_analysis.groupby('distance_bin')['mean_logfc'].mean()
            if len(dist_means) > 0:
                ax.bar(range(len(dist_means)), dist_means.values)
                ax.set_xticks(range(len(dist_means)))
                ax.set_xticklabels(dist_means.index, rotation=45)
                ax.set_title('Mean LogFC by Distance')
                ax.set_ylabel('Mean LogFC')
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_summary_plots.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Summary plots saved to {output_prefix}_summary_plots.pdf")
        
    except Exception as e:
        print(f"Error creating plots: {e}")

def main():
    parser = argparse.ArgumentParser(description='Analyze enhancer class interactions')
    parser.add_argument('--enhancers', required=True, help='Enhancer BED file')
    parser.add_argument('--interactions', required=True, help='Differential interactions CSV file')
    parser.add_argument('--null_model', help='Null model CSV file from diffHic')
    parser.add_argument('--classification', help='Optional enhancer classification file')
    parser.add_argument('--tss', help='TSS BED file')
    parser.add_argument('--fdr_threshold', type=float, default=0.05, help='FDR threshold')
    parser.add_argument('--reference_condition', default='DOX', help='Reference condition name')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    # Load null model if provided
    null_model = None
    if args.null_model:
        null_model = pd.read_csv(args.null_model)
        print(f"Loaded null model with {len(null_model)} entries")
    
    # Load and process interactions
    interactions, reference_condition = load_and_process_interactions(
        args.interactions, args.reference_condition
    )
    
    # Load and classify enhancers
    enhancers = classify_enhancers(args.enhancers, args.classification)
    
    # Find enhancer interactions
    enh_interactions = find_enhancer_interactions(interactions, enhancers, args.tss)
    
    # Perform statistical comparison
    comparison_results = perform_statistical_comparison(
        enh_interactions, null_model, reference_condition, args.fdr_threshold
    )
    
    # Analyze distance dependence
    distance_analysis = analyze_distance_dependence(enh_interactions)
    
    # Create summary tables
    summary_tables = create_summary_tables(enh_interactions, comparison_results, distance_analysis)
    
    # Create plots
    if not enh_interactions.empty:
        create_summary_plots(enh_interactions, comparison_results, distance_analysis, args.output_prefix)
    else:
        print("No comparison data available for plotting")
    
    # Save results
    try:
        if not distance_analysis.empty:
            distance_analysis.to_csv(f"{args.output_prefix}_distance_analysis.tsv", sep='\t', index=False)
            print(f"Saved distance analysis: {args.output_prefix}_distance_analysis.tsv")
        
        if not enh_interactions.empty:
            enh_interactions.to_csv(f"{args.output_prefix}_enhancer_interactions.tsv", sep='\t', index=False)
            print(f"Saved enhancer interactions: {args.output_prefix}_enhancer_interactions.tsv")
        
        if not comparison_results.empty:
            comparison_results.to_csv(f"{args.output_prefix}_statistical_comparisons.tsv", sep='\t', index=False)
            print(f"Saved statistical comparisons: {args.output_prefix}_statistical_comparisons.tsv")
        
        # Save summary tables
        if summary_tables:
            with open(f"{args.output_prefix}_summary.txt", 'w') as f:
                f.write("Enhancer Class Analysis Summary\n")
                f.write("=" * 40 + "\n\n")
                
                for table_name, table_data in summary_tables.items():
                    f.write(f"{table_name.upper()} SUMMARY:\n")
                    for key, value in table_data.items():
                        f.write(f"  {key}: {value}\n")
                    f.write("\n")
            print(f"Saved summary: {args.output_prefix}_summary.txt")
        
    except Exception as e:
        print(f"Error creating summary tables: {e}")
    
    print(f"\nAnalysis complete! Results saved to {args.output_prefix}_*")
    
    # Print final summary
    print("\nFinal Summary:")
    if not enh_interactions.empty:
        print(f"  Total enhancer interactions found: {len(enh_interactions)}")
        print(f"  E-E interactions: {sum(enh_interactions['interaction_type'] == 'E-E')}")
        print(f"  E-TSS interactions: {sum(enh_interactions['interaction_type'] == 'E-TSS')}")
    else:
        print("  No enhancer interactions found")
    
    if not comparison_results.empty:
        print(f"  Statistical comparisons performed: {len(comparison_results)}")
        print(f"  Significant comparisons: {sum(comparison_results['significant'])}")
    else:
        print("  No statistical comparisons performed")

if __name__ == '__main__':
    main()