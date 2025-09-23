#!/usr/bin/env python3
"""
Analyze enhancer classes and their chromatin interactions across conditions.
Uses pre-identified differential interactions from diffHic analysis.
Based on Zabidi et al. 2015 classification of housekeeping vs developmental enhancers.
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
import warnings
warnings.filterwarnings('ignore')

def classify_enhancers(enhancer_file, classification_file=None):
    """
    Classify enhancers as housekeeping or developmental.
    If classification_file provided, use it; otherwise use heuristics.
    """
    print("Loading enhancer annotations...")
    
    # Try to read the enhancer file - handle different formats
    try:
        enhancers = pd.read_csv(enhancer_file, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    except:
        # Try with header
        enhancers = pd.read_csv(enhancer_file, sep='\t')
        # Rename columns if they have different names
        column_mapping = {
            'chr': 'chrom', 'chromosome': 'chrom',
            'begin': 'start', 'pos': 'start',
            'stop': 'end', 'finish': 'end',
            'id': 'name', 'element_name': 'name'
        }
        enhancers.rename(columns=column_mapping, inplace=True)
    
    # Ensure we have required columns
    required_cols = ['chrom', 'start', 'end']
    missing_cols = [col for col in required_cols if col not in enhancers.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Add name column if missing
    if 'name' not in enhancers.columns:
        enhancers['name'] = enhancers.apply(
            lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1
        )
    
    # Add score and strand if missing
    if 'score' not in enhancers.columns:
        enhancers['score'] = 0
    if 'strand' not in enhancers.columns:
        enhancers['strand'] = '.'
    
    if classification_file and Path(classification_file).exists():
        # Load actual classification from file
        print(f"Loading enhancer classification from {classification_file}")
        classification = pd.read_csv(classification_file, sep='\t')
        enhancers = enhancers.merge(classification[['name', 'class']], on='name', how='left')
        enhancers['class'].fillna('unknown', inplace=True)
    else:
        # Use heuristic based on enhancer names or nearby genes
        print("Using heuristic classification (consider using actual Zabidi et al. data)")
        # Housekeeping enhancers are typically associated with broadly expressed genes
        housekeeping_markers = ['ubiq', 'house', 'const', 'rp', 'ef1', 'gapdh', 'actb', 'rpl', 'rps']
        
        enhancers['class'] = enhancers['name'].apply(
            lambda x: 'housekeeping' if any(marker in str(x).lower() for marker in housekeeping_markers)
            else 'developmental'
        )
    
    print(f"Classified {sum(enhancers['class'] == 'housekeeping')} housekeeping enhancers")
    print(f"Classified {sum(enhancers['class'] == 'developmental')} developmental enhancers")
    print(f"Unknown classification: {sum(enhancers['class'] == 'unknown')}")
    
    return enhancers

def load_differential_interactions(interactions_file, fdr_threshold=0.05):
    """
    Load pre-identified differential interactions from diffHic analysis.
    """
    print(f"Loading differential interactions from {interactions_file}")
    
    # Load the data
    interactions = pd.read_csv(interactions_file)
    
    # Filter by significance if FDR column exists
    if 'FDR' in interactions.columns:
        significant_interactions = interactions[interactions['FDR'] < fdr_threshold].copy()
        print(f"Found {len(significant_interactions)} significant interactions (FDR < {fdr_threshold})")
        print(f"Total interactions in file: {len(interactions)}")
    else:
        print("No FDR column found - using all interactions")
        significant_interactions = interactions.copy()
    
    # Check required columns
    required_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
    missing_cols = [col for col in required_cols if col not in significant_interactions.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in interactions file: {missing_cols}")
    
    # Add logFC column if missing (set to 0)
    if 'logFC' not in significant_interactions.columns:
        significant_interactions['logFC'] = 0
        print("Warning: No logFC column found - setting to 0")
    
    # Add interaction_distance if missing
    if 'interaction_distance' not in significant_interactions.columns:
        significant_interactions['interaction_distance'] = np.where(
            significant_interactions['chr1'] == significant_interactions['chr2'],
            np.abs(significant_interactions['start2'] - significant_interactions['start1']),
            -1  # Trans interactions
        )
    
    # Add interaction_type if missing
    if 'interaction_type' not in significant_interactions.columns:
        significant_interactions['interaction_type'] = np.where(
            significant_interactions['chr1'] == significant_interactions['chr2'],
            'cis', 'trans'
        )
    
    # Add infection/condition column if missing
    if 'infection' not in significant_interactions.columns:
        # Try to infer from filename
        if 'filename' in significant_interactions.columns:
            significant_interactions['infection'] = significant_interactions['filename'].str.extract(
                r'(wMel|wRi|wWil|DOX)'
            )[0]
        else:
            significant_interactions['infection'] = 'unknown'
    
    print(f"Interaction types: {significant_interactions['interaction_type'].value_counts().to_dict()}")
    if 'infection' in significant_interactions.columns:
        print(f"Conditions: {significant_interactions['infection'].value_counts().to_dict()}")
    
    return significant_interactions

def find_enhancer_interactions(enhancers, interactions, interaction_type='both'):
    """
    Find interactions involving enhancers.
    interaction_type: 'E-TSS', 'E-E', or 'both'
    """
    print(f"Finding enhancer interactions ({interaction_type})...")
    
    # Create BedTool objects
    enhancers_bed = pybedtools.BedTool.from_dataframe(
        enhancers[['chrom', 'start', 'end', 'name', 'score', 'strand', 'class']]
    )
    
    enhancer_interactions = []
    
    for _, interaction in interactions.iterrows():
        # Create BedTool objects for interaction anchors
        anchor1_bed = pybedtools.BedTool(
            f"{interaction['chr1']}\t{interaction['start1']}\t{interaction['end1']}", 
            from_string=True
        )
        anchor2_bed = pybedtools.BedTool(
            f"{interaction['chr2']}\t{interaction['start2']}\t{interaction['end2']}", 
            from_string=True
        )
        
        # Find overlaps with enhancers
        anchor1_overlaps = anchor1_bed.intersect(enhancers_bed, wa=True, wb=True)
        anchor2_overlaps = anchor2_bed.intersect(enhancers_bed, wa=True, wb=True)
        
        anchor1_enhancers = []
        anchor2_enhancers = []
        
        # Extract enhancer information for anchor 1
        for overlap in anchor1_overlaps:
            if len(overlap.fields) >= 10:  # Ensure we have enhancer info
                anchor1_enhancers.append({
                    'name': overlap.fields[6],
                    'class': overlap.fields[9]
                })
        
        # Extract enhancer information for anchor 2
        for overlap in anchor2_overlaps:
            if len(overlap.fields) >= 10:
                anchor2_enhancers.append({
                    'name': overlap.fields[6],
                    'class': overlap.fields[9]
                })
        
        # Classify interaction type
        if len(anchor1_enhancers) > 0 and len(anchor2_enhancers) > 0:
            # E-E interaction
            for enh1 in anchor1_enhancers:
                for enh2 in anchor2_enhancers:
                    if enh1['class'] == enh2['class']:
                        contact_class = enh1['class']
                    else:
                        contact_class = 'cross_class'
                    
                    enhancer_interactions.append({
                        'interaction_type': 'E-E',
                        'contact_class': contact_class,
                        'enhancer1': enh1['name'],
                        'enhancer2': enh2['name'],
                        'enhancer1_class': enh1['class'],
                        'enhancer2_class': enh2['class'],
                        'chr1': interaction['chr1'],
                        'start1': interaction['start1'],
                        'end1': interaction['end1'],
                        'chr2': interaction['chr2'],
                        'start2': interaction['start2'],
                        'end2': interaction['end2'],
                        'logFC': interaction['logFC'],
                        'FDR': interaction.get('FDR', 1.0),
                        'distance': interaction['interaction_distance'],
                        'infection': interaction.get('infection', 'unknown')
                    })
        
        elif len(anchor1_enhancers) > 0 or len(anchor2_enhancers) > 0:
            # Potential E-TSS interaction (one anchor is enhancer, other might be gene/TSS)
            enhancer_anchor = anchor1_enhancers if len(anchor1_enhancers) > 0 else anchor2_enhancers
            non_enhancer_anchor = 2 if len(anchor1_enhancers) > 0 else 1
            
            for enh in enhancer_anchor:
                enhancer_interactions.append({
                    'interaction_type': 'E-TSS',
                    'contact_class': enh['class'],
                    'enhancer': enh['name'],
                    'enhancer_class': enh['class'],
                    'chr1': interaction['chr1'],
                    'start1': interaction['start1'],
                    'end1': interaction['end1'],
                    'chr2': interaction['chr2'],
                    'start2': interaction['start2'],
                    'end2': interaction['end2'],
                    'logFC': interaction['logFC'],
                    'FDR': interaction.get('FDR', 1.0),
                    'distance': interaction['interaction_distance'],
                    'infection': interaction.get('infection', 'unknown')
                })
    
    enhancer_interactions_df = pd.DataFrame(enhancer_interactions)
    
    if len(enhancer_interactions_df) > 0:
        print(f"Found {len(enhancer_interactions_df)} enhancer interactions:")
        print(f"  E-E interactions: {sum(enhancer_interactions_df['interaction_type'] == 'E-E')}")
        print(f"  E-TSS interactions: {sum(enhancer_interactions_df['interaction_type'] == 'E-TSS')}")
        
        if 'contact_class' in enhancer_interactions_df.columns:
            print("  By class:")
            print(enhancer_interactions_df['contact_class'].value_counts().to_string())
    else:
        print("No enhancer interactions found!")
    
    return enhancer_interactions_df

def compare_conditions(enhancer_interactions, reference_condition='DOX'):
    """
    Compare enhancer interactions between conditions.
    """
    print(f"Comparing conditions (reference: {reference_condition})...")
    
    if 'infection' not in enhancer_interactions.columns:
        print("Warning: No infection/condition column found")
        return pd.DataFrame()
    
    conditions = enhancer_interactions['infection'].unique()
    print(f"Available conditions: {conditions}")
    
    if reference_condition not in conditions:
        print(f"Warning: Reference condition '{reference_condition}' not found in data")
        # Use the first available condition as reference
        reference_condition = conditions[0]
        print(f"Using '{reference_condition}' as reference")
    
    comparisons = []
    
    # Group by interaction type and contact class
    for interaction_type in enhancer_interactions['interaction_type'].unique():
        type_data = enhancer_interactions[enhancer_interactions['interaction_type'] == interaction_type]
        
        for contact_class in type_data['contact_class'].unique():
            class_data = type_data[type_data['contact_class'] == contact_class]
            
            # Get reference data
            ref_data = class_data[class_data['infection'] == reference_condition]
            
            for condition in conditions:
                if condition == reference_condition:
                    continue
                
                condition_data = class_data[class_data['infection'] == condition]
                
                if len(ref_data) == 0 or len(condition_data) == 0:
                    continue
                
                # Compare logFC distributions
                ref_logfc = ref_data['logFC'].values
                cond_logfc = condition_data['logFC'].values
                
                # Perform statistical test
                try:
                    statistic, p_value = stats.mannwhitneyu(
                        cond_logfc, ref_logfc, alternative='two-sided'
                    )
                except:
                    p_value = 1.0
                    statistic = 0
                
                # Calculate effect size (Cohen's d)
                pooled_std = np.sqrt(((len(ref_logfc) - 1) * np.var(ref_logfc, ddof=1) + 
                                     (len(cond_logfc) - 1) * np.var(cond_logfc, ddof=1)) / 
                                    (len(ref_logfc) + len(cond_logfc) - 2))
                
                if pooled_std > 0:
                    cohens_d = (np.mean(cond_logfc) - np.mean(ref_logfc)) / pooled_std
                else:
                    cohens_d = 0
                
                comparisons.append({
                    'interaction_type': interaction_type,
                    'contact_class': contact_class,
                    'comparison': f'{reference_condition}_vs_{condition}',
                    'condition': condition,
                    'reference': reference_condition,
                    'n_reference': len(ref_data),
                    'n_condition': len(condition_data),
                    'mean_logFC_reference': np.mean(ref_logfc),
                    'mean_logFC_condition': np.mean(cond_logfc),
                    'log2_fold_change': np.log2(np.mean(cond_logfc) / np.mean(ref_logfc)) 
                                       if np.mean(ref_logfc) != 0 else 0,
                    'p_value': p_value,
                    'effect_size': cohens_d,
                    'statistic': statistic
                })
    
    comparisons_df = pd.DataFrame(comparisons)
    
    if len(comparisons_df) > 0:
        # Apply FDR correction
        if len(comparisons_df) > 1:
            _, comparisons_df['fdr'], _, _ = multipletests(
                comparisons_df['p_value'], method='fdr_bh'
            )
        else:
            comparisons_df['fdr'] = comparisons_df['p_value']
        
        comparisons_df['significant'] = comparisons_df['fdr'] < 0.05
        
        print(f"Performed {len(comparisons_df)} comparisons")
        print(f"Significant results: {sum(comparisons_df['significant'])}")
    
    return comparisons_df

def analyze_distance_dependence(enhancer_interactions):
    """
    Analyze distance-dependent changes in enhancer interactions.
    """
    print("Analyzing distance dependence...")
    
    # Define distance bins
    distance_bins = [0, 10000, 50000, 100000, 500000, 1000000, np.inf]
    distance_labels = ['<10kb', '10-50kb', '50-100kb', '100-500kb', '500kb-1Mb', '>1Mb']
    
    # Only analyze cis interactions
    cis_interactions = enhancer_interactions[
        (enhancer_interactions['distance'] > 0) & 
        (enhancer_interactions['distance'] != -1)
    ].copy()
    
    if len(cis_interactions) == 0:
        print("No cis interactions found for distance analysis")
        return pd.DataFrame()
    
    # Add distance bins
    cis_interactions['distance_bin'] = pd.cut(
        cis_interactions['distance'], 
        bins=distance_bins, 
        labels=distance_labels,
        include_lowest=True
    )
    
    distance_analysis = []
    
    # Group by interaction type, contact class, and distance bin
    for interaction_type in cis_interactions['interaction_type'].unique():
        type_data = cis_interactions[cis_interactions['interaction_type'] == interaction_type]
        
        for contact_class in type_data['contact_class'].unique():
            class_data = type_data[type_data['contact_class'] == contact_class]
            
            for dist_bin in distance_labels:
                bin_data = class_data[class_data['distance_bin'] == dist_bin]
                
                if len(bin_data) < 5:  # Minimum samples
                    continue
                
                # Group by condition
                for condition in bin_data['infection'].unique():
                    condition_data = bin_data[bin_data['infection'] == condition]
                    
                    distance_analysis.append({
                        'interaction_type': interaction_type,
                        'contact_class': contact_class,
                        'condition': condition,
                        'distance_bin': dist_bin,
                        'n_interactions': len(condition_data),
                        'mean_logFC': condition_data['logFC'].mean(),
                        'median_logFC': condition_data['logFC'].median(),
                        'std_logFC': condition_data['logFC'].std(),
                        'mean_distance': condition_data['distance'].mean()
                    })
    
    distance_df = pd.DataFrame(distance_analysis)
    
    if len(distance_df) > 0:
        print(f"Distance analysis completed for {len(distance_df)} combinations")
    
    return distance_df

def create_comprehensive_plots(comparisons_df, distance_df, enhancer_interactions, output_prefix):
    """Create comprehensive visualizations for enhancer analysis"""
    
    if len(comparisons_df) == 0:
        print("No comparison data available for plotting")
        return
    
    fig = plt.figure(figsize=(20, 16))
    gs = fig.add_gridspec(4, 3, hspace=0.4, wspace=0.3)
    
    # Plot 1: Overview of interaction changes
    ax1 = fig.add_subplot(gs[0, :])
    
    # Create a pivot table for visualization
    if len(comparisons_df) > 0:
        pivot_data = comparisons_df.pivot_table(
            index='contact_class', 
            columns=['interaction_type', 'condition'], 
            values='log2_fold_change',
            aggfunc='mean'
        )
        
        if not pivot_data.empty:
            im = ax1.imshow(pivot_data.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
            ax1.set_xticks(range(len(pivot_data.columns)))
            ax1.set_xticklabels([f"{col[0]}\n{col[1]}" for col in pivot_data.columns], rotation=45)
            ax1.set_yticks(range(len(pivot_data.index)))
            ax1.set_yticklabels(pivot_data.index)
            ax1.set_title('Log2 Fold Change Heatmap')
            plt.colorbar(im, ax=ax1)
            
            # Add significance markers
            for i, contact_class in enumerate(pivot_data.index):
                for j, (int_type, condition) in enumerate(pivot_data.columns):
                    sig_row = comparisons_df[
                        (comparisons_df['contact_class'] == contact_class) &
                        (comparisons_df['interaction_type'] == int_type) &
                        (comparisons_df['condition'] == condition) &
                        (comparisons_df['significant'])
                    ]
                    if len(sig_row) > 0:
                        ax1.text(j, i, '*', ha='center', va='center', 
                                color='white', fontsize=16, fontweight='bold')
    
    # Plot 2: Bar plot of significant changes
    ax2 = fig.add_subplot(gs[1, 0])
    if len(comparisons_df) > 0:
        sig_counts = comparisons_df.groupby(['interaction_type', 'condition'])['significant'].sum()
        sig_counts.plot(kind='bar', ax=ax2)
        ax2.set_title('Number of Significant Changes')
        ax2.set_ylabel('Count')
        ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: Effect size distribution
    ax3 = fig.add_subplot(gs[1, 1])
    if len(comparisons_df) > 0:
        ax3.hist(comparisons_df['effect_size'], bins=20, alpha=0.7, edgecolor='black')
        ax3.axvline(0, color='red', linestyle='--', alpha=0.7)
        ax3.set_xlabel('Effect Size (Cohen\'s d)')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Effect Size Distribution')
    
    # Plot 4: P-value distribution
    ax4 = fig.add_subplot(gs[1, 2])
    if len(comparisons_df) > 0:
        ax4.hist(comparisons_df['p_value'], bins=20, alpha=0.7, edgecolor='black')
        ax4.axvline(0.05, color='red', linestyle='--', alpha=0.7, label='p=0.05')
        ax4.set_xlabel('P-value')
        ax4.set_ylabel('Frequency')
        ax4.set_title('P-value Distribution')
        ax4.legend()
    
    # Plot 5: Distance-dependent changes
    ax5 = fig.add_subplot(gs[2, :])
    if len(distance_df) > 0:
        for interaction_type in distance_df['interaction_type'].unique():
            type_data = distance_df[distance_df['interaction_type'] == interaction_type]
            
            for contact_class in type_data['contact_class'].unique():
                class_data = type_data[type_data['contact_class'] == contact_class]
                
                # Average across conditions
                avg_data = class_data.groupby('distance_bin')['mean_logFC'].mean()
                
                if len(avg_data) > 0:
                    x_pos = range(len(avg_data))
                    ax5.plot(x_pos, avg_data.values, marker='o', 
                            label=f'{interaction_type} {contact_class}', linewidth=2)
        
        if len(distance_df['distance_bin'].unique()) > 0:
            ax5.set_xticks(range(len(distance_df['distance_bin'].unique())))
            ax5.set_xticklabels(sorted(distance_df['distance_bin'].unique()))
            ax5.set_xlabel('Distance Bin')
            ax5.set_ylabel('Mean Log2 Fold Change')
            ax5.set_title('Distance-Dependent Changes')
            ax5.axhline(0, color='black', linestyle='--', alpha=0.3)
            ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Plot 6: Volcano plot for most significant comparison
    ax6 = fig.add_subplot(gs[3, :])
    if len(comparisons_df) > 0:
        # Plot effect size vs -log10(p-value)
        comparisons_df['neg_log10_p'] = -np.log10(comparisons_df['p_value'].clip(lower=1e-10))
        
        # Non-significant points
        non_sig = comparisons_df[~comparisons_df['significant']]
        ax6.scatter(non_sig['effect_size'], non_sig['neg_log10_p'], 
                   c='gray', alpha=0.6, s=50)
        
        # Significant points
        sig = comparisons_df[comparisons_df['significant']]
        if len(sig) > 0:
            colors = {'E-TSS': 'red', 'E-E': 'blue'}
            for int_type in sig['interaction_type'].unique():
                type_data = sig[sig['interaction_type'] == int_type]
                color = colors.get(int_type, 'green')
                ax6.scatter(type_data['effect_size'], type_data['neg_log10_p'],
                           c=color, alpha=0.8, s=100, label=int_type)
            
            # Label significant points
            for _, row in sig.iterrows():
                ax6.annotate(f"{row['contact_class']}\n{row['condition']}", 
                           (row['effect_size'], row['neg_log10_p']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, alpha=0.8)
        
        ax6.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
        ax6.axvline(0, color='black', linestyle='-', alpha=0.3)
        ax6.set_xlabel('Effect Size (Cohen\'s d)')
        ax6.set_ylabel('-Log10 P-value')
        ax6.set_title('Volcano Plot: Effect Size vs Significance')
        ax6.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_enhancer_analysis_comprehensive.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved comprehensive plot: {output_prefix}_enhancer_analysis_comprehensive.pdf")

def create_summary_tables(comparisons_df, distance_df, enhancer_interactions, output_prefix):
    """Create summary tables"""
    
    # Save main results
    if len(comparisons_df) > 0:
        comparisons_df.to_csv(f"{output_prefix}_comparisons.tsv", sep='\t', index=False)
        print(f"Saved comparisons: {output_prefix}_comparisons.tsv")
        
        # Save significant results
        sig_results = comparisons_df[comparisons_df['significant']]
        if len(sig_results) > 0:
            sig_results.to_csv(f"{output_prefix}_significant_changes.tsv", sep='\t', index=False)
            print(f"Saved significant changes: {output_prefix}_significant_changes.tsv")
    
    # Save distance analysis
    if len(distance_df) > 0:
        distance_df.to_csv(f"{output_prefix}_distance_analysis.tsv", sep='\t', index=False)
        print(f"Saved distance analysis: {output_prefix}_distance_analysis.tsv")
    
    # Save detailed enhancer interactions
    if len(enhancer_interactions) > 0:
        enhancer_interactions.to_csv(f"{output_prefix}_enhancer_interactions.tsv", sep='\t', index=False)
        print(f"Saved enhancer interactions: {output_prefix}_enhancer_interactions.tsv")
    
    # Create summary statistics
    summary_stats = []
    
    if len(enhancer_interactions) > 0:
        summary_stats.append(f"Total enhancer interactions: {len(enhancer_interactions)}")
        summary_stats.append(f"E-E interactions: {sum(enhancer_interactions['interaction_type'] == 'E-E')}")
        summary_stats.append(f"E-TSS interactions: {sum(enhancer_interactions['interaction_type'] == 'E-TSS')}")
        
        if 'contact_class' in enhancer_interactions.columns:
            class_counts = enhancer_interactions['contact_class'].value_counts()
            for class_name, count in class_counts.items():
                summary_stats.append(f"{class_name} interactions: {count}")
    
    if len(comparisons_df) > 0:
        summary_stats.append(f"\nTotal comparisons: {len(comparisons_df)}")
        summary_stats.append(f"Significant comparisons: {sum(comparisons_df['significant'])}")
        
        if sum(comparisons_df['significant']) > 0:
            sig_df = comparisons_df[comparisons_df['significant']]
            summary_stats.append(f"Mean effect size (significant): {sig_df['effect_size'].mean():.3f}")
            summary_stats.append(f"Mean log2FC (significant): {sig_df['log2_fold_change'].mean():.3f}")
    
    # Save summary
    with open(f"{output_prefix}_summary.txt", 'w') as f:
        f.write("Enhancer Class Analysis Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write("\n".join(summary_stats))
    
    print(f"Saved summary: {output_prefix}_summary.txt")

def main():
    parser = argparse.ArgumentParser(description='Analyze enhancer class interactions from differential interaction data')
    parser.add_argument('--enhancers', required=True, help='Enhancer BED file')
    parser.add_argument('--classification', help='Optional enhancer classification file')
    parser.add_argument('--interactions', required=True, help='CSV file with differential interactions')
    parser.add_argument('--fdr_threshold', type=float, default=0.05, help='FDR threshold for significant interactions')
    parser.add_argument('--reference_condition', default='DOX', help='Reference condition for comparisons')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    print("Starting enhancer class analysis...")
    print(f"Input files:")
    print(f"  Enhancers: {args.enhancers}")
    print(f"  Interactions: {args.interactions}")
    print(f"  Classification: {args.classification}")
    print(f"  Reference condition: {args.reference_condition}")
    print(f"  FDR threshold: {args.fdr_threshold}")
    
    # Step 1: Classify enhancers
    enhancers = classify_enhancers(args.enhancers, args.classification)
    
    # Step 2: Load differential interactions
    interactions = load_differential_interactions(args.interactions, args.fdr_threshold)
    
    if len(interactions) == 0:
        print("No significant interactions found! Exiting...")
        return
    
    # Step 3: Find enhancer interactions
    enhancer_interactions = find_enhancer_interactions(enhancers, interactions)
    
    if len(enhancer_interactions) == 0:
        print("No enhancer interactions found! Check your enhancer annotations.")
        return
    
    # Step 4: Compare conditions
    comparisons_df = compare_conditions(enhancer_interactions, args.reference_condition)
    
    # Step 5: Analyze distance dependence
    distance_df = analyze_distance_dependence(enhancer_interactions)
    
    # Step 6: Create visualizations
    create_comprehensive_plots(comparisons_df, distance_df, enhancer_interactions, args.output_prefix)
    
    # Step 7: Save results
    create_summary_tables(comparisons_df, distance_df, enhancer_interactions, args.output_prefix)
    
    print(f"\nAnalysis complete! Results saved to {args.output_prefix}_*")
    
    # Print final summary
    print("\nFinal Summary:")
    print(f"  Total enhancer interactions found: {len(enhancer_interactions)}")
    if len(comparisons_df) > 0:
        print(f"  Total comparisons performed: {len(comparisons_df)}")
        print(f"  Significant changes: {sum(comparisons_df['significant'])}")
        
        if sum(comparisons_df['significant']) > 0:
            sig_results = comparisons_df[comparisons_df['significant']]
            print(f"  Most significant change: {sig_results.loc[sig_results['p_value'].idxmin(), 'comparison']} "
                  f"({sig_results.loc[sig_results['p_value'].idxmin(), 'contact_class']} "
                  f"{sig_results.loc[sig_results['p_value'].idxmin(), 'interaction_type']})")

if __name__ == '__main__':
    main()