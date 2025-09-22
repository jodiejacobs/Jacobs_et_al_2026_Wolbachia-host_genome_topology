#!/usr/bin/env python3
"""
Complete enhancer class analysis script for Wolbachia chromatin interactions.
Analyzes enhancer interactions using pre-identified differential interactions from diffHic,
with null model comparison for robust FDR correction.
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
    """
    print("Loading enhancer annotations...")
    
    # Try to read the enhancer file - handle different formats
    try:
        # First try to read all columns to see what we have
        enhancers = pd.read_csv(enhancer_file, sep='\t', header=None)
        print(f"Enhancer file has {len(enhancers.columns)} columns")
        print(f"First few rows:")
        print(enhancers.head())
        
        # Check if we have 7 columns (extended BED with classification)
        if len(enhancers.columns) >= 7:
            enhancers.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'class']
            print("Detected extended BED format with classification in 7th column")
        elif len(enhancers.columns) >= 6:
            enhancers.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        elif len(enhancers.columns) >= 4:
            enhancers.columns = ['chrom', 'start', 'end', 'name'] + [f'col_{i}' for i in range(4, len(enhancers.columns))]
        else:
            enhancers.columns = ['chrom', 'start', 'end'] + [f'col_{i}' for i in range(3, len(enhancers.columns))]
            
        print(f"Assigned column names: {list(enhancers.columns)}")
        
    except Exception as e:
        print(f"Error reading enhancer file: {e}")
        return pd.DataFrame()
    
    # Ensure we have required columns
    required_cols = ['chrom', 'start', 'end']
    missing_cols = [col for col in required_cols if col not in enhancers.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Fix data types for coordinates
    enhancers['start'] = pd.to_numeric(enhancers['start'], errors='coerce')
    enhancers['end'] = pd.to_numeric(enhancers['end'], errors='coerce')
    
    # Add name column if missing
    if 'name' not in enhancers.columns:
        enhancers['name'] = enhancers.apply(
            lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1
        )
    else:
        # Ensure name is string
        enhancers['name'] = enhancers['name'].astype(str)
    
    # Add score column if missing or fix it
    if 'score' not in enhancers.columns:
        enhancers['score'] = 0
    else:
        # Handle score column - convert to numeric or set to 0
        enhancers['score'] = pd.to_numeric(enhancers['score'], errors='coerce').fillna(0)
    
    # Add strand if missing
    if 'strand' not in enhancers.columns:
        enhancers['strand'] = '.'
    else:
        enhancers['strand'] = enhancers['strand'].astype(str).fillna('.')
    
    print(f"Data types after processing:")
    print(enhancers.dtypes)
    
    # Check if classification is already in the BED file
    if 'class' in enhancers.columns:
        print("Found classification column in BED file - using those annotations")
        # Clean up the classification values
        enhancers['class'] = enhancers['class'].astype(str).str.lower().str.strip()
        # Map common variations
        class_mapping = {
            'hk': 'housekeeping',
            'housekeeping': 'housekeeping',
            'dev': 'developmental', 
            'developmental': 'developmental',
            'tissue': 'developmental',
            'tissue_specific': 'developmental'
        }
        enhancers['class'] = enhancers['class'].map(class_mapping).fillna('unknown')
        
    elif classification_file and Path(classification_file).exists():
        # Load actual classification from separate file
        print(f"Loading enhancer classification from {classification_file}")
        try:
            # Try to read with header first
            classification = pd.read_csv(classification_file, sep='\t')
            
            # Check if it has proper column names
            if 'name' in classification.columns and 'class' in classification.columns:
                enhancers = enhancers.merge(classification[['name', 'class']], on='name', how='left')
            else:
                # Assume first column is name, second is class (no header)
                print("Classification file appears to have no header - assuming first column is name, second is class")
                classification = pd.read_csv(classification_file, sep='\t', header=None, names=['name', 'class'])
                
                # Check if we have the right number of enhancers
                print(f"Found {len(classification)} classifications for {len(enhancers)} enhancers")
                
                # If the classification file has the same number of rows as enhancers, assume it's in the same order
                if len(classification) == len(enhancers):
                    print("Same number of entries - assuming same order")
                    enhancers['class'] = classification['class'].values
                else:
                    # Try to merge by name
                    enhancers = enhancers.merge(classification[['name', 'class']], on='name', how='left')
            
            enhancers['class'].fillna('unknown', inplace=True)
            
        except Exception as e:
            print(f"Error loading classification file: {e}")
            print("Using heuristic classification instead")
            # Fall back to heuristic
            housekeeping_markers = ['ubiq', 'house', 'const', 'rp', 'ef1', 'gapdh', 'actb', 'rpl', 'rps']
            enhancers['class'] = enhancers['name'].apply(
                lambda x: 'housekeeping' if any(marker in str(x).lower() for marker in housekeeping_markers)
                else 'developmental'
            )
    else:
        # Use heuristic based on enhancer names
        print("Using heuristic classification")
        housekeeping_markers = ['ubiq', 'house', 'const', 'rp', 'ef1', 'gapdh', 'actb', 'rpl', 'rps']
        
        enhancers['class'] = enhancers['name'].apply(
            lambda x: 'housekeeping' if any(marker in str(x).lower() for marker in housekeeping_markers)
            else 'developmental'
        )
    
    print(f"Final data types:")
    print(enhancers.dtypes)
    print(f"Sample of processed data:")
    print(enhancers.head())
    
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
    
    # Add logFC column if missing
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

def load_null_model(null_file):
    """
    Load null model results for comparison.
    """
    print(f"Loading null model from {null_file}")
    
    null_data = pd.read_csv(null_file)
    
    # Remove the unnamed index column if it exists
    if null_data.columns[0].startswith('Unnamed') or null_data.columns[0] == '':
        null_data = null_data.drop(columns=null_data.columns[0])
    
    print(f"Loaded {len(null_data)} null interactions")
    print(f"Null model columns: {list(null_data.columns)}")
    
    # Check for required columns
    required_cols = ['logFC', 'logCPM', 'PValue', 'FDR']
    missing_cols = [col for col in required_cols if col not in null_data.columns]
    if missing_cols:
        print(f"Warning: Missing columns in null model: {missing_cols}")
    
    return null_data

def find_enhancer_interactions(enhancers, interactions, interaction_type='both'):
    """
    Find interactions involving enhancers.
    """
    print(f"Finding enhancer interactions...")
    
    # Create proper BED format for bedtools (6 columns: chr, start, end, name, score, strand)
    enhancers_bed_df = enhancers[['chrom', 'start', 'end', 'name', 'score', 'strand']].copy()
    
    # Ensure numeric columns are correct types
    enhancers_bed_df['start'] = enhancers_bed_df['start'].astype(int)
    enhancers_bed_df['end'] = enhancers_bed_df['end'].astype(int)
    
    # Handle score column - set to 0 if it's not numeric
    try:
        enhancers_bed_df['score'] = pd.to_numeric(enhancers_bed_df['score'], errors='coerce').fillna(0).astype(int)
    except:
        enhancers_bed_df['score'] = 0
    
    # Ensure name is string
    enhancers_bed_df['name'] = enhancers_bed_df['name'].astype(str)
    
    # Ensure strand is string and set default if needed
    if enhancers_bed_df['strand'].isna().any():
        enhancers_bed_df['strand'] = enhancers_bed_df['strand'].fillna('.')
    enhancers_bed_df['strand'] = enhancers_bed_df['strand'].astype(str)
    
    # Create BedTool object
    try:
        enhancers_bed = pybedtools.BedTool.from_dataframe(enhancers_bed_df)
        print(f"Created BedTool with {len(enhancers_bed)} enhancers")
    except Exception as e:
        print(f"Error creating BedTool from enhancers: {e}")
        print("Enhancer dataframe info:")
        print(enhancers_bed_df.dtypes)
        print(enhancers_bed_df.head())
        return pd.DataFrame()
    
    enhancer_interactions = []
    
    print(f"Processing {len(interactions)} interactions...")
    
    for idx, interaction in interactions.iterrows():
        if idx % 10000 == 0:
            print(f"  Processed {idx}/{len(interactions)} interactions...")
        
        try:
            # Create BedTool objects for interaction anchors (ensure proper format)
            anchor1_str = f"{interaction['chr1']}\t{int(interaction['start1'])}\t{int(interaction['end1'])}"
            anchor2_str = f"{interaction['chr2']}\t{int(interaction['start2'])}\t{int(interaction['end2'])}"
            
            anchor1_bed = pybedtools.BedTool(anchor1_str, from_string=True)
            anchor2_bed = pybedtools.BedTool(anchor2_str, from_string=True)
            
            # Find overlaps with enhancers
            anchor1_overlaps = anchor1_bed.intersect(enhancers_bed, wa=True, wb=True)
            anchor2_overlaps = anchor2_bed.intersect(enhancers_bed, wa=True, wb=True)
            
            anchor1_enhancers = []
            anchor2_enhancers = []
            
            # Extract enhancer information for anchor 1
            for overlap in anchor1_overlaps:
                if len(overlap.fields) >= 9:  # 3 (anchor) + 6 (enhancer) fields
                    enh_name = overlap.fields[6]  # 4th field of enhancer (name)
                    # Get enhancer class from original dataframe
                    enh_class = enhancers[enhancers['name'] == enh_name]['class'].iloc[0] if len(enhancers[enhancers['name'] == enh_name]) > 0 else 'unknown'
                    anchor1_enhancers.append({
                        'name': enh_name,
                        'class': enh_class
                    })
            
            # Extract enhancer information for anchor 2
            for overlap in anchor2_overlaps:
                if len(overlap.fields) >= 9:
                    enh_name = overlap.fields[6]
                    enh_class = enhancers[enhancers['name'] == enh_name]['class'].iloc[0] if len(enhancers[enhancers['name'] == enh_name]) > 0 else 'unknown'
                    anchor2_enhancers.append({
                        'name': enh_name,
                        'class': enh_class
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
                # E-TSS interaction
                enhancer_anchor = anchor1_enhancers if len(anchor1_enhancers) > 0 else anchor2_enhancers
                
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
        
        except Exception as e:
            if idx < 10:  # Only print first few errors to avoid spam
                print(f"Warning: Error processing interaction {idx}: {e}")
            continue
    
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

def compare_conditions_with_null(enhancer_interactions, null_data, reference_condition='DOX'):
    """
    Compare enhancer interactions between conditions using null model for FDR correction.
    """
    print(f"Comparing conditions with null model (reference: {reference_condition})...")
    
    if 'infection' not in enhancer_interactions.columns:
        print("Warning: No infection/condition column found")
        return pd.DataFrame()
    
    conditions = enhancer_interactions['infection'].unique()
    print(f"Available conditions: {conditions}")
    
    if reference_condition not in conditions:
        print(f"Warning: Reference condition '{reference_condition}' not found in data")
        reference_condition = conditions[0]
        print(f"Using '{reference_condition}' as reference")
    
    comparisons = []
    
    # Prepare null distribution
    null_logfc = null_data['logFC'].values if 'logFC' in null_data.columns else np.zeros(len(null_data))
    null_pvalues = null_data['PValue'].values if 'PValue' in null_data.columns else np.ones(len(null_data))
    
    print(f"Null model stats: mean logFC = {np.mean(null_logfc):.4f}, mean p-value = {np.mean(null_pvalues):.4f}")
    
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
                
                # Calculate fold change
                mean_ref = np.mean(ref_logfc)
                mean_cond = np.mean(cond_logfc)
                log2_fold_change = mean_cond - mean_ref
                
                # Statistical tests
                try:
                    stat_cond_ref, p_cond_ref = stats.mannwhitneyu(
                        cond_logfc, ref_logfc, alternative='two-sided'
                    )
                except:
                    p_cond_ref = 1.0
                    stat_cond_ref = 0
                
                try:
                    stat_cond_null, p_cond_null = stats.mannwhitneyu(
                        cond_logfc, null_logfc, alternative='two-sided'
                    )
                except:
                    p_cond_null = 1.0
                    stat_cond_null = 0
                
                try:
                    stat_ref_null, p_ref_null = stats.mannwhitneyu(
                        ref_logfc, null_logfc, alternative='two-sided'
                    )
                except:
                    p_ref_null = 1.0
                    stat_ref_null = 0
                
                # Calculate effect sizes
                def cohens_d(x1, x2):
                    if len(x1) == 0 or len(x2) == 0:
                        return 0
                    pooled_std = np.sqrt(((len(x1) - 1) * np.var(x1, ddof=1) + 
                                         (len(x2) - 1) * np.var(x2, ddof=1)) / 
                                        (len(x1) + len(x2) - 2))
                    if pooled_std > 0:
                        return (np.mean(x1) - np.mean(x2)) / pooled_std
                    return 0
                
                effect_size_cond_ref = cohens_d(cond_logfc, ref_logfc)
                effect_size_cond_null = cohens_d(cond_logfc, null_logfc)
                effect_size_ref_null = cohens_d(ref_logfc, null_logfc)
                
                # Calculate empirical p-values
                emp_p_cond = np.mean(np.abs(null_logfc) >= np.abs(mean_cond)) if len(null_logfc) > 0 else 1.0
                emp_p_ref = np.mean(np.abs(null_logfc) >= np.abs(mean_ref)) if len(null_logfc) > 0 else 1.0
                
                comparisons.append({
                    'interaction_type': interaction_type,
                    'contact_class': contact_class,
                    'comparison': f'{reference_condition}_vs_{condition}',
                    'condition': condition,
                    'reference': reference_condition,
                    'n_reference': len(ref_data),
                    'n_condition': len(condition_data),
                    'n_null': len(null_logfc),
                    'mean_logFC_reference': mean_ref,
                    'mean_logFC_condition': mean_cond,
                    'mean_logFC_null': np.mean(null_logfc),
                    'log2_fold_change': log2_fold_change,
                    'p_value_cond_vs_ref': p_cond_ref,
                    'p_value_cond_vs_null': p_cond_null,
                    'p_value_ref_vs_null': p_ref_null,
                    'empirical_p_condition': emp_p_cond,
                    'empirical_p_reference': emp_p_ref,
                    'effect_size_cond_ref': effect_size_cond_ref,
                    'effect_size_cond_null': effect_size_cond_null,
                    'effect_size_ref_null': effect_size_ref_null,
                    'statistic_cond_ref': stat_cond_ref,
                    'statistic_cond_null': stat_cond_null,
                    'statistic_ref_null': stat_ref_null
                })
    
    comparisons_df = pd.DataFrame(comparisons)
    
    if len(comparisons_df) > 0:
        # Apply FDR correction
        p_value_columns = ['p_value_cond_vs_ref', 'p_value_cond_vs_null', 
                          'p_value_ref_vs_null', 'empirical_p_condition', 'empirical_p_reference']
        
        for col in p_value_columns:
            if len(comparisons_df) > 1:
                _, fdr_col, _, _ = multipletests(comparisons_df[col], method='fdr_bh')
                comparisons_df[f'fdr_{col}'] = fdr_col
            else:
                comparisons_df[f'fdr_{col}'] = comparisons_df[col]
        
        # Determine significance
        comparisons_df['significant_vs_ref'] = comparisons_df['fdr_p_value_cond_vs_ref'] < 0.05
        comparisons_df['significant_vs_null'] = comparisons_df['fdr_p_value_cond_vs_null'] < 0.05
        comparisons_df['significant_empirical'] = comparisons_df['fdr_empirical_p_condition'] < 0.05
        
        # Overall significance
        comparisons_df['significant'] = (
            comparisons_df['significant_vs_null'] & 
            (np.abs(comparisons_df['effect_size_cond_null']) > 0.2)
        )
        
        print(f"Performed {len(comparisons_df)} comparisons")
        print(f"Significant vs reference: {sum(comparisons_df['significant_vs_ref'])}")
        print(f"Significant vs null: {sum(comparisons_df['significant_vs_null'])}")
        print(f"Overall significant: {sum(comparisons_df['significant'])}")
    
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
                
                if len(bin_data) < 5:
                    continue
                
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

def create_comprehensive_plots(comparisons_df, distance_df, enhancer_interactions, null_data, output_prefix):
    """Create comprehensive visualizations"""
    
    if len(comparisons_df) == 0:
        print("No comparison data available for plotting")
        return
    
    fig = plt.figure(figsize=(24, 20))
    gs = fig.add_gridspec(5, 4, hspace=0.4, wspace=0.3)
    
    # Plot 1: Null distribution vs observed distributions
    ax1 = fig.add_subplot(gs[0, :2])
    if len(null_data) > 0 and 'logFC' in null_data.columns:
        null_logfc = null_data['logFC'].values
        ax1.hist(null_logfc, bins=50, alpha=0.7, label='Null', color='gray', density=True)
        
        colors = {'DOX': 'blue', 'wMel': 'red', 'wRi': 'green', 'wWil': 'orange'}
        for condition in enhancer_interactions['infection'].unique():
            cond_data = enhancer_interactions[enhancer_interactions['infection'] == condition]
            if len(cond_data) > 0:
                ax1.hist(cond_data['logFC'], bins=30, alpha=0.6, 
                        label=condition, color=colors.get(condition, 'black'), 
                        density=True, histtype='step', linewidth=2)
        
        ax1.axvline(0, color='black', linestyle='--', alpha=0.5)
        ax1.set_xlabel('Log Fold Change')
        ax1.set_ylabel('Density')
        ax1.set_title('Distribution Comparison: Null vs Observed')
        ax1.legend()
    
    # Plot 2: Significance comparison
    ax2 = fig.add_subplot(gs[0, 2:])
    if len(comparisons_df) > 0:
        p_columns = ['p_value_cond_vs_ref', 'p_value_cond_vs_null', 'empirical_p_condition']
        p_labels = ['vs Reference', 'vs Null', 'Empirical']
        
        sig_counts = []
        for col in p_columns:
            if f'fdr_{col}' in comparisons_df.columns:
                sig_counts.append(sum(comparisons_df[f'fdr_{col}'] < 0.05))
            else:
                sig_counts.append(0)
        
        bars = ax2.bar(p_labels, sig_counts, color=['skyblue', 'lightcoral', 'lightgreen'])
        ax2.set_ylabel('Number of Significant Results')
        ax2.set_title('Significance by Test Type (FDR < 0.05)')
        
        for bar, count in zip(bars, sig_counts):
            if count > 0:
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                        str(count), ha='center', va='bottom')
    
    # Plot 3: Effect size comparison
    ax3 = fig.add_subplot(gs[1, :2])
    if len(comparisons_df) > 0:
        effect_columns = ['effect_size_cond_ref', 'effect_size_cond_null', 'effect_size_ref_null']
        effect_labels = ['Cond vs Ref', 'Cond vs Null', 'Ref vs Null']
        
        effect_data = []
        for col in effect_columns:
            if col in comparisons_df.columns:
                effect_data.append(comparisons_df[col].values)
            else:
                effect_data.append([])
        
        bp = ax3.boxplot([data for data in effect_data if len(data) > 0], 
                        labels=[label for i, label in enumerate(effect_labels) if len(effect_data[i]) > 0],
                        patch_artist=True)
        
        colors = ['lightblue', 'lightcoral', 'lightgreen']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        
        ax3.axhline(0, color='black', linestyle='--', alpha=0.5)
        ax3.set_ylabel('Effect Size (Cohen\'s d)')
        ax3.set_title('Effect Size Distributions')
    
    # Plot 4: Volcano plot
    ax4 = fig.add_subplot(gs[1, 2:])
    if len(comparisons_df) > 0 and 'empirical_p_condition' in comparisons_df.columns:
        comparisons_df['neg_log10_emp_p'] = -np.log10(
            comparisons_df['empirical_p_condition'].clip(lower=1e-10)
        )
        
        # Plot non-significant
        non_sig = comparisons_df[~comparisons_df['significant_empirical']]
        if len(non_sig) > 0:
            ax4.scatter(non_sig['effect_size_cond_null'], non_sig['neg_log10_emp_p'], 
                       c='gray', alpha=0.6, s=50)
        
        # Plot significant
        sig = comparisons_df[comparisons_df['significant_empirical']]
        if len(sig) > 0:
            colors = {'E-TSS': 'red', 'E-E': 'blue'}
            for int_type in sig['interaction_type'].unique():
                type_data = sig[sig['interaction_type'] == int_type]
                color = colors.get(int_type, 'green')
                ax4.scatter(type_data['effect_size_cond_null'], type_data['neg_log10_emp_p'],
                           c=color, alpha=0.8, s=100, label=int_type)
            
            # Label significant points
            for _, row in sig.iterrows():
                ax4.annotate(f"{row['contact_class']}\n{row['condition']}", 
                           (row['effect_size_cond_null'], row['neg_log10_emp_p']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, alpha=0.8)
        
        ax4.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
        ax4.axvline(0, color='black', linestyle='-', alpha=0.3)
        ax4.set_xlabel('Effect Size vs Null (Cohen\'s d)')
        ax4.set_ylabel('-Log10 Empirical P-value')
        ax4.set_title('Volcano Plot: Effect Size vs Null Model')
        ax4.legend()
    
    # Plot 5: Heatmap
    ax5 = fig.add_subplot(gs[2, :])
    if len(comparisons_df) > 0:
        heatmap_data = comparisons_df.pivot_table(
            index=['interaction_type', 'contact_class'], 
            columns='condition', 
            values='effect_size_cond_null',
            aggfunc='mean'
        )
        
        if not heatmap_data.empty:
            im = ax5.imshow(heatmap_data.values, aspect='auto', cmap='RdBu_r', 
                           vmin=-2, vmax=2)
            
            ax5.set_xticks(range(len(heatmap_data.columns)))
            ax5.set_xticklabels(heatmap_data.columns)
            ax5.set_yticks(range(len(heatmap_data.index)))
            ax5.set_yticklabels([f"{idx[0]}\n{idx[1]}" for idx in heatmap_data.index])
            ax5.set_title('Effect Size vs Null Model Heatmap')
            
            cbar = plt.colorbar(im, ax=ax5)
            cbar.set_label('Effect Size (Cohen\'s d)')
            
            # Add significance markers
            sig_pivot = comparisons_df.pivot_table(
                index=['interaction_type', 'contact_class'], 
                columns='condition', 
                values='significant',
                aggfunc='any'
            )
            
            for i in range(len(heatmap_data.index)):
                for j in range(len(heatmap_data.columns)):
                    if not sig_pivot.empty and sig_pivot.iloc[i, j]:
                        ax5.text(j, i, '*', ha='center', va='center', 
                                color='white', fontsize=16, fontweight='bold')
    
    # Plot 6: Distance-dependent analysis
    ax6 = fig.add_subplot(gs[3, :])
    if len(distance_df) > 0:
        for interaction_type in distance_df['interaction_type'].unique():
            type_data = distance_df[distance_df['interaction_type'] == interaction_type]
            
            for contact_class in type_data['contact_class'].unique():
                class_data = type_data[type_data['contact_class'] == contact_class]
                
                # Average across conditions
                avg_data = class_data.groupby('distance_bin')['mean_logFC'].mean()
                
                if len(avg_data) > 0:
                    x_pos = range(len(avg_data))
                    ax6.plot(x_pos, avg_data.values, marker='o', 
                            label=f'{interaction_type} {contact_class}', linewidth=2)
        
        if len(distance_df['distance_bin'].unique()) > 0:
            ax6.set_xticks(range(len(distance_df['distance_bin'].unique())))
            ax6.set_xticklabels(sorted(distance_df['distance_bin'].unique()))
            ax6.set_xlabel('Distance Bin')
            ax6.set_ylabel('Mean Log2 Fold Change')
            ax6.set_title('Distance-Dependent Changes')
            ax6.axhline(0, color='black', linestyle='--', alpha=0.3)
            ax6.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Plot 7: Summary statistics
    ax7 = fig.add_subplot(gs[4, :])
    ax7.axis('off')
    
    summary_text = "Analysis Summary with Null Model Comparison\n" + "=" * 50 + "\n\n"
    
    # Add basic statistics
    if len(enhancer_interactions) > 0:
        summary_text += f"Total enhancer interactions: {len(enhancer_interactions)}\n"
        summary_text += f"E-E interactions: {sum(enhancer_interactions['interaction_type'] == 'E-E')}\n"
        summary_text += f"E-TSS interactions: {sum(enhancer_interactions['interaction_type'] == 'E-TSS')}\n\n"
    
    # Add null model statistics
    if len(null_data) > 0:
        null_logfc = null_data['logFC'].values if 'logFC' in null_data.columns else []
        if len(null_logfc) > 0:
            summary_text += f"Null model interactions: {len(null_logfc)}\n"
            summary_text += f"Null mean logFC: {np.mean(null_logfc):.4f} ± {np.std(null_logfc):.4f}\n\n"
    
    # Add comparison statistics
    if len(comparisons_df) > 0:
        summary_text += f"Total comparisons: {len(comparisons_df)}\n"
        
        sig_cols = ['significant_vs_ref', 'significant_vs_null', 'significant_empirical', 'significant']
        for col in sig_cols:
            if col in comparisons_df.columns:
                n_sig = sum(comparisons_df[col])
                summary_text += f"{col.replace('_', ' ').title()}: {n_sig}\n"
        
        # Most significant results
        if 'empirical_p_condition' in comparisons_df.columns:
            most_sig = comparisons_df.loc[comparisons_df['empirical_p_condition'].idxmin()]
            summary_text += f"\nMost significant change:\n"
            summary_text += f"  {most_sig['interaction_type']} {most_sig['contact_class']} "
            summary_text += f"({most_sig['condition']})\n"
            summary_text += f"  Effect size: {most_sig.get('effect_size_cond_null', 0):.3f}\n"
            summary_text += f"  Empirical p-value: {most_sig['empirical_p_condition']:.2e}\n"
    
    ax7.text(0.05, 0.95, summary_text, transform=ax7.transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_enhancer_analysis_comprehensive.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved comprehensive plot: {output_prefix}_enhancer_analysis_comprehensive.pdf")

def create_summary_tables(comparisons_df, distance_df, enhancer_interactions, null_data, output_prefix):
    """Create summary tables including null model comparisons"""
    
    # Save main results with null comparisons
    if len(comparisons_df) > 0:
        comparisons_df.to_csv(f"{output_prefix}_comparisons_with_null.tsv", sep='\t', index=False)
        print(f"Saved comparisons: {output_prefix}_comparisons_with_null.tsv")
        
        # Save significant results for different criteria
        sig_criteria = {
            'vs_reference': 'significant_vs_ref',
            'vs_null': 'significant_vs_null', 
            'empirical': 'significant_empirical',
            'overall': 'significant'
        }
        
        for suffix, col in sig_criteria.items():
            if col in comparisons_df.columns:
                sig_results = comparisons_df[comparisons_df[col]]
                if len(sig_results) > 0:
                    sig_results.to_csv(f"{output_prefix}_significant_{suffix}.tsv", sep='\t', index=False)
                    print(f"Saved significant {suffix} results: {output_prefix}_significant_{suffix}.tsv")
    
    # Save distance analysis
    if len(distance_df) > 0:
        distance_df.to_csv(f"{output_prefix}_distance_analysis.tsv", sep='\t', index=False)
        print(f"Saved distance analysis: {output_prefix}_distance_analysis.tsv")
    
    # Save detailed enhancer interactions
    if len(enhancer_interactions) > 0:
        enhancer_interactions.to_csv(f"{output_prefix}_enhancer_interactions.tsv", sep='\t', index=False)
        print(f"Saved enhancer interactions: {output_prefix}_enhancer_interactions.tsv")
    
    # Save null model summary
    if len(null_data) > 0:
        null_summary = {
            'n_null_interactions': len(null_data),
            'null_mean_logFC': null_data['logFC'].mean() if 'logFC' in null_data.columns else 0,
            'null_std_logFC': null_data['logFC'].std() if 'logFC' in null_data.columns else 0,
            'null_mean_pvalue': null_data['PValue'].mean() if 'PValue' in null_data.columns else 1,
            'null_mean_logCPM': null_data['logCPM'].mean() if 'logCPM' in null_data.columns else 0
        }
        
        null_summary_df = pd.DataFrame([null_summary])
        null_summary_df.to_csv(f"{output_prefix}_null_model_summary.tsv", sep='\t', index=False)
        print(f"Saved null model summary: {output_prefix}_null_model_summary.tsv")
    
    # Create comprehensive summary statistics
    summary_stats = []
    
    summary_stats.append("=" * 60)
    summary_stats.append("ENHANCER CLASS ANALYSIS WITH NULL MODEL COMPARISON")
    summary_stats.append("=" * 60)
    summary_stats.append("")
    
    # Data summary
    if len(enhancer_interactions) > 0:
        summary_stats.append("DATA SUMMARY:")
        summary_stats.append(f"  Total enhancer interactions: {len(enhancer_interactions)}")
        summary_stats.append(f"  E-E interactions: {sum(enhancer_interactions['interaction_type'] == 'E-E')}")
        summary_stats.append(f"  E-TSS interactions: {sum(enhancer_interactions['interaction_type'] == 'E-TSS')}")
        
        if 'contact_class' in enhancer_interactions.columns:
            class_counts = enhancer_interactions['contact_class'].value_counts()
            for class_name, count in class_counts.items():
                summary_stats.append(f"  {class_name} interactions: {count}")
        
        if 'infection' in enhancer_interactions.columns:
            condition_counts = enhancer_interactions['infection'].value_counts()
            summary_stats.append(f"\n  By condition:")
            for condition, count in condition_counts.items():
                summary_stats.append(f"    {condition}: {count}")
        
        summary_stats.append("")
    
    # Null model summary
    if len(null_data) > 0:
        summary_stats.append("NULL MODEL SUMMARY:")
        summary_stats.append(f"  Null interactions: {len(null_data)}")
        if 'logFC' in null_data.columns:
            null_logfc = null_data['logFC'].values
            summary_stats.append(f"  Null mean logFC: {np.mean(null_logfc):.4f} ± {np.std(null_logfc):.4f}")
            summary_stats.append(f"  Null logFC range: [{np.min(null_logfc):.4f}, {np.max(null_logfc):.4f}]")
        if 'PValue' in null_data.columns:
            summary_stats.append(f"  Null mean p-value: {null_data['PValue'].mean():.4f}")
        summary_stats.append("")
    
    # Statistical results
    if len(comparisons_df) > 0:
        summary_stats.append("STATISTICAL RESULTS:")
        summary_stats.append(f"  Total comparisons performed: {len(comparisons_df)}")
        
        # Count significant results for each test type
        sig_types = {
            'significant_vs_ref': 'Significant vs Reference',
            'significant_vs_null': 'Significant vs Null Model',
            'significant_empirical': 'Significant (Empirical)',
            'significant': 'Overall Significant'
        }
        
        for col, description in sig_types.items():
            if col in comparisons_df.columns:
                n_sig = sum(comparisons_df[col])
                pct_sig = (n_sig / len(comparisons_df)) * 100
                summary_stats.append(f"  {description}: {n_sig} ({pct_sig:.1f}%)")
        
        # Effect size summary
        if 'effect_size_cond_null' in comparisons_df.columns:
            effect_sizes = comparisons_df['effect_size_cond_null'].values
            summary_stats.append(f"\n  Effect size vs null (Cohen's d):")
            summary_stats.append(f"    Mean: {np.mean(effect_sizes):.3f} ± {np.std(effect_sizes):.3f}")
            summary_stats.append(f"    Range: [{np.min(effect_sizes):.3f}, {np.max(effect_sizes):.3f}]")
            
            # Large effect sizes
            large_effects = sum(np.abs(effect_sizes) > 0.8)
            summary_stats.append(f"    Large effects (|d| > 0.8): {large_effects}")
        
        # Most significant results
        if 'empirical_p_condition' in comparisons_df.columns:
            most_sig_idx = comparisons_df['empirical_p_condition'].idxmin()
            most_sig = comparisons_df.loc[most_sig_idx]
            
            summary_stats.append(f"\n  Most significant result:")
            summary_stats.append(f"    {most_sig['interaction_type']} {most_sig['contact_class']} ({most_sig['condition']})")
            summary_stats.append(f"    Effect size vs null: {most_sig.get('effect_size_cond_null', 0):.3f}")
            summary_stats.append(f"    Empirical p-value: {most_sig['empirical_p_condition']:.2e}")
            if f"fdr_empirical_p_condition" in comparisons_df.columns:
                summary_stats.append(f"    FDR: {most_sig[f'fdr_empirical_p_condition']:.2e}")
        
        summary_stats.append("")
    
    # Distance analysis summary
    if len(distance_df) > 0:
        summary_stats.append("DISTANCE ANALYSIS:")
        distance_bins = distance_df['distance_bin'].value_counts()
        summary_stats.append(f"  Distance bins analyzed: {len(distance_bins)}")
        for bin_name, count in distance_bins.items():
            bin_data = distance_df[distance_df['distance_bin'] == bin_name]
            mean_logfc = bin_data['mean_logFC'].mean()
            summary_stats.append(f"    {bin_name}: {count} comparisons, mean logFC = {mean_logfc:.3f}")
        summary_stats.append("")
    
    # Interpretation
    summary_stats.append("INTERPRETATION:")
    if len(comparisons_df) > 0 and 'significant' in comparisons_df.columns:
        n_significant = sum(comparisons_df['significant'])
        if n_significant > 0:
            summary_stats.append(f"  {n_significant} enhancer interaction changes are statistically significant")
            summary_stats.append("  after controlling for null expectations and multiple testing.")
            
            # Identify patterns
            sig_df = comparisons_df[comparisons_df['significant']]
            if len(sig_df) > 0:
                # Most affected condition
                condition_counts = sig_df['condition'].value_counts()
                if len(condition_counts) > 0:
                    top_condition = condition_counts.index[0]
                    summary_stats.append(f"  Most affected condition: {top_condition} ({condition_counts.iloc[0]} changes)")
                
                # Most affected interaction type
                type_counts = sig_df['interaction_type'].value_counts()
                if len(type_counts) > 0:
                    top_type = type_counts.index[0]
                    summary_stats.append(f"  Most affected interaction type: {top_type} ({type_counts.iloc[0]} changes)")
                
                # Direction of changes
                if 'effect_size_cond_null' in sig_df.columns:
                    positive_effects = sum(sig_df['effect_size_cond_null'] > 0)
                    negative_effects = sum(sig_df['effect_size_cond_null'] < 0)
                    summary_stats.append(f"  Direction: {positive_effects} increases, {negative_effects} decreases")
        else:
            summary_stats.append("  No statistically significant changes detected after")
            summary_stats.append("  controlling for null expectations.")
    
    # Save comprehensive summary
    with open(f"{output_prefix}_comprehensive_summary.txt", 'w') as f:
        f.write("\n".join(summary_stats))
    
    print(f"Saved comprehensive summary: {output_prefix}_comprehensive_summary.txt")

def main():
    parser = argparse.ArgumentParser(description='Complete enhancer class interaction analysis with null model')
    parser.add_argument('--enhancers', required=True, help='Enhancer BED file')
    parser.add_argument('--classification', help='Optional enhancer classification file')
    parser.add_argument('--interactions', required=True, help='CSV file with differential interactions')
    parser.add_argument('--null_model', help='CSV file with null model results for FDR correction')
    parser.add_argument('--fdr_threshold', type=float, default=0.05, help='FDR threshold for significant interactions')
    parser.add_argument('--reference_condition', default='DOX', help='Reference condition for comparisons')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    print("Starting complete enhancer class analysis with null model comparison...")
    print(f"Input files:")
    print(f"  Enhancers: {args.enhancers}")
    print(f"  Interactions: {args.interactions}")
    print(f"  Null model: {args.null_model}")
    print(f"  Classification: {args.classification}")
    print(f"  Reference condition: {args.reference_condition}")
    print(f"  FDR threshold: {args.fdr_threshold}")
    
    # Step 1: Classify enhancers
    try:
        enhancers = classify_enhancers(args.enhancers, args.classification)
    except Exception as e:
        print(f"Error loading enhancers: {e}")
        return
    
    # Step 2: Load differential interactions
    try:
        interactions = load_differential_interactions(args.interactions, args.fdr_threshold)
        if len(interactions) == 0:
            print("No significant interactions found! Exiting...")
            return
    except Exception as e:
        print(f"Error loading interactions: {e}")
        return
    
    # Step 3: Load null model if provided
    null_data = pd.DataFrame()
    if args.null_model:
        try:
            null_data = load_null_model(args.null_model)
            print("Successfully loaded null model for robust statistical comparison")
        except Exception as e:
            print(f"Warning: Could not load null model: {e}")
            print("Proceeding without null model comparison...")
    else:
        print("No null model provided - will perform standard comparison only")
    
    # Step 4: Find enhancer interactions
    try:
        enhancer_interactions = find_enhancer_interactions(enhancers, interactions)
        if len(enhancer_interactions) == 0:
            print("No enhancer interactions found! Check your enhancer annotations.")
            return
    except Exception as e:
        print(f"Error finding enhancer interactions: {e}")
        return
    
    # Step 5: Compare conditions
    try:
        if len(null_data) > 0:
            print("Performing robust statistical comparison with null model...")
            comparisons_df = compare_conditions_with_null(enhancer_interactions, null_data, args.reference_condition)
        else:
            print("Performing standard comparison without null model...")
            # Create a simple comparison function for no null model case
            comparisons_df = pd.DataFrame()  # Would need to implement simple version
            print("Warning: Standard comparison not implemented in this version - please provide null model")
    except Exception as e:
        print(f"Error in statistical comparison: {e}")
        comparisons_df = pd.DataFrame()
    
    # Step 6: Analyze distance dependence
    try:
        distance_df = analyze_distance_dependence(enhancer_interactions)
    except Exception as e:
        print(f"Error in distance analysis: {e}")
        distance_df = pd.DataFrame()
    
    # Step 7: Create visualizations
    try:
        create_comprehensive_plots(comparisons_df, distance_df, enhancer_interactions, null_data, args.output_prefix)
    except Exception as e:
        print(f"Error creating plots: {e}")
    
    # Step 8: Save results
    try:
        create_summary_tables(comparisons_df, distance_df, enhancer_interactions, null_data, args.output_prefix)
    except Exception as e:
        print(f"Error creating summary tables: {e}")
    
    print(f"\nAnalysis complete! Results saved to {args.output_prefix}_*")
    
    # Print final summary
    print("\nFinal Summary:")
    print(f"  Total enhancer interactions found: {len(enhancer_interactions)}")
    if len(comparisons_df) > 0:
        print(f"  Total comparisons performed: {len(comparisons_df)}")
        
        # Print results for different significance criteria
        sig_types = ['significant_vs_ref', 'significant_vs_null', 'significant_empirical', 'significant']
        for sig_type in sig_types:
            if sig_type in comparisons_df.columns:
                n_sig = sum(comparisons_df[sig_type])
                print(f"  {sig_type.replace('_', ' ').title()}: {n_sig}")
        
        # Print most significant result
        if len(null_data) > 0 and 'empirical_p_condition' in comparisons_df.columns and len(comparisons_df) > 0:
            most_sig_idx = comparisons_df['empirical_p_condition'].idxmin()
            most_sig = comparisons_df.loc[most_sig_idx]
            print(f"  Most significant change: {most_sig['interaction_type']} "
                  f"{most_sig['contact_class']} ({most_sig['condition']}) "
                  f"p={most_sig['empirical_p_condition']:.2e}")
    else:
        print("  No statistical comparisons performed")

if __name__ == '__main__':
    main()