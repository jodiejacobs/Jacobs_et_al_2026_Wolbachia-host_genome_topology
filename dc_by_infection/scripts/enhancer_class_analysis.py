#!/usr/bin/env python3
"""
Fixed version of enhancer class analysis script.
Addresses:
1. Proper handling of condition-specific interaction files
2. Correct null model comparison
3. Better error handling for missing condition information
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
import sys
import glob

def load_interactions_by_condition(interactions_pattern, reference_condition='DOX'):
    """
    Load differential interactions from multiple condition-specific files.
    
    Parameters:
    -----------
    interactions_pattern : str
        Pattern to match condition-specific files, e.g., 
        'results/*_vs_DOX_interactions.csv' or a single combined file
    reference_condition : str
        Name of reference condition (default: 'DOX')
    
    Returns:
    --------
    interactions : pd.DataFrame
        Combined interactions with 'condition' column
    reference_condition : str
        Validated reference condition name
    """
    print(f"Loading interactions from pattern: {interactions_pattern}")
    
    # Check if it's a pattern or single file
    if '*' in interactions_pattern:
        # Multiple files
        files = glob.glob(interactions_pattern)
        if not files:
            raise ValueError(f"No files found matching pattern: {interactions_pattern}")
        
        print(f"Found {len(files)} interaction files:")
        for f in files:
            print(f"  - {f}")
        
        all_interactions = []
        for file in files:
            # Extract condition from filename
            # Expects format like: wMel_vs_DOX_interactions.csv
            basename = Path(file).stem
            
            # Try to extract condition name
            if '_vs_' in basename:
                condition = basename.split('_vs_')[0]
            elif 'infection' in basename.lower():
                # Handle format like: infectionJW18wMel_interactions.csv
                condition = basename.split('infection')[-1].split('_')[0]
                condition = condition.replace('JW18', '')
            else:
                print(f"Warning: Could not extract condition from filename: {file}")
                condition = basename
            
            df = pd.read_csv(file)
            df['condition'] = condition
            df['source_file'] = file
            all_interactions.append(df)
            print(f"  Loaded {len(df)} interactions for condition: {condition}")
        
        interactions = pd.concat(all_interactions, ignore_index=True)
        
    else:
        # Single file - must have condition column
        print(f"Loading single file: {interactions_pattern}")
        interactions = pd.read_csv(interactions_pattern)
        
        if 'condition' not in interactions.columns:
            print("WARNING: No 'condition' column found in interaction file!")
            print("Available columns:", interactions.columns.tolist())
            
            # Try to infer condition from other columns
            condition_cols = [col for col in interactions.columns 
                            if any(x in col.lower() for x in ['infection', 'condition', 'sample', 'strain'])]
            
            if condition_cols:
                print(f"Found potential condition columns: {condition_cols}")
                print(f"Using column: {condition_cols[0]}")
                interactions['condition'] = interactions[condition_cols[0]]
            else:
                print("ERROR: Cannot determine condition information!")
                print("Options:")
                print("  1. Add a 'condition' column to your CSV file")
                print("  2. Use separate files for each condition")
                print("  3. Use a wildcard pattern to load multiple files")
                sys.exit(1)
    
    print(f"\nTotal interactions loaded: {len(interactions)}")
    print("Available columns:", interactions.columns.tolist())
    
    # Clean condition names
    interactions['condition_clean'] = interactions['condition'].str.replace('JW18', '').str.strip()
    
    # Standardize condition names
    condition_mapping = {
        'wMel': 'wMel',
        'wRi': 'wRi', 
        'wWil': 'wWil',
        'DOX': 'DOX',
        'dox': 'DOX',
        'Dox': 'DOX'
    }
    
    interactions['condition_clean'] = interactions['condition_clean'].replace(condition_mapping)
    
    available_conditions = interactions['condition_clean'].unique()
    print(f"Conditions found: {available_conditions}")
    
    if reference_condition not in available_conditions:
        print(f"Warning: Reference condition '{reference_condition}' not found")
        print(f"Available conditions: {available_conditions}")
        if len(available_conditions) > 0:
            reference_condition = available_conditions[0]
            print(f"Using '{reference_condition}' as reference")
    
    return interactions, reference_condition


def classify_enhancers(enhancer_file, classification_file=None):
    """
    Load and classify enhancers as housekeeping or developmental.
    """
    print("Loading enhancer annotations...")
    
    try:
        enhancers = pd.read_csv(enhancer_file, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'class'])
        print("Loaded enhancer file with 7 columns (including class)")
    except:
        try:
            enhancers = pd.read_csv(enhancer_file, sep='\t')
            print(f"Loaded enhancer file with header: {enhancers.columns.tolist()}")
        except Exception as e:
            print(f"Error reading enhancer file: {e}")
            return pd.DataFrame()
    
    print(f"Enhancer columns: {enhancers.columns.tolist()}")
    print(f"First few rows:")
    print(enhancers.head())
    
    # Ensure required columns exist
    required_cols = ['chrom', 'start', 'end']
    for col in required_cols:
        if col not in enhancers.columns:
            raise ValueError(f"Required column '{col}' not found in enhancer file")
    
    # Convert coordinates to proper types
    print("Converting coordinate columns...")
    enhancers['start'] = pd.to_numeric(enhancers['start'], errors='coerce')
    enhancers['end'] = pd.to_numeric(enhancers['end'], errors='coerce')
    
    # Remove rows with invalid coordinates
    before_clean = len(enhancers)
    enhancers = enhancers.dropna(subset=['start', 'end'])
    
    # Convert to int
    enhancers['start'] = enhancers['start'].astype(int)
    enhancers['end'] = enhancers['end'].astype(int)
    
    # Remove invalid ranges
    enhancers = enhancers[enhancers['start'] < enhancers['end']]
    enhancers = enhancers[(enhancers['start'] >= 0) & (enhancers['end'] >= 0)]
    
    print(f"Coordinate cleaning: {before_clean} -> {len(enhancers)} enhancers")
    
    if len(enhancers) == 0:
        print("No valid enhancers after coordinate cleaning!")
        return pd.DataFrame()
    
    # Create name column if it doesn't exist
    if 'name' not in enhancers.columns:
        enhancers['name'] = enhancers.apply(lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1)
        print("Created name column from coordinates")
    
    # Handle classification
    if 'class' in enhancers.columns and not enhancers['class'].isna().all():
        print("Using classification from BED file")
        print(f"Unique classes in BED file: {enhancers['class'].unique()}")
        
        # Standardize class names
        enhancers['class'] = enhancers['class'].str.lower().str.strip()
        enhancers['class'] = enhancers['class'].replace({
            'housekeeping': 'housekeeping',
            'hk': 'housekeeping',
            'developmental': 'developmental',
            'dev': 'developmental',
            'tissue_specific': 'developmental',
            'tissue-specific': 'developmental'
        })
        
        enhancers['class'] = enhancers['class'].fillna('developmental')
        
    elif classification_file and Path(classification_file).exists():
        print(f"Loading enhancer classification from {classification_file}")
        try:
            classification = pd.read_csv(classification_file, sep='\t', header=None, 
                                       names=['name', 'class'])
            print(f"Classification file head:")
            print(classification.head())
            
            classification['class'] = classification['class'].str.lower().str.strip()
            classification['class'] = classification['class'].replace({
                'housekeeping': 'housekeeping',
                'hk': 'housekeeping', 
                'developmental': 'developmental',
                'dev': 'developmental',
                'tissue_specific': 'developmental',
                'tissue-specific': 'developmental'
            })
            
            print(f"Classification distribution:")
            print(classification['class'].value_counts())
            
            enhancers = enhancers.drop(columns=['class'], errors='ignore')
            enhancers = enhancers.merge(classification, on='name', how='left')
            merged_count = sum(enhancers['class'].notna())
            print(f"Successfully merged {merged_count} classifications out of {len(enhancers)} enhancers")
            
        except Exception as e:
            print(f"Error loading classification file: {e}")
            if 'class' not in enhancers.columns:
                enhancers = apply_heuristic_classification(enhancers)
    else:
        print("No classification file provided and no class in BED file, using heuristic classification")
        enhancers = apply_heuristic_classification(enhancers)
    
    # Fill missing classifications
    enhancers['class'] = enhancers['class'].fillna('developmental')
    
    print(f"Final classification counts:")
    print(f"  Housekeeping: {sum(enhancers['class'] == 'housekeeping')}")
    print(f"  Developmental: {sum(enhancers['class'] == 'developmental')}")
    print(f"  Other: {sum(~enhancers['class'].isin(['housekeeping', 'developmental']))}")
    
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
    Find enhancer-enhancer and enhancer-TSS interactions using manual overlap detection.
    This is more robust than pybedtools for large datasets.
    """
    print("Finding enhancer interactions...")
    print(f"Interactions: {len(interactions)}")
    print(f"Enhancers: {len(enhancers)}")
    
    # Validate coordinates
    interactions_clean = interactions.copy()
    coord_cols = ['start1', 'end1', 'start2', 'end2']
    for col in coord_cols:
        interactions_clean[col] = pd.to_numeric(interactions_clean[col], errors='coerce')
    
    interactions_clean = interactions_clean.dropna(subset=coord_cols)
    for col in coord_cols:
        interactions_clean[col] = interactions_clean[col].astype(int)
    
    # Ensure start <= end
    mask1 = interactions_clean['start1'] > interactions_clean['end1']
    interactions_clean.loc[mask1, ['start1', 'end1']] = interactions_clean.loc[mask1, ['end1', 'start1']].values
    
    mask2 = interactions_clean['start2'] > interactions_clean['end2']
    interactions_clean.loc[mask2, ['start2', 'end2']] = interactions_clean.loc[mask2, ['end2', 'start2']].values
    
    print(f"Cleaned interactions: {len(interactions)} -> {len(interactions_clean)}")
    
    # Clean enhancers
    enhancers_clean = enhancers.copy()
    enhancers_clean['start'] = pd.to_numeric(enhancers_clean['start'], errors='coerce')
    enhancers_clean['end'] = pd.to_numeric(enhancers_clean['end'], errors='coerce')
    enhancers_clean = enhancers_clean.dropna(subset=['start', 'end'])
    enhancers_clean['start'] = enhancers_clean['start'].astype(int)
    enhancers_clean['end'] = enhancers_clean['end'].astype(int)
    
    mask = enhancers_clean['start'] > enhancers_clean['end']
    enhancers_clean.loc[mask, ['start', 'end']] = enhancers_clean.loc[mask, ['end', 'start']].values
    
    print(f"Cleaned enhancers: {len(enhancers)} -> {len(enhancers_clean)}")
    
    if interactions_clean.empty or enhancers_clean.empty:
        print("No valid data after cleaning!")
        return pd.DataFrame()
    
    # Manual overlap detection
    return find_overlaps_manual(interactions_clean, enhancers_clean)


def find_overlaps_manual(interactions, enhancers):
    """
    Manual overlap detection - more reliable than pybedtools for this use case.
    """
    print("Using manual overlap detection...")
    
    enh_interactions = []
    
    # Index enhancers by chromosome for faster lookup
    enh_by_chr = defaultdict(list)
    for idx, enh in enhancers.iterrows():
        enh_by_chr[enh['chrom']].append({
            'start': enh['start'],
            'end': enh['end'],
            'name': enh['name'],
            'class': enh['class']
        })
    
    # Sort enhancers by start position
    for chrom in enh_by_chr:
        enh_by_chr[chrom].sort(key=lambda x: x['start'])
    
    def find_overlapping_enhancers(chrom, start, end):
        """Find enhancers overlapping with given region"""
        if chrom not in enh_by_chr:
            return []
        
        overlapping = []
        for enh in enh_by_chr[chrom]:
            # Check for overlap: not (end1 <= start2 or start1 >= end2)
            if not (end <= enh['start'] or start >= enh['end']):
                overlapping.append(enh)
        
        return overlapping
    
    print(f"Processing {len(interactions)} interactions...")
    progress_step = max(len(interactions) // 20, 1)
    
    for idx, interaction in interactions.iterrows():
        if idx % progress_step == 0:
            print(f"  Progress: {idx}/{len(interactions)} ({100*idx/len(interactions):.1f}%)")
        
        # Find overlapping enhancers for each anchor
        anchor1_enhancers = find_overlapping_enhancers(
            interaction['chr1'], interaction['start1'], interaction['end1']
        )
        anchor2_enhancers = find_overlapping_enhancers(
            interaction['chr2'], interaction['start2'], interaction['end2']
        )
        
        # Process E-E interactions (both anchors have enhancers)
        if anchor1_enhancers and anchor2_enhancers:
            for enh1 in anchor1_enhancers:
                for enh2 in anchor2_enhancers:
                    if enh1['class'] == enh2['class']:
                        contact_class = enh1['class']
                    else:
                        contact_class = 'cross_class'
                    
                    interaction_data = interaction.copy()
                    interaction_data['interaction_type'] = 'E-E'
                    interaction_data['contact_class'] = contact_class
                    interaction_data['enh1_name'] = enh1['name']
                    interaction_data['enh1_class'] = enh1['class']
                    interaction_data['enh2_name'] = enh2['name']
                    interaction_data['enh2_class'] = enh2['class']
                    
                    enh_interactions.append(interaction_data)
        
        # Process E-TSS interactions (only one anchor has enhancer)
        elif anchor1_enhancers or anchor2_enhancers:
            for enh in (anchor1_enhancers + anchor2_enhancers):
                interaction_data = interaction.copy()
                interaction_data['interaction_type'] = 'E-TSS'
                interaction_data['contact_class'] = enh['class']
                interaction_data['enh_name'] = enh['name']
                interaction_data['enh_class'] = enh['class']
                enh_interactions.append(interaction_data)
    
    if enh_interactions:
        enh_df = pd.DataFrame(enh_interactions)
        enh_df = enh_df.drop_duplicates()
        
        print(f"\nFound {len(enh_df)} enhancer interactions:")
        print(f"  E-E interactions: {sum(enh_df['interaction_type'] == 'E-E')}")
        print(f"  E-TSS interactions: {sum(enh_df['interaction_type'] == 'E-TSS')}")
        print("  By class:")
        print(enh_df['contact_class'].value_counts())
        
        return enh_df
    else:
        print("No enhancer interactions found")
        return pd.DataFrame()


def perform_statistical_comparison(enh_interactions, null_model, reference_condition='DOX', fdr_threshold=0.05):
    """
    Compare enhancer interactions between conditions using the null model.
    
    Key concept: The null model represents interactions under the null hypothesis
    (no real biological difference). We compare observed logFC values against this
    null distribution to determine statistical significance.
    """
    print("\n" + "="*70)
    print("STATISTICAL COMPARISON WITH NULL MODEL")
    print("="*70)
    
    if enh_interactions.empty:
        print("No enhancer interactions to analyze")
        return pd.DataFrame()
    
    # Load and characterize null model
    null_logfc_mean = 0
    null_logfc_std = 1
    
    if null_model is not None and not null_model.empty:
        print("\nNull Model Characterization:")
        null_logfc_values = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
        
        if len(null_logfc_values) > 0:
            null_logfc_mean = null_logfc_values.mean()
            null_logfc_std = null_logfc_values.std()
            
            print(f"  Null model size: {len(null_logfc_values)} interactions")
            print(f"  Mean logFC: {null_logfc_mean:.4f}")
            print(f"  Std logFC: {null_logfc_std:.4f}")
            print(f"  Range: [{null_logfc_values.min():.4f}, {null_logfc_values.max():.4f}]")
            
            # Check if null model looks correct (should be centered near 0)
            if abs(null_logfc_mean) > 0.5:
                print(f"  WARNING: Null model mean ({null_logfc_mean:.4f}) is far from 0!")
                print(f"  This may indicate an issue with the null model generation.")
        else:
            print("  WARNING: No valid logFC values in null model!")
    else:
        print("\nWARNING: No null model provided!")
        print("Using default null distribution (mean=0, std=1)")
        print("Results may not be properly calibrated.")
    
    # Get available conditions
    available_conditions = enh_interactions['condition_clean'].unique()
    print(f"\nAvailable conditions: {list(available_conditions)}")
    
    if reference_condition not in available_conditions:
        print(f"WARNING: Reference condition '{reference_condition}' not found!")
        reference_condition = available_conditions[0]
        print(f"Using '{reference_condition}' as reference")
    
    # Perform comparisons
    comparison_results = []
    infected_conditions = [c for c in ['wMel', 'wRi', 'wWil'] if c in available_conditions]
    
    print(f"\nPerforming comparisons:")
    print(f"  Reference: {reference_condition}")
    print(f"  Test conditions: {infected_conditions}")
    
    for condition in infected_conditions:
        print(f"\n{'='*70}")
        print(f"Comparing: {condition} vs {reference_condition}")
        print(f"{'='*70}")
        
        # Get data for this comparison
        ref_data = enh_interactions[enh_interactions['condition_clean'] == reference_condition]
        cond_data = enh_interactions[enh_interactions['condition_clean'] == condition]
        
        print(f"Reference interactions: {len(ref_data)}")
        print(f"Test interactions: {len(cond_data)}")
        
        if ref_data.empty or cond_data.empty:
            print("Skipping - insufficient data")
            continue
        
        # Analyze by interaction type and enhancer class
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
                
                # Get logFC values
                ref_logfc = pd.to_numeric(ref_class['logFC'], errors='coerce').dropna()
                cond_logfc = pd.to_numeric(cond_class['logFC'], errors='coerce').dropna()
                
                if len(ref_logfc) < 3 or len(cond_logfc) < 3:
                    continue
                
                print(f"\n  {int_type} - {contact_class}:")
                print(f"    Ref: n={len(ref_logfc)}, mean={ref_logfc.mean():.3f}")
                print(f"    Test: n={len(cond_logfc)}, mean={cond_logfc.mean():.3f}")
                
                try:
                    # Calculate observed difference
                    mean_diff = cond_logfc.mean() - ref_logfc.mean()
                    print(f"    Observed difference: {mean_diff:.3f}")
                    
                    # Method 1: Z-test against null distribution
                    if null_logfc_std > 0:
                        se_diff = np.sqrt(null_logfc_std**2/len(ref_logfc) + 
                                         null_logfc_std**2/len(cond_logfc))
                        z_score = (mean_diff - null_logfc_mean) / se_diff
                        p_value_z = 2 * (1 - stats.norm.cdf(abs(z_score)))
                        print(f"    Z-test: z={z_score:.3f}, p={p_value_z:.4e}")
                    else:
                        p_value_z = 1.0
                        z_score = 0
                    
                    # Method 2: Permutation test
                    n_permutations = 1000
                    null_diffs = []
                    
                    for _ in range(n_permutations):
                        perm_ref = np.random.normal(null_logfc_mean, null_logfc_std, len(ref_logfc))
                        perm_cond = np.random.normal(null_logfc_mean, null_logfc_std, len(cond_logfc))
                        null_diff = perm_cond.mean() - perm_ref.mean()
                        null_diffs.append(null_diff)
                    
                    null_diffs = np.array(null_diffs)
                    p_value_perm = np.sum(np.abs(null_diffs) >= np.abs(mean_diff)) / n_permutations
                    print(f"    Permutation: p={p_value_perm:.4f}")
                    
                    # Method 3: Bootstrap confidence interval
                    n_bootstrap = 1000
                    bootstrap_diffs = []
                    
                    for _ in range(n_bootstrap):
                        boot_ref = np.random.choice(ref_logfc, len(ref_logfc), replace=True)
                        boot_cond = np.random.choice(cond_logfc, len(cond_logfc), replace=True)
                        bootstrap_diffs.append(boot_cond.mean() - boot_ref.mean())
                    
                    ci_lower = np.percentile(bootstrap_diffs, 2.5)
                    ci_upper = np.percentile(bootstrap_diffs, 97.5)
                    ci_excludes_zero = not (ci_lower <= 0 <= ci_upper)
                    print(f"    95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
                    
                    # Effect size
                    effect_size = abs(mean_diff - null_logfc_mean) / null_logfc_std if null_logfc_std > 0 else 0
                    print(f"    Effect size: {effect_size:.3f} std")
                    
                    # Biological significance (>0.5 log2FC = ~1.4-fold change)
                    biologically_significant = abs(mean_diff) > 0.5
                    
                    comparison_results.append({
                        'comparison': f'{condition}_vs_{reference_condition}',
                        'interaction_type': int_type,
                        'contact_class': contact_class,
                        'n_ref': len(ref_logfc),
                        'n_test': len(cond_logfc),
                        'mean_logfc_ref': ref_logfc.mean(),
                        'mean_logfc_test': cond_logfc.mean(),
                        'logfc_difference': mean_diff,
                        'null_expectation': null_logfc_mean,
                        'deviation_from_null': mean_diff - null_logfc_mean,
                        'effect_size': effect_size,
                        'z_score': z_score,
                        'p_value_z_test': p_value_z,
                        'p_value_permutation': p_value_perm,
                        'ci_lower': ci_lower,
                        'ci_upper': ci_upper,
                        'ci_excludes_zero': ci_excludes_zero,
                        'biologically_significant': biologically_significant
                    })
                    
                except Exception as e:
                    print(f"    Error in analysis: {e}")
                    continue
    
    if not comparison_results:
        print("\nNo comparisons could be performed!")
        return pd.DataFrame()
    
    # Convert to DataFrame and apply FDR correction
    results_df = pd.DataFrame(comparison_results)
    
    print(f"\n{'='*70}")
    print(f"FDR CORRECTION")
    print(f"{'='*70}")
    
    for p_col in ['p_value_z_test', 'p_value_permutation']:
        if p_col in results_df.columns:
            _, fdr_values, _, _ = multipletests(results_df[p_col], method='fdr_bh')
            results_df[f'{p_col}_fdr'] = fdr_values
            results_df[f'{p_col}_significant'] = fdr_values < fdr_threshold
            
            n_sig = sum(fdr_values < fdr_threshold)
            print(f"{p_col}: {n_sig}/{len(results_df)} significant (FDR < {fdr_threshold})")
    
    # Combined significance
    sig_cols = [col for col in results_df.columns if col.endswith('_significant')]
    if sig_cols:
        results_df['any_significant'] = results_df[sig_cols].any(axis=1)
        print(f"\nOverall: {sum(results_df['any_significant'])}/{len(results_df)} significant by at least one method")
    
    # Add interpretation
    results_df['interpretation'] = results_df.apply(interpret_result, axis=1)
    
    return results_df


def interpret_result(row):
    """Interpret the statistical and biological significance."""
    interpretations = []
    
    if row.get('any_significant', False):
        interpretations.append("statistically_significant")
    
    if row.get('biologically_significant', False):
        interpretations.append("biologically_significant")
    
    effect_size = row.get('effect_size', 0)
    if effect_size > 1:
        interpretations.append("large_effect")
    elif effect_size > 0.5:
        interpretations.append("medium_effect")
    elif effect_size > 0.2:
        interpretations.append("small_effect")
    else:
        interpretations.append("negligible_effect")
    
    diff = row.get('logfc_difference', 0)
    if diff > 0:
        interpretations.append("increased_in_test")
    elif diff < 0:
        interpretations.append("decreased_in_test")
    else:
        interpretations.append("no_change")
    
    return "; ".join(interpretations)


def analyze_distance_dependence(enh_interactions):
    """Analyze how enhancer interactions change with genomic distance."""
    print("\nAnalyzing distance dependence...")
    
    if enh_interactions.empty:
        return pd.DataFrame()
    
    # Filter for cis interactions
    cis_interactions = enh_interactions[enh_interactions['chr1'] == enh_interactions['chr2']].copy()
    
    if cis_interactions.empty:
        return pd.DataFrame()
    
    # Calculate distance
    cis_interactions['distance'] = abs(
        pd.to_numeric(cis_interactions['start2'], errors='coerce') - 
        pd.to_numeric(cis_interactions['start1'], errors='coerce')
    )
    
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
    
    # Analyze by condition, type, and distance
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
                    
                    if len(bin_data) < 3:
                        continue
                    
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
    print(f"Distance analysis completed: {len(distance_df)} combinations")
    
    return distance_df


def create_comprehensive_plots(enh_interactions, comparison_results, null_model, output_prefix):
    """Create comprehensive visualization of results."""
    
    if comparison_results.empty:
        print("No comparison results to plot")
        return
    
    print("\nCreating plots...")
    
    # Create main figure
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    
    # Plot 1: Effect sizes
    ax1 = fig.add_subplot(gs[0, 0])
    scatter_data = comparison_results[comparison_results['any_significant'] == True]
    if not scatter_data.empty:
        ax1.scatter(scatter_data['logfc_difference'], scatter_data['effect_size'],
                   c='red', alpha=0.6, s=50, label='Significant')
    
    scatter_data = comparison_results[comparison_results['any_significant'] == False]
    if not scatter_data.empty:
        ax1.scatter(scatter_data['logfc_difference'], scatter_data['effect_size'],
                   c='gray', alpha=0.3, s=30, label='Not significant')
    
    ax1.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    ax1.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5, label='Medium effect')
    ax1.axhline(y=1.0, color='red', linestyle=':', alpha=0.5, label='Large effect')
    ax1.set_xlabel('LogFC Difference')
    ax1.set_ylabel('Effect Size (std)')
    ax1.set_title('Effect Sizes vs Null Model')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: P-value comparison
    ax2 = fig.add_subplot(gs[0, 1])
    if 'p_value_z_test' in comparison_results.columns and 'p_value_permutation' in comparison_results.columns:
        ax2.scatter(comparison_results['p_value_z_test'],
                   comparison_results['p_value_permutation'],
                   alpha=0.6, s=40)
        ax2.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        ax2.axhline(y=0.05, color='red', linestyle=':', alpha=0.5)
        ax2.axvline(x=0.05, color='red', linestyle=':', alpha=0.5)
        ax2.set_xlabel('Z-test P-value')
        ax2.set_ylabel('Permutation P-value')
        ax2.set_title('P-value Method Comparison')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.grid(True, alpha=0.3)
    
    # Plot 3: Confidence intervals
    ax3 = fig.add_subplot(gs[0, 2])
    y_pos = range(len(comparison_results))
    colors = ['red' if sig else 'gray' for sig in comparison_results.get('any_significant', [False]*len(comparison_results))]
    ax3.errorbar(comparison_results['logfc_difference'], y_pos,
               xerr=[comparison_results['logfc_difference'] - comparison_results['ci_lower'],
                     comparison_results['ci_upper'] - comparison_results['logfc_difference']],
               fmt='o', alpha=0.6, c='gray', ecolor=colors, markersize=4)
    ax3.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    ax3.set_xlabel('LogFC Difference')
    ax3.set_ylabel('Comparison Index')
    ax3.set_title('95% Confidence Intervals')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Deviation from null
    ax4 = fig.add_subplot(gs[0, 3])
    ax4.hist(comparison_results['deviation_from_null'], bins=20, alpha=0.7, edgecolor='black', color='steelblue')
    ax4.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Null expectation')
    ax4.set_xlabel('Deviation from Null')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Distribution of Deviations')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Effect sizes by interaction type
    ax5 = fig.add_subplot(gs[1, 0])
    if 'interaction_type' in comparison_results.columns:
        for int_type in comparison_results['interaction_type'].unique():
            type_data = comparison_results[comparison_results['interaction_type'] == int_type]
            ax5.scatter([int_type] * len(type_data), type_data['effect_size'], 
                       alpha=0.6, s=50, label=int_type)
        ax5.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5)
        ax5.axhline(y=1.0, color='red', linestyle=':', alpha=0.5)
        ax5.set_ylabel('Effect Size')
        ax5.set_title('Effect Sizes by Interaction Type')
        ax5.grid(True, alpha=0.3, axis='y')
    
    # Plot 6: Effect sizes by contact class
    ax6 = fig.add_subplot(gs[1, 1])
    if 'contact_class' in comparison_results.columns:
        for contact_class in comparison_results['contact_class'].unique():
            class_data = comparison_results[comparison_results['contact_class'] == contact_class]
            ax6.scatter([contact_class] * len(class_data), class_data['effect_size'],
                       alpha=0.6, s=50, label=contact_class)
        ax6.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5)
        ax6.axhline(y=1.0, color='red', linestyle=':', alpha=0.5)
        ax6.set_ylabel('Effect Size')
        ax6.set_title('Effect Sizes by Contact Class')
        ax6.grid(True, alpha=0.3, axis='y')
    
    # Plot 7: Null model distribution
    ax7 = fig.add_subplot(gs[1, 2])
    if null_model is not None and not null_model.empty:
        null_logfc = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
        if len(null_logfc) > 0:
            ax7.hist(null_logfc, bins=30, alpha=0.7, color='lightblue', 
                   edgecolor='black', label='Null distribution')
            ax7.axvline(x=null_logfc.mean(), color='red', linestyle='--', 
                      linewidth=2, label=f'Mean: {null_logfc.mean():.3f}')
            ax7.axvline(x=0, color='black', linestyle=':', linewidth=1, alpha=0.5)
            ax7.set_xlabel('LogFC')
            ax7.set_ylabel('Frequency')
            ax7.set_title('Null Model Distribution')
            ax7.legend()
            ax7.grid(True, alpha=0.3)
    
    # Plot 8: Significance summary
    ax8 = fig.add_subplot(gs[1, 3])
    if 'comparison' in comparison_results.columns and 'any_significant' in comparison_results.columns:
        sig_counts = comparison_results.groupby('comparison')['any_significant'].sum()
        total_counts = comparison_results.groupby('comparison').size()
        
        x = range(len(sig_counts))
        ax8.bar(x, total_counts, alpha=0.3, label='Total', color='gray')
        ax8.bar(x, sig_counts, alpha=0.8, label='Significant', color='red')
        ax8.set_xticks(x)
        ax8.set_xticklabels(sig_counts.index, rotation=45, ha='right')
        ax8.set_ylabel('Count')
        ax8.set_title('Significant Comparisons')
        ax8.legend()
        ax8.grid(True, alpha=0.3, axis='y')
    
    # Plot 9: LogFC comparison by condition
    ax9 = fig.add_subplot(gs[2, :2])
    if 'comparison' in comparison_results.columns:
        comparison_results_sorted = comparison_results.sort_values('comparison')
        x_pos = range(len(comparison_results_sorted))
        colors_sig = ['red' if sig else 'gray' 
                     for sig in comparison_results_sorted.get('any_significant', [False]*len(comparison_results_sorted))]
        
        ax9.scatter(x_pos, comparison_results_sorted['logfc_difference'], 
                   c=colors_sig, alpha=0.6, s=60)
        ax9.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        ax9.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5, label='Bio. sig. threshold')
        ax9.axhline(y=-0.5, color='orange', linestyle=':', alpha=0.5)
        
        # Add labels for significant points
        for i, (idx, row) in enumerate(comparison_results_sorted.iterrows()):
            if row.get('any_significant', False):
                label = f"{row['interaction_type']}-{row['contact_class']}"
                ax9.annotate(label, (i, row['logfc_difference']), 
                           fontsize=6, rotation=45, ha='right')
        
        ax9.set_xlabel('Comparison Index')
        ax9.set_ylabel('LogFC Difference')
        ax9.set_title('LogFC Differences Across All Comparisons')
        ax9.legend()
        ax9.grid(True, alpha=0.3)
    
    # Plot 10: Effect size vs biological significance
    ax10 = fig.add_subplot(gs[2, 2:])
    if 'biologically_significant' in comparison_results.columns:
        bio_sig = comparison_results[comparison_results['biologically_significant'] == True]
        bio_not_sig = comparison_results[comparison_results['biologically_significant'] == False]
        
        if not bio_sig.empty:
            ax10.scatter(bio_sig['effect_size'], bio_sig['logfc_difference'],
                       c='red', alpha=0.6, s=60, label='Biologically significant')
        if not bio_not_sig.empty:
            ax10.scatter(bio_not_sig['effect_size'], bio_not_sig['logfc_difference'],
                       c='gray', alpha=0.3, s=40, label='Not biologically significant')
        
        ax10.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5)
        ax10.axhline(y=-0.5, color='orange', linestyle=':', alpha=0.5)
        ax10.axvline(x=0.5, color='blue', linestyle=':', alpha=0.5)
        ax10.set_xlabel('Effect Size (std)')
        ax10.set_ylabel('LogFC Difference')
        ax10.set_title('Statistical vs Biological Significance')
        ax10.legend()
        ax10.grid(True, alpha=0.3)
    
    plt.savefig(f"{output_prefix}_comprehensive_analysis.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Comprehensive plots saved to {output_prefix}_comprehensive_analysis.pdf")


def create_summary_report(enh_interactions, comparison_results, null_model, output_prefix):
    """Create a detailed text summary report."""
    
    with open(f"{output_prefix}_analysis_report.txt", 'w') as f:
        f.write("="*80 + "\n")
        f.write("ENHANCER CLASS ANALYSIS REPORT\n")
        f.write("="*80 + "\n\n")
        
        # Null model summary
        f.write("NULL MODEL SUMMARY\n")
        f.write("-"*80 + "\n")
        if null_model is not None and not null_model.empty:
            null_logfc = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
            if len(null_logfc) > 0:
                f.write(f"Null model size: {len(null_logfc)} interactions\n")
                f.write(f"Mean logFC: {null_logfc.mean():.4f}\n")
                f.write(f"Std logFC: {null_logfc.std():.4f}\n")
                f.write(f"Range: [{null_logfc.min():.4f}, {null_logfc.max():.4f}]\n")
                f.write(f"Median: {null_logfc.median():.4f}\n")
        else:
            f.write("No null model provided\n")
        f.write("\n")
        
        # Enhancer interaction summary
        f.write("ENHANCER INTERACTIONS SUMMARY\n")
        f.write("-"*80 + "\n")
        if not enh_interactions.empty:
            f.write(f"Total enhancer interactions: {len(enh_interactions)}\n")
            f.write(f"E-E interactions: {sum(enh_interactions['interaction_type'] == 'E-E')}\n")
            f.write(f"E-TSS interactions: {sum(enh_interactions['interaction_type'] == 'E-TSS')}\n\n")
            
            f.write("By condition:\n")
            for cond in enh_interactions['condition_clean'].unique():
                n = sum(enh_interactions['condition_clean'] == cond)
                f.write(f"  {cond}: {n} interactions\n")
            
            f.write("\nBy contact class:\n")
            for cls in enh_interactions['contact_class'].unique():
                n = sum(enh_interactions['contact_class'] == cls)
                f.write(f"  {cls}: {n} interactions\n")
        else:
            f.write("No enhancer interactions found\n")
        f.write("\n")
        
        # Comparison results summary
        f.write("STATISTICAL COMPARISONS SUMMARY\n")
        f.write("-"*80 + "\n")
        if not comparison_results.empty:
            f.write(f"Total comparisons performed: {len(comparison_results)}\n\n")
            
            if 'any_significant' in comparison_results.columns:
                n_sig = sum(comparison_results['any_significant'])
                f.write(f"Significant comparisons (any method): {n_sig}/{len(comparison_results)}\n")
            
            if 'p_value_z_test_significant' in comparison_results.columns:
                n_z = sum(comparison_results['p_value_z_test_significant'])
                f.write(f"Significant by Z-test: {n_z}/{len(comparison_results)}\n")
            
            if 'p_value_permutation_significant' in comparison_results.columns:
                n_perm = sum(comparison_results['p_value_permutation_significant'])
                f.write(f"Significant by permutation: {n_perm}/{len(comparison_results)}\n")
            
            if 'biologically_significant' in comparison_results.columns:
                n_bio = sum(comparison_results['biologically_significant'])
                f.write(f"Biologically significant (>0.5 log2FC): {n_bio}/{len(comparison_results)}\n")
            
            f.write("\nEffect size distribution:\n")
            if 'effect_size' in comparison_results.columns:
                large = sum(comparison_results['effect_size'] > 1)
                medium = sum((comparison_results['effect_size'] > 0.5) & (comparison_results['effect_size'] <= 1))
                small = sum((comparison_results['effect_size'] > 0.2) & (comparison_results['effect_size'] <= 0.5))
                f.write(f"  Large effects (>1 std): {large}\n")
                f.write(f"  Medium effects (0.5-1 std): {medium}\n")
                f.write(f"  Small effects (0.2-0.5 std): {small}\n")
            
            # Detailed results for significant comparisons
            if 'any_significant' in comparison_results.columns:
                sig_results = comparison_results[comparison_results['any_significant']]
                if not sig_results.empty:
                    f.write("\n" + "="*80 + "\n")
                    f.write("SIGNIFICANT COMPARISONS DETAILS\n")
                    f.write("="*80 + "\n\n")
                    
                    for idx, row in sig_results.iterrows():
                        f.write(f"\nComparison: {row['comparison']}\n")
                        f.write(f"  Interaction type: {row['interaction_type']}\n")
                        f.write(f"  Contact class: {row['contact_class']}\n")
                        f.write(f"  Sample sizes: n_ref={row['n_ref']}, n_test={row['n_test']}\n")
                        f.write(f"  Mean logFC (ref): {row['mean_logfc_ref']:.3f}\n")
                        f.write(f"  Mean logFC (test): {row['mean_logfc_test']:.3f}\n")
                        f.write(f"  LogFC difference: {row['logfc_difference']:.3f}\n")
                        f.write(f"  Effect size: {row['effect_size']:.3f} std\n")
                        f.write(f"  P-value (Z-test): {row['p_value_z_test']:.4e}\n")
                        f.write(f"  P-value (permutation): {row['p_value_permutation']:.4f}\n")
                        f.write(f"  95% CI: [{row['ci_lower']:.3f}, {row['ci_upper']:.3f}]\n")
                        f.write(f"  Interpretation: {row['interpretation']}\n")
                        f.write("-"*80 + "\n")
        else:
            f.write("No comparison results\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"Analysis report saved to {output_prefix}_analysis_report.txt")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze enhancer class interactions with proper null model comparison',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single combined file with condition column
  python script.py --enhancers enhancers.bed --interactions all_interactions.csv \\
                   --null_model null_results.csv --output_prefix results/analysis

  # Multiple condition-specific files
  python script.py --enhancers enhancers.bed \\
                   --interactions "results/*_vs_DOX_interactions.csv" \\
                   --null_model null_results.csv --output_prefix results/analysis
        """
    )
    
    parser.add_argument('--enhancers', required=True, 
                       help='Enhancer BED file with classifications')
    parser.add_argument('--interactions', required=True,
                       help='Interaction CSV file or pattern (e.g., "*_vs_DOX_interactions.csv")')
    parser.add_argument('--null_model', required=True,
                       help='Null model CSV file from diffHic')
    parser.add_argument('--classification', 
                       help='Optional separate enhancer classification file')
    parser.add_argument('--tss', 
                       help='TSS BED file (optional)')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significance (default: 0.05)')
    parser.add_argument('--reference_condition', default='DOX',
                       help='Reference condition name (default: DOX)')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    print("\n" + "="*80)
    print("ENHANCER CLASS ANALYSIS")
    print("="*80 + "\n")
    
    # Load null model
    print("Loading null model...")
    null_model = pd.read_csv(args.null_model)
    print(f"Loaded null model with {len(null_model)} entries\n")
    
    # Load interactions
    interactions, reference_condition = load_interactions_by_condition(
        args.interactions, args.reference_condition
    )
    
    # Load and classify enhancers
    enhancers = classify_enhancers(args.enhancers, args.classification)
    
    if enhancers.empty:
        print("ERROR: No valid enhancers loaded!")
        sys.exit(1)
    
    # Find enhancer interactions
    enh_interactions = find_enhancer_interactions(interactions, enhancers, args.tss)
    
    if enh_interactions.empty:
        print("WARNING: No enhancer interactions found!")
    else:
        # Perform statistical comparison
        comparison_results = perform_statistical_comparison(
            enh_interactions, null_model, reference_condition, args.fdr_threshold
        )
        
        # Analyze distance dependence
        distance_analysis = analyze_distance_dependence(enh_interactions)
        
        # Create visualizations
        if not comparison_results.empty:
            create_comprehensive_plots(enh_interactions, comparison_results, null_model, args.output_prefix)
            create_summary_report(enh_interactions, comparison_results, null_model, args.output_prefix)
        
        # Save results
        print("\nSaving results...")
        
        enh_interactions.to_csv(f"{args.output_prefix}_enhancer_interactions.tsv", 
                               sep='\t', index=False)
        print(f"  Saved: {args.output_prefix}_enhancer_interactions.tsv")
        
        if not comparison_results.empty:
            comparison_results.to_csv(f"{args.output_prefix}_statistical_comparisons.tsv",
                                    sep='\t', index=False)
            print(f"  Saved: {args.output_prefix}_statistical_comparisons.tsv")
            
            # Save significant results
            if 'any_significant' in comparison_results.columns:
                sig_results = comparison_results[comparison_results['any_significant']]
                if not sig_results.empty:
                    sig_results.to_csv(f"{args.output_prefix}_significant_results.tsv",
                                     sep='\t', index=False)
                    print(f"  Saved: {args.output_prefix}_significant_results.tsv")
        
        if not distance_analysis.empty:
            distance_analysis.to_csv(f"{args.output_prefix}_distance_analysis.tsv",
                                   sep='\t', index=False)
            print(f"  Saved: {args.output_prefix}_distance_analysis.tsv")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)


if __name__ == '__main__':
    main()
