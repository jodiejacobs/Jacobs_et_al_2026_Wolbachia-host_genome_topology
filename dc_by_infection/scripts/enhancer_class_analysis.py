#!/usr/bin/env python3
"""
This is a working version as of 06/19/25
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
    
    # Your BED file has 7 columns: chrom start end name score strand class
    try:
        enhancers = pd.read_csv(enhancer_file, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'class'])
        print("Loaded enhancer file with 7 columns (including class)")
    except:
        try:
            # Fallback: try with header
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
        # Classification is already in the BED file
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
        
        # Fill any missing with developmental
        enhancers['class'] = enhancers['class'].fillna('developmental')
        
    elif classification_file and Path(classification_file).exists():
        print(f"Loading enhancer classification from {classification_file}")
        try:
            classification = pd.read_csv(classification_file, sep='\t', header=None, 
                                       names=['name', 'class'])
            print(f"Classification file head:")
            print(classification.head())
            
            # Clean classification data
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
            
            # Merge with enhancers (this will override BED file classification if both exist)
            enhancers = enhancers.drop(columns=['class'], errors='ignore')  # Remove BED class first
            enhancers = enhancers.merge(classification, on='name', how='left')
            merged_count = sum(enhancers['class'].notna())
            print(f"Successfully merged {merged_count} classifications out of {len(enhancers)} enhancers")
            
        except Exception as e:
            print(f"Error loading classification file: {e}")
            # Fallback to heuristic if both BED and file fail
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
    Find enhancer-enhancer and enhancer-TSS interactions.
    """
    print("Finding enhancer interactions...")
    
    # Validate and clean data first
    print("Validating interaction coordinates...")
    
    # Clean interactions data
    interactions_clean = interactions.copy()
    
    # Convert coordinate columns to numeric and handle any issues
    coord_cols = ['start1', 'end1', 'start2', 'end2']
    for col in coord_cols:
        interactions_clean[col] = pd.to_numeric(interactions_clean[col], errors='coerce')
    
    # Remove rows with invalid coordinates
    interactions_clean = interactions_clean.dropna(subset=coord_cols)
    
    # Ensure coordinates are integers
    for col in coord_cols:
        interactions_clean[col] = interactions_clean[col].astype(int)
    
    # Ensure start <= end for each region
    interactions_clean.loc[interactions_clean['start1'] > interactions_clean['end1'], ['start1', 'end1']] = \
        interactions_clean.loc[interactions_clean['start1'] > interactions_clean['end1'], ['end1', 'start1']].values
    
    interactions_clean.loc[interactions_clean['start2'] > interactions_clean['end2'], ['start2', 'end2']] = \
        interactions_clean.loc[interactions_clean['start2'] > interactions_clean['end2'], ['end2', 'start2']].values
    
    print(f"Cleaned interactions: {len(interactions)} -> {len(interactions_clean)}")
    
    # Clean enhancers data
    enhancers_clean = enhancers.copy()
    enhancers_clean['start'] = pd.to_numeric(enhancers_clean['start'], errors='coerce')
    enhancers_clean['end'] = pd.to_numeric(enhancers_clean['end'], errors='coerce')
    enhancers_clean = enhancers_clean.dropna(subset=['start', 'end'])
    enhancers_clean['start'] = enhancers_clean['start'].astype(int)
    enhancers_clean['end'] = enhancers_clean['end'].astype(int)
    
    # Ensure start <= end
    enhancers_clean.loc[enhancers_clean['start'] > enhancers_clean['end'], ['start', 'end']] = \
        enhancers_clean.loc[enhancers_clean['start'] > enhancers_clean['end'], ['end', 'start']].values
    
    print(f"Cleaned enhancers: {len(enhancers)} -> {len(enhancers_clean)}")
    
    # Convert interactions to BedTool format with proper validation
    interactions_bed1 = interactions_clean[['chr1', 'start1', 'end1']].copy()
    interactions_bed1.columns = ['chrom', 'start', 'end']
    interactions_bed1['interaction_idx'] = interactions_clean.index
    
    interactions_bed2 = interactions_clean[['chr2', 'start2', 'end2']].copy() 
    interactions_bed2.columns = ['chrom', 'start', 'end']
    interactions_bed2['interaction_idx'] = interactions_clean.index
    
    # Validate BED format before creating BedTool objects
    def validate_bed_format(df, name):
        """Validate BED format data"""
        print(f"Validating {name} BED format...")
        
        # Check for required columns
        if not all(col in df.columns for col in ['chrom', 'start', 'end']):
            raise ValueError(f"Missing required columns in {name}")
        
        # Check data types
        if not pd.api.types.is_integer_dtype(df['start']) or not pd.api.types.is_integer_dtype(df['end']):
            print(f"Warning: Non-integer coordinates in {name}")
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
        
        # Check for negative coordinates
        if (df['start'] < 0).any() or (df['end'] < 0).any():
            print(f"Warning: Negative coordinates found in {name}")
            df = df[(df['start'] >= 0) & (df['end'] >= 0)]
        
        # Check for valid ranges
        if (df['start'] >= df['end']).any():
            print(f"Warning: Invalid ranges (start >= end) in {name}")
            df = df[df['start'] < df['end']]
        
        print(f"  {name} validated: {len(df)} regions")
        return df
    
    # Validate all datasets
    interactions_bed1 = validate_bed_format(interactions_bed1, "interactions_bed1")
    interactions_bed2 = validate_bed_format(interactions_bed2, "interactions_bed2")
    
    enhancers_bed = enhancers_clean[['chrom', 'start', 'end', 'name', 'class']].copy()
    enhancers_bed = validate_bed_format(enhancers_bed, "enhancers")
    
    if len(interactions_bed1) == 0 or len(interactions_bed2) == 0 or len(enhancers_bed) == 0:
        print("No valid data after cleaning - cannot proceed with intersection")
        return pd.DataFrame()
    
    try:
        # Convert to BedTools with error handling
        print("Creating BedTool objects...")
        enh_bt = pybedtools.BedTool.from_dataframe(enhancers_bed)
        int1_bt = pybedtools.BedTool.from_dataframe(interactions_bed1)
        int2_bt = pybedtools.BedTool.from_dataframe(interactions_bed2)
        
        print("Performing intersections...")
        
    except Exception as e:
        print(f"Error creating BedTool objects: {e}")
        # Fallback to manual overlap detection
        return find_overlaps_manual(interactions_clean, enhancers_clean)
        
    try:
        # Find overlaps
        enh_int1 = int1_bt.intersect(enh_bt, wa=True, wb=True)
        enh_int2 = int2_bt.intersect(enh_bt, wa=True, wb=True)
        
    except Exception as e:
        print(f"Error with bedtools intersect: {e}")
        print("Falling back to manual overlap detection...")
        return find_overlaps_manual(interactions_clean, enhancers_clean)
    
    # Process overlaps
    enh_interactions = []
    
    # Process anchor 1 overlaps
    anchor1_overlaps = {}
    try:
        for overlap in enh_int1:
            fields = str(overlap).strip().split('\t')
            if len(fields) >= 9:
                idx = int(fields[3])  # interaction_idx
                enh_name = fields[7]  # enhancer name
                enh_class = fields[8]  # enhancer class
                anchor1_overlaps[idx] = {'name': enh_name, 'class': enh_class}
    except Exception as e:
        print(f"Error processing anchor1 overlaps: {e}")
    
    # Process anchor 2 overlaps
    anchor2_overlaps = {}
    try:
        for overlap in enh_int2:
            fields = str(overlap).strip().split('\t')
            if len(fields) >= 9:
                idx = int(fields[3])  # interaction_idx
                enh_name = fields[7]  # enhancer name  
                enh_class = fields[8]  # enhancer class
                anchor2_overlaps[idx] = {'name': enh_name, 'class': enh_class}
    except Exception as e:
        print(f"Error processing anchor2 overlaps: {e}")
    
    # Combine overlaps to find enhancer interactions
    for idx in interactions_clean.index:
        if idx not in interactions_clean.index:
            continue
            
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
        try:
            interaction_data = interactions_clean.iloc[interactions_clean.index.get_loc(idx)].copy()
            interaction_data['interaction_type'] = interaction_type
            interaction_data['contact_class'] = contact_class
            
            if anchor1_enh:
                interaction_data['enh1_name'] = anchor1_enh['name']
                interaction_data['enh1_class'] = anchor1_enh['class']
            if anchor2_enh:
                interaction_data['enh2_name'] = anchor2_enh['name']
                interaction_data['enh2_class'] = anchor2_enh['class']
            
            enh_interactions.append(interaction_data)
        except Exception as e:
            print(f"Error processing interaction {idx}: {e}")
            continue
    
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

def find_overlaps_manual(interactions, enhancers):
    """
    Manual overlap detection as fallback when bedtools fails.
    """
    print("Using manual overlap detection...")
    
    enh_interactions = []
    
    # Create spatial indices for faster lookup
    from collections import defaultdict
    
    # Index enhancers by chromosome
    enh_by_chr = defaultdict(list)
    for idx, enh in enhancers.iterrows():
        enh_by_chr[enh['chrom']].append({
            'start': enh['start'],
            'end': enh['end'],
            'name': enh['name'],
            'class': enh['class']
        })
    
    # Sort enhancers by start position for faster searching
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
    
    for idx, interaction in interactions.iterrows():
        if idx % 100000 == 0:
            print(f"  Processed {idx} interactions...")
        
        # Find overlapping enhancers for each anchor
        anchor1_enhancers = find_overlapping_enhancers(
            interaction['chr1'], interaction['start1'], interaction['end1']
        )
        anchor2_enhancers = find_overlapping_enhancers(
            interaction['chr2'], interaction['start2'], interaction['end2']
        )
        
        # Process all combinations of overlapping enhancers
        for enh1 in anchor1_enhancers:
            for enh2 in anchor2_enhancers:
                # E-E interaction
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
        
        # Single enhancer interactions (E-TSS)
        for enh in anchor1_enhancers:
            if not anchor2_enhancers:  # Only anchor1 has enhancer
                interaction_data = interaction.copy()
                interaction_data['interaction_type'] = 'E-TSS'
                interaction_data['contact_class'] = enh['class']
                interaction_data['enh1_name'] = enh['name']
                interaction_data['enh1_class'] = enh['class']
                enh_interactions.append(interaction_data)
        
        for enh in anchor2_enhancers:
            if not anchor1_enhancers:  # Only anchor2 has enhancer
                interaction_data = interaction.copy()
                interaction_data['interaction_type'] = 'E-TSS'
                interaction_data['contact_class'] = enh['class']
                interaction_data['enh2_name'] = enh['name']
                interaction_data['enh2_class'] = enh['class']
                enh_interactions.append(interaction_data)
    
    if enh_interactions:
        enh_df = pd.DataFrame(enh_interactions)
        # Remove duplicates
        enh_df = enh_df.drop_duplicates()
        
        print(f"Found {len(enh_df)} enhancer interactions:")
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
    Compare enhancer interactions between conditions using diffHic null model as baseline.
    
    The null model represents the expected background distribution of logFC values
    when there are no real biological differences.
    """
    print("Performing statistical comparison with proper null model handling...")
    
    if enh_interactions.empty:
        print("No enhancer interactions to analyze")
        return pd.DataFrame()
    
    # Load null model parameters
    null_logfc_mean = 0  # Expected under null hypothesis
    null_logfc_std = 1   # Default, should be estimated from null model
    
    if null_model is not None and not null_model.empty:
        # Extract null distribution parameters
        null_logfc_values = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
        if len(null_logfc_values) > 0:
            null_logfc_mean = null_logfc_values.mean()
            null_logfc_std = null_logfc_values.std()
            print(f"Null model: mean logFC = {null_logfc_mean:.3f}, std = {null_logfc_std:.3f}")
        
        # Also get p-value distribution under null (should be uniform)
        null_pvalues = pd.to_numeric(null_model['PValue'], errors='coerce').dropna()
        if len(null_pvalues) > 0:
            print(f"Null p-values: mean = {null_pvalues.mean():.3f}, expected = 0.5")
    
    available_conditions = enh_interactions['condition_clean'].unique()
    print(f"Available conditions: {available_conditions}")
    
    if reference_condition not in available_conditions:
        reference_condition = available_conditions[0]
        print(f"Using '{reference_condition}' as reference")
    
    comparison_results = []
    ref_data = enh_interactions[enh_interactions['condition_clean'] == reference_condition]
    
    if ref_data.empty:
        print(f"No reference data found for condition: {reference_condition}")
        return pd.DataFrame()
    
    # Compare each infected condition to reference
    for condition in ['wMel', 'wRi', 'wWil']:
        if condition not in available_conditions:
            continue
            
        cond_data = enh_interactions[enh_interactions['condition_clean'] == condition]
        if cond_data.empty:
            continue
        
        print(f"Comparing {condition} vs {reference_condition}")
        
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
                
                try:
                    # Get logFC values
                    ref_logfc = pd.to_numeric(ref_class['logFC'], errors='coerce').dropna()
                    cond_logfc = pd.to_numeric(cond_class['logFC'], errors='coerce').dropna()
                    
                    if len(ref_logfc) == 0 or len(cond_logfc) == 0:
                        continue
                    
                    # Calculate the difference in logFC between conditions
                    mean_diff = cond_logfc.mean() - ref_logfc.mean()
                    
                    # Test against null model expectation
                    # H0: The observed difference could arise from the null distribution
                    # H1: The observed difference is significantly different from null
                    
                    # Method 1: Z-test against null distribution
                    if null_logfc_std > 0:
                        # Standard error for difference of means
                        se_diff = np.sqrt(null_logfc_std**2/len(ref_logfc) + 
                                         null_logfc_std**2/len(cond_logfc))
                        z_score = (mean_diff - null_logfc_mean) / se_diff
                        p_value_z = 2 * (1 - stats.norm.cdf(abs(z_score)))
                    else:
                        p_value_z = 1.0
                        z_score = 0
                    
                    # Method 2: Permutation test with null model as baseline
                    n_permutations = 1000
                    null_diffs = []
                    
                    for _ in range(n_permutations):
                        # Sample from null distribution
                        perm_ref = np.random.normal(null_logfc_mean, null_logfc_std, len(ref_logfc))
                        perm_cond = np.random.normal(null_logfc_mean, null_logfc_std, len(cond_logfc))
                        null_diff = perm_cond.mean() - perm_ref.mean()
                        null_diffs.append(null_diff)
                    
                    # Calculate p-value from permutation
                    null_diffs = np.array(null_diffs)
                    p_value_perm = np.sum(np.abs(null_diffs) >= np.abs(mean_diff)) / n_permutations
                    
                    # Method 3: Bootstrap confidence interval
                    n_bootstrap = 1000
                    bootstrap_diffs = []
                    
                    for _ in range(n_bootstrap):
                        boot_ref = np.random.choice(ref_logfc, len(ref_logfc), replace=True)
                        boot_cond = np.random.choice(cond_logfc, len(cond_logfc), replace=True)
                        bootstrap_diffs.append(boot_cond.mean() - boot_ref.mean())
                    
                    # 95% confidence interval
                    ci_lower = np.percentile(bootstrap_diffs, 2.5)
                    ci_upper = np.percentile(bootstrap_diffs, 97.5)
                    
                    # Effect size: how many standard deviations from null expectation
                    effect_size = abs(mean_diff - null_logfc_mean) / null_logfc_std if null_logfc_std > 0 else 0
                    
                    # Biological significance: is the change meaningful?
                    # Use a threshold based on your domain knowledge (e.g., >0.5 log2FC)
                    biologically_significant = abs(mean_diff) > 0.5
                    
                    comparison_results.append({
                        'comparison': f'{reference_condition}_vs_{condition}',
                        'interaction_type': int_type,
                        'contact_class': contact_class,
                        'n_ref': len(ref_logfc),
                        'n_cond': len(cond_logfc),
                        'mean_logfc_ref': ref_logfc.mean(),
                        'mean_logfc_cond': cond_logfc.mean(),
                        'logfc_difference': mean_diff,
                        'null_expectation': null_logfc_mean,
                        'deviation_from_null': mean_diff - null_logfc_mean,
                        'effect_size': effect_size,
                        'z_score': z_score,
                        'p_value_z_test': p_value_z,
                        'p_value_permutation': p_value_perm,
                        'ci_lower': ci_lower,
                        'ci_upper': ci_upper,
                        'biologically_significant': biologically_significant
                    })
                    
                except Exception as e:
                    print(f"Error in comparison for {int_type} {contact_class}: {e}")
                    continue
    
    if not comparison_results:
        return pd.DataFrame()
    
    # Convert to DataFrame
    results_df = pd.DataFrame(comparison_results)
    
    # Apply FDR correction to each p-value method
    for p_col in ['p_value_z_test', 'p_value_permutation']:
        if p_col in results_df.columns:
            _, fdr_values, _, _ = multipletests(results_df[p_col], method='fdr_bh')
            results_df[f'{p_col}_fdr'] = fdr_values
            results_df[f'{p_col}_significant'] = fdr_values < fdr_threshold
    
    # Combined significance: significant by any method
    sig_cols = [col for col in results_df.columns if col.endswith('_significant')]
    if sig_cols:
        results_df['any_significant'] = results_df[sig_cols].any(axis=1)
    
    # Add interpretation
    results_df['interpretation'] = results_df.apply(interpret_result, axis=1)
    
    n_significant = results_df['any_significant'].sum() if 'any_significant' in results_df.columns else 0
    print(f"Found {n_significant} comparisons significant by at least one method")
    
    return results_df

def interpret_result(row):
    """
    Interpret the statistical and biological significance of results.
    """
    interpretations = []
    
    # Statistical significance
    if row.get('any_significant', False):
        interpretations.append("statistically_significant")
    
    # Biological significance
    if row.get('biologically_significant', False):
        interpretations.append("biologically_significant")
    
    # Effect size interpretation
    effect_size = row.get('effect_size', 0)
    if effect_size > 1:
        interpretations.append("large_effect")
    elif effect_size > 0.5:
        interpretations.append("medium_effect")
    elif effect_size > 0.2:
        interpretations.append("small_effect")
    else:
        interpretations.append("negligible_effect")
    
    # Direction of change
    diff = row.get('logfc_difference', 0)
    if diff > 0:
        interpretations.append("increased_in_infection")
    elif diff < 0:
        interpretations.append("decreased_in_infection")
    else:
        interpretations.append("no_change")
    
    return "; ".join(interpretations)


def create_summary_with_null(comparison_results, null_model):
    """
    Create a summary table including null model information.
    """
    summary = {}
    
    if not comparison_results.empty:
        # Basic statistics
        summary['total_comparisons'] = len(comparison_results)
        summary['statistically_significant'] = comparison_results.get('any_significant', pd.Series([False]*len(comparison_results))).sum()
        summary['biologically_significant'] = comparison_results.get('biologically_significant', pd.Series([False]*len(comparison_results))).sum()
        summary['both_significant'] = (
            comparison_results.get('any_significant', pd.Series([False]*len(comparison_results))) & 
            comparison_results.get('biologically_significant', pd.Series([False]*len(comparison_results)))
        ).sum()
        
        # Effect size distribution
        if 'effect_size' in comparison_results.columns:
            summary['mean_effect_size'] = comparison_results['effect_size'].mean()
            summary['large_effects'] = (comparison_results['effect_size'] > 1).sum()
            summary['medium_effects'] = ((comparison_results['effect_size'] > 0.5) & 
                                        (comparison_results['effect_size'] <= 1)).sum()
            summary['small_effects'] = ((comparison_results['effect_size'] > 0.2) & 
                                       (comparison_results['effect_size'] <= 0.5)).sum()
    
    # Null model statistics
    if null_model is not None and not null_model.empty:
        null_logfc = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
        if len(null_logfc) > 0:
            summary['null_model_mean'] = null_logfc.mean()
            summary['null_model_std'] = null_logfc.std()
            summary['null_model_size'] = len(null_logfc)
    
    return summary

def create_null_comparison_plots(comparison_results, null_model, output_prefix):
    """
    Create plots showing comparison against null model.
    """
    if comparison_results.empty:
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot 1: Effect sizes vs null expectation
    ax = axes[0, 0]
    ax.scatter(comparison_results['logfc_difference'], 
               comparison_results['effect_size'],
               c=['red' if sig else 'gray' for sig in comparison_results.get('any_significant', [False]*len(comparison_results))],
               alpha=0.6)
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.5, label='Null expectation')
    ax.set_xlabel('LogFC Difference (Infected - Uninfected)')
    ax.set_ylabel('Effect Size (std deviations from null)')
    ax.set_title('Effect Sizes vs Null Model')
    ax.legend()
    
    # Plot 2: P-value comparison
    ax = axes[0, 1]
    if 'p_value_z_test' in comparison_results.columns and 'p_value_permutation' in comparison_results.columns:
        ax.scatter(comparison_results['p_value_z_test'],
                   comparison_results['p_value_permutation'],
                   alpha=0.6)
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        ax.set_xlabel('Z-test P-value')
        ax.set_ylabel('Permutation P-value')
        ax.set_title('P-value Method Comparison')
    
    # Plot 3: Confidence intervals
    ax = axes[0, 2]
    if 'ci_lower' in comparison_results.columns and 'ci_upper' in comparison_results.columns:
        y_pos = range(len(comparison_results))
        ax.errorbar(comparison_results['logfc_difference'], y_pos,
                   xerr=[comparison_results['logfc_difference'] - comparison_results['ci_lower'],
                         comparison_results['ci_upper'] - comparison_results['logfc_difference']],
                   fmt='o', alpha=0.6)
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel('LogFC Difference')
        ax.set_ylabel('Comparison Index')
        ax.set_title('95% Confidence Intervals')
    
    # Plot 4: Distribution of deviations from null
    ax = axes[1, 0]
    ax.hist(comparison_results['deviation_from_null'], bins=20, alpha=0.7, edgecolor='black')
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Null expectation')
    ax.set_xlabel('Deviation from Null Expectation')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of Deviations from Null')
    ax.legend()
    
    # Plot 5: Effect sizes by interaction type
    ax = axes[1, 1]
    if 'interaction_type' in comparison_results.columns:
        int_types = comparison_results['interaction_type'].unique()
        for i, int_type in enumerate(int_types):
            type_data = comparison_results[comparison_results['interaction_type'] == int_type]
            ax.scatter([i] * len(type_data), type_data['effect_size'], alpha=0.6, label=int_type)
        ax.set_xticks(range(len(int_types)))
        ax.set_xticklabels(int_types)
        ax.set_ylabel('Effect Size')
        ax.set_title('Effect Sizes by Interaction Type')
        ax.legend()
    
    # Plot 6: Null model distribution (if available)
    ax = axes[1, 2]
    if null_model is not None and not null_model.empty:
        null_logfc = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
        if len(null_logfc) > 0:
            ax.hist(null_logfc, bins=20, alpha=0.7, color='lightblue', 
                   edgecolor='black', label='Null distribution')
            ax.axvline(x=null_logfc.mean(), color='red', linestyle='--', 
                      linewidth=2, label=f'Null mean: {null_logfc.mean():.3f}')
            ax.set_xlabel('LogFC')
            ax.set_ylabel('Frequency')
            ax.set_title('Null Model Distribution')
            ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_null_model_comparison.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Null model comparison plots saved to {output_prefix}_null_model_comparison.pdf")

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

# def main():
#     parser = argparse.ArgumentParser(description='Analyze enhancer class interactions')
#     parser.add_argument('--enhancers', required=True, help='Enhancer BED file')
#     parser.add_argument('--interactions', required=True, help='Differential interactions CSV file')
#     parser.add_argument('--null_model', help='Null model CSV file from diffHic')
#     parser.add_argument('--classification', help='Optional enhancer classification file')
#     parser.add_argument('--tss', help='TSS BED file')
#     parser.add_argument('--fdr_threshold', type=float, default=0.05, help='FDR threshold')
#     parser.add_argument('--reference_condition', default='DOX', help='Reference condition name')
#     parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    
#     args = parser.parse_args()
    
#     # Load null model if provided
#     null_model = None
#     if args.null_model:
#         null_model = pd.read_csv(args.null_model)
#         print(f"Loaded null model with {len(null_model)} entries")
    
#     # Load and process interactions
#     interactions, reference_condition = load_and_process_interactions(
#         args.interactions, args.reference_condition
#     )
    
#     # Load and classify enhancers
#     enhancers = classify_enhancers(args.enhancers, args.classification)
    
#     # Find enhancer interactions
#     enh_interactions = find_enhancer_interactions(interactions, enhancers, args.tss)
    
#     # Perform statistical comparison
#     comparison_results = perform_statistical_comparison(
#         enh_interactions, null_model, reference_condition, args.fdr_threshold
#     )
    
#     # Analyze distance dependence
#     distance_analysis = analyze_distance_dependence(enh_interactions)
    
#     # Create summary tables
#     summary_tables = create_summary_tables(enh_interactions, comparison_results, distance_analysis)
    
#     # Create plots
#     if not enh_interactions.empty:
#         create_summary_plots(enh_interactions, comparison_results, distance_analysis, args.output_prefix)
#     else:
#         print("No comparison data available for plotting")
    
#     # Save results
#     try:
#         if not distance_analysis.empty:
#             distance_analysis.to_csv(f"{args.output_prefix}_distance_analysis.tsv", sep='\t', index=False)
#             print(f"Saved distance analysis: {args.output_prefix}_distance_analysis.tsv")
        
#         if not enh_interactions.empty:
#             enh_interactions.to_csv(f"{args.output_prefix}_enhancer_interactions.tsv", sep='\t', index=False)
#             print(f"Saved enhancer interactions: {args.output_prefix}_enhancer_interactions.tsv")
        
#         if not comparison_results.empty:
#             comparison_results.to_csv(f"{args.output_prefix}_statistical_comparisons.tsv", sep='\t', index=False)
#             print(f"Saved statistical comparisons: {args.output_prefix}_statistical_comparisons.tsv")
        
#         # Save summary tables
#         if summary_tables:
#             with open(f"{args.output_prefix}_summary.txt", 'w') as f:
#                 f.write("Enhancer Class Analysis Summary\n")
#                 f.write("=" * 40 + "\n\n")
                
#                 for table_name, table_data in summary_tables.items():
#                     f.write(f"{table_name.upper()} SUMMARY:\n")
#                     for key, value in table_data.items():
#                         f.write(f"  {key}: {value}\n")
#                     f.write("\n")
#             print(f"Saved summary: {args.output_prefix}_summary.txt")
        
#     except Exception as e:
#         print(f"Error creating summary tables: {e}")
    
#     print(f"\nAnalysis complete! Results saved to {args.output_prefix}_*")
    
#     # Print final summary
#     print("\nFinal Summary:")
#     if not enh_interactions.empty:
#         print(f"  Total enhancer interactions found: {len(enh_interactions)}")
#         print(f"  E-E interactions: {sum(enh_interactions['interaction_type'] == 'E-E')}")
#         print(f"  E-TSS interactions: {sum(enh_interactions['interaction_type'] == 'E-TSS')}")
#     else:
#         print("  No enhancer interactions found")
    
#     if not comparison_results.empty:
#         print(f"  Statistical comparisons performed: {len(comparison_results)}")
#         print(f"  Significant comparisons: {sum(comparison_results['significant'])}")
#     else:
#         print("  No statistical comparisons performed")

# if __name__ == '__main__':
#     main()

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
    
    # Perform statistical comparison with proper null model handling
    comparison_results = perform_statistical_comparison(
        enh_interactions, null_model, reference_condition, args.fdr_threshold
    )
    
    # Analyze distance dependence
    distance_analysis = analyze_distance_dependence(enh_interactions)
    
    # Create summary tables
    summary_tables = create_summary_tables(enh_interactions, comparison_results, distance_analysis)
    
    # Create plots
    if not enh_interactions.empty:
        # Original summary plots
        create_summary_plots(enh_interactions, comparison_results, distance_analysis, args.output_prefix)
        
        # NEW: Null model comparison plots
        if not comparison_results.empty:
            create_null_comparison_plots(comparison_results, null_model, args.output_prefix)
            
            # Create summary with null model info
            null_summary = create_summary_with_null(comparison_results, null_model)
            
            # Save null model summary
            with open(f"{args.output_prefix}_null_model_summary.txt", 'w') as f:
                f.write("Null Model Comparison Summary\n")
                f.write("=" * 40 + "\n\n")
                for key, value in null_summary.items():
                    f.write(f"{key}: {value}\n")
            
            print(f"Saved null model summary: {args.output_prefix}_null_model_summary.txt")
        else:
            print("No comparison results for null model plotting")
    else:
        print("No enhancer interactions found for plotting")
    
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
            
            # NEW: Save significant results separately
            if 'any_significant' in comparison_results.columns:
                sig_results = comparison_results[comparison_results['any_significant']]
                if not sig_results.empty:
                    sig_results.to_csv(f"{args.output_prefix}_significant_null_comparisons.tsv", sep='\t', index=False)
                    print(f"Saved significant null comparisons: {args.output_prefix}_significant_null_comparisons.tsv")
                else:
                    print("No significant null comparisons found")
            
            # NEW: Save results by statistical method
            for method in ['z_test', 'permutation']:
                sig_col = f'p_value_{method}_significant'
                if sig_col in comparison_results.columns:
                    method_sig = comparison_results[comparison_results[sig_col]]
                    if not method_sig.empty:
                        method_sig.to_csv(f"{args.output_prefix}_significant_{method}_results.tsv", sep='\t', index=False)
                        print(f"Saved significant {method} results: {args.output_prefix}_significant_{method}_results.tsv")
        
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
        print(f"Error saving results: {e}")
        import traceback
        traceback.print_exc()
    
    print(f"\nAnalysis complete! Results saved to {args.output_prefix}_*")
    
    # Enhanced final summary with null model information
    print("\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    
    if not enh_interactions.empty:
        print(f"ENHANCER INTERACTIONS:")
        print(f"  Total enhancer interactions found: {len(enh_interactions)}")
        print(f"  E-E interactions: {sum(enh_interactions['interaction_type'] == 'E-E')}")
        print(f"  E-TSS interactions: {sum(enh_interactions['interaction_type'] == 'E-TSS')}")
        
        # Breakdown by enhancer class
        for contact_class in enh_interactions['contact_class'].unique():
            class_count = sum(enh_interactions['contact_class'] == contact_class)
            print(f"  {contact_class} interactions: {class_count}")
    else:
        print("ENHANCER INTERACTIONS: None found")
    
    if not comparison_results.empty:
        print(f"\nSTATISTICAL COMPARISONS:")
        print(f"  Total comparisons performed: {len(comparison_results)}")
        
        # Check for different significance columns
        if 'any_significant' in comparison_results.columns:
            any_sig = sum(comparison_results['any_significant'])
            print(f"  Significant by any method: {any_sig}")
            
        if 'p_value_z_test_significant' in comparison_results.columns:
            z_sig = sum(comparison_results['p_value_z_test_significant'])
            print(f"  Significant by Z-test: {z_sig}")
            
        if 'p_value_permutation_significant' in comparison_results.columns:
            perm_sig = sum(comparison_results['p_value_permutation_significant'])
            print(f"  Significant by permutation test: {perm_sig}")
            
        if 'biologically_significant' in comparison_results.columns:
            bio_sig = sum(comparison_results['biologically_significant'])
            print(f"  Biologically significant (>0.5 log2FC): {bio_sig}")
            
        # Effect size summary
        if 'effect_size' in comparison_results.columns:
            large_effects = sum(comparison_results['effect_size'] > 1)
            medium_effects = sum((comparison_results['effect_size'] > 0.5) & (comparison_results['effect_size'] <= 1))
            small_effects = sum((comparison_results['effect_size'] > 0.2) & (comparison_results['effect_size'] <= 0.5))
            print(f"  Large effects (>1 std): {large_effects}")
            print(f"  Medium effects (0.5-1 std): {medium_effects}")
            print(f"  Small effects (0.2-0.5 std): {small_effects}")
            
        # By comparison type
        if 'comparison' in comparison_results.columns:
            print(f"\nBY COMPARISON:")
            for comparison in comparison_results['comparison'].unique():
                comp_data = comparison_results[comparison_results['comparison'] == comparison]
                if 'any_significant' in comp_data.columns:
                    sig_count = sum(comp_data['any_significant'])
                    print(f"  {comparison}: {sig_count}/{len(comp_data)} significant")
                    
        # By interaction type
        if 'interaction_type' in comparison_results.columns:
            print(f"\nBY INTERACTION TYPE:")
            for int_type in comparison_results['interaction_type'].unique():
                type_data = comparison_results[comparison_results['interaction_type'] == int_type]
                if 'any_significant' in type_data.columns:
                    sig_count = sum(type_data['any_significant'])
                    print(f"  {int_type}: {sig_count}/{len(type_data)} significant")
                    
    else:
        print("\nSTATISTICAL COMPARISONS: None performed")
    
    # Null model summary
    if null_model is not None and not null_model.empty:
        print(f"\nNULL MODEL:")
        print(f"  Null model entries: {len(null_model)}")
        if 'logFC' in null_model.columns:
            null_logfc = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
            if len(null_logfc) > 0:
                print(f"  Null logFC mean: {null_logfc.mean():.3f}")
                print(f"  Null logFC std: {null_logfc.std():.3f}")
    else:
        print("\nNULL MODEL: Not provided")
    
    print("\n" + "="*60)
    print("FILES GENERATED:")
    expected_files = [
        f"{args.output_prefix}_enhancer_interactions.tsv",
        f"{args.output_prefix}_statistical_comparisons.tsv", 
        f"{args.output_prefix}_distance_analysis.tsv",
        f"{args.output_prefix}_summary_plots.pdf",
        f"{args.output_prefix}_null_model_comparison.pdf",
        f"{args.output_prefix}_null_model_summary.txt",
        f"{args.output_prefix}_summary.txt"
    ]
    
    for file in expected_files:
        print(f"  {file}")
    
    print("="*60)

if __name__ == '__main__':
    main()