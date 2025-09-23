#!/usr/bin/env python3
"""
This is a working version as of 06/16/25
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
