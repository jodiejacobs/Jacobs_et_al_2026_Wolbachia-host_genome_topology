#!/usr/bin/env python3
"""
Analyze X-chromosome regulation changes in response to Wolbachia infection.
Focus on dosage compensation complex (DCC) high-affinity sites (HAS/CES).
Based on Schauer et al. 2017 methods.

Modified to include statistical comparisons against null model.
FIXED: Handles different column name formats between real and null data.
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

def detect_column_names(df):
    """
    Detect column names for chromosome coordinates from a dataframe.
    Returns dict with standardized column names.
    """
    cols = list(df.columns)
    
    # Detect chromosome columns
    if 'chr1' in cols and 'chr2' in cols:
        chr1_col, chr2_col = 'chr1', 'chr2'
        start1_col, end1_col = 'start1', 'end1'
        start2_col, end2_col = 'start2', 'end2'
    elif 'anchor1.chr' in cols and 'anchor2.chr' in cols:
        chr1_col, chr2_col = 'anchor1.chr', 'anchor2.chr'
        start1_col, end1_col = 'anchor1.start', 'anchor1.end'
        start2_col, end2_col = 'anchor2.start', 'anchor2.end'
    elif 'seqnames1' in cols and 'seqnames2' in cols:
        chr1_col, chr2_col = 'seqnames1', 'seqnames2'
        start1_col, end1_col = 'start1', 'end1'
        start2_col, end2_col = 'start2', 'end2'
    else:
        # Try to find by pattern
        chr_cols = [c for c in cols if 'chr' in c.lower() or 'seqname' in c.lower()]
        if len(chr_cols) >= 2:
            chr1_col, chr2_col = sorted(chr_cols)[:2]
            # Try to infer start/end columns
            start_cols = [c for c in cols if 'start' in c.lower()]
            end_cols = [c for c in cols if 'end' in c.lower()]
            if len(start_cols) >= 2 and len(end_cols) >= 2:
                start1_col, start2_col = sorted(start_cols)[:2]
                end1_col, end2_col = sorted(end_cols)[:2]
            else:
                raise ValueError(f"Could not identify start/end columns. Available: {cols}")
        else:
            raise ValueError(f"Could not identify chromosome columns. Available: {cols}")
    
    # Detect statistical columns
    fdr_col = None
    for col in ['FDR', 'padj', 'p.adjust', 'PValue', 'adj.P.Val']:
        if col in cols:
            fdr_col = col
            break
    
    logfc_col = None
    for col in ['logFC', 'log2FoldChange', 'logFC.up', 'log2FC']:
        if col in cols:
            logfc_col = col
            break
    
    return {
        'chr1_col': chr1_col,
        'chr2_col': chr2_col,
        'start1_col': start1_col,
        'end1_col': end1_col,
        'start2_col': start2_col,
        'end2_col': end2_col,
        'fdr_col': fdr_col,
        'logfc_col': logfc_col
    }

def analyze_x_interactions(interaction_file):
    """
    Analyze differential interactions specifically on the X chromosome.
    Handles different column name formats.
    """
    print(f"\nAnalyzing X interactions from {interaction_file}...")
    
    # Load interactions
    interactions = pd.read_csv(interaction_file)
    
    # Detect column names
    try:
        col_info = detect_column_names(interactions)
    except ValueError as e:
        print(f"Error: {e}")
        print(f"Available columns: {list(interactions.columns)}")
        raise
    
    print(f"Detected columns: chr={col_info['chr1_col']}, fdr={col_info['fdr_col']}, logfc={col_info['logfc_col']}")
    
    chr1_col = col_info['chr1_col']
    chr2_col = col_info['chr2_col']
    fdr_col = col_info['fdr_col']
    logfc_col = col_info['logfc_col']
    
    # Filter for X chromosome interactions
    if fdr_col and logfc_col:
        # With statistical filtering
        x_cis = interactions[
            (interactions[chr1_col] == 'X') & 
            (interactions[chr2_col] == 'X') &
            (interactions[fdr_col] < 0.01) &
            (abs(interactions[logfc_col]) > 1)
        ]
        
        x_trans = interactions[
            ((interactions[chr1_col] == 'X') | (interactions[chr2_col] == 'X')) &
            (interactions[chr1_col] != interactions[chr2_col]) &
            (interactions[fdr_col] < 0.01) &
            (abs(interactions[logfc_col]) > 1)
        ]
        
        total_sig = len(interactions[
            (interactions[fdr_col] < 0.01) &
            (abs(interactions[logfc_col]) > 1)
        ])
        
        mean_logfc_cis = x_cis[logfc_col].mean() if len(x_cis) > 0 else 0
        mean_logfc_trans = x_trans[logfc_col].mean() if len(x_trans) > 0 else 0
        upregulated_cis = sum(x_cis[logfc_col] > 0) if len(x_cis) > 0 else 0
        downregulated_cis = sum(x_cis[logfc_col] < 0) if len(x_cis) > 0 else 0
        upregulated_trans = sum(x_trans[logfc_col] > 0) if len(x_trans) > 0 else 0
        downregulated_trans = sum(x_trans[logfc_col] < 0) if len(x_trans) > 0 else 0
    else:
        # Without statistical filtering
        print(f"Warning: Missing statistical columns, counting all X interactions")
        x_cis = interactions[
            (interactions[chr1_col] == 'X') & 
            (interactions[chr2_col] == 'X')
        ]
        
        x_trans = interactions[
            ((interactions[chr1_col] == 'X') | (interactions[chr2_col] == 'X')) &
            (interactions[chr1_col] != interactions[chr2_col])
        ]
        
        total_sig = len(interactions)
        mean_logfc_cis = 0
        mean_logfc_trans = 0
        upregulated_cis = 0
        downregulated_cis = 0
        upregulated_trans = 0
        downregulated_trans = 0
    
    results = {
        'n_cis': len(x_cis),
        'n_trans': len(x_trans),
        'total_sig': total_sig,
        'mean_logfc_cis': mean_logfc_cis,
        'mean_logfc_trans': mean_logfc_trans,
        'upregulated_cis': upregulated_cis,
        'downregulated_cis': downregulated_cis,
        'upregulated_trans': upregulated_trans,
        'downregulated_trans': downregulated_trans,
        'col_info': col_info  # Store for use in other functions
    }
    
    print(f"Found {results['n_cis']} X-cis and {results['n_trans']} X-trans differential interactions")
    print(f"Total significant interactions: {results['total_sig']}")
    
    return results

def analyze_has_proximity(has_sites, interaction_file, window_size=50000, col_info=None):
    """
    Analyze differential interactions near HAS sites.
    """
    print(f"\nAnalyzing interactions near HAS sites (within {window_size}bp)...")
    
    # Load interactions
    interactions = pd.read_csv(interaction_file)
    
    # Detect or use provided column info
    if col_info is None:
        col_info = detect_column_names(interactions)
    
    chr1_col = col_info['chr1_col']
    chr2_col = col_info['chr2_col']
    start1_col = col_info['start1_col']
    end1_col = col_info['end1_col']
    start2_col = col_info['start2_col']
    end2_col = col_info['end2_col']
    fdr_col = col_info['fdr_col']
    logfc_col = col_info['logfc_col']
    
    # Filter for significant interactions if we have the columns
    if fdr_col and logfc_col:
        sig_interactions = interactions[
            (interactions[fdr_col] < 0.01) &
            (abs(interactions[logfc_col]) > 1)
        ].copy()
    else:
        sig_interactions = interactions.copy()
    
    # Convert to BED format for both anchors
    anchor1_bed = sig_interactions[[chr1_col, start1_col, end1_col]].copy()
    anchor1_bed.columns = ['chrom', 'start', 'end']
    anchor1_bed['interaction_idx'] = sig_interactions.index
    
    anchor2_bed = sig_interactions[[chr2_col, start2_col, end2_col]].copy()
    anchor2_bed.columns = ['chrom', 'start', 'end']
    anchor2_bed['interaction_idx'] = sig_interactions.index
    
    # Create HAS windows
    has_windows = has_sites.copy()
    has_windows['start'] = has_windows['start'] - window_size
    has_windows['end'] = has_windows['end'] + window_size
    has_windows = has_windows[has_windows['start'] >= 0]
    
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
            idx = int(overlap.fields[3])
            has_proximal_indices.add(idx)
        
        for overlap in anchor2_overlaps:
            idx = int(overlap.fields[3])
            has_proximal_indices.add(idx)
        
        # Get proximal interactions
        proximal_interactions = sig_interactions.loc[list(has_proximal_indices)]
        
        if logfc_col and len(proximal_interactions) > 0:
            mean_logfc = proximal_interactions[logfc_col].mean()
            upregulated = sum(proximal_interactions[logfc_col] > 0)
            downregulated = sum(proximal_interactions[logfc_col] < 0)
        else:
            mean_logfc = 0
            upregulated = 0
            downregulated = 0
        
        results = {
            'n_has_proximal': len(proximal_interactions),
            'has_proximal_upregulated': upregulated,
            'has_proximal_downregulated': downregulated,
            'mean_logfc_has_proximal': mean_logfc
        }
        
    except Exception as e:
        print(f"Warning: Could not analyze HAS proximity: {e}")
        results = {
            'n_has_proximal': 0,
            'has_proximal_upregulated': 0,
            'has_proximal_downregulated': 0,
            'mean_logfc_has_proximal': 0
        }
    
    print(f"Found {results['n_has_proximal']} interactions near HAS sites")
    
    return results

def test_has_enrichment_vs_null(has_proximity_real, has_proximity_null, 
                                 total_real, total_null):
    """
    Test if HAS-proximal interactions are enriched in real vs null data.
    Uses Fisher's exact test.
    """
    from scipy.stats import fisher_exact
    
    has_real = has_proximity_real['n_has_proximal']
    not_has_real = total_real - has_real
    
    has_null = has_proximity_null['n_has_proximal']
    not_has_null = total_null - has_null
    
    prop_real = has_real / total_real if total_real > 0 else 0
    prop_null = has_null / total_null if total_null > 0 else 0
    
    fold_enrichment = prop_real / prop_null if prop_null > 0 else np.inf
    
    contingency_table = [
        [has_real, not_has_real],
        [has_null, not_has_null]
    ]
    
    odds_ratio, pval = fisher_exact(contingency_table, alternative='greater')
    
    return {
        'real_has_proximal': has_real,
        'real_total': total_real,
        'real_proportion': prop_real,
        'null_has_proximal': has_null,
        'null_total': total_null,
        'null_proportion': prop_null,
        'fold_enrichment': fold_enrichment,
        'odds_ratio': odds_ratio,
        'pvalue': pval,
        'significant': pval < 0.05
    }

def test_x_enrichment_vs_null(x_interactions_real, x_interactions_null):
    """
    Test if X chromosome has more differential interactions than expected by null.
    """
    from scipy.stats import fisher_exact
    
    results = {}
    
    # Test X-cis
    x_cis_real = x_interactions_real['n_cis']
    total_real = x_interactions_real['total_sig']
    not_x_cis_real = total_real - x_cis_real
    
    x_cis_null = x_interactions_null['n_cis']
    total_null = x_interactions_null['total_sig']
    not_x_cis_null = total_null - x_cis_null
    
    contingency_cis = [
        [x_cis_real, not_x_cis_real],
        [x_cis_null, not_x_cis_null]
    ]
    
    odds_ratio_cis, pval_cis = fisher_exact(contingency_cis, alternative='greater')
    
    results['x_cis'] = {
        'real_count': x_cis_real,
        'real_total': total_real,
        'real_proportion': x_cis_real / total_real if total_real > 0 else 0,
        'null_count': x_cis_null,
        'null_total': total_null,
        'null_proportion': x_cis_null / total_null if total_null > 0 else 0,
        'odds_ratio': odds_ratio_cis,
        'pvalue': pval_cis,
        'significant': pval_cis < 0.05
    }
    
    # Test X-trans
    x_trans_real = x_interactions_real['n_trans']
    not_x_trans_real = total_real - x_trans_real
    
    x_trans_null = x_interactions_null['n_trans']
    not_x_trans_null = total_null - x_trans_null
    
    contingency_trans = [
        [x_trans_real, not_x_trans_real],
        [x_trans_null, not_x_trans_null]
    ]
    
    odds_ratio_trans, pval_trans = fisher_exact(contingency_trans, alternative='greater')
    
    results['x_trans'] = {
        'real_count': x_trans_real,
        'real_total': total_real,
        'real_proportion': x_trans_real / total_real if total_real > 0 else 0,
        'null_count': x_trans_null,
        'null_total': total_null,
        'null_proportion': x_trans_null / total_null if total_null > 0 else 0,
        'odds_ratio': odds_ratio_trans,
        'pvalue': pval_trans,
        'significant': pval_trans < 0.05
    }
    
    return results

def create_x_regulation_plots(has_results, compartment_results, interaction_results, 
                               has_proximity_results, output_prefix, 
                               interaction_results_null=None, has_proximity_results_null=None,
                               has_enrichment_stats=None, x_enrichment_stats=None):
    """
    Create visualizations for X chromosome regulation analysis.
    """
    # Determine layout
    if interaction_results_null is not None:
        fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    else:
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
    
    # Additional plots if null data present
    if interaction_results_null is not None and axes.shape[0] > 2:
        # Plot 7: Real vs Null - HAS proximity
        ax = axes[2, 0]
        if has_enrichment_stats:
            labels = ['Real', 'Null']
            proportions = [
                has_enrichment_stats['real_proportion'] * 100,
                has_enrichment_stats['null_proportion'] * 100
            ]
            colors = ['steelblue', 'lightgray']
            ax.bar(labels, proportions, color=colors)
            
            if has_enrichment_stats['significant']:
                max_y = max(proportions) * 1.1
                pval = has_enrichment_stats['pvalue']
                stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*'
                ax.text(0.5, max_y, stars, ha='center', fontsize=16)
            
            ax.set_ylabel('% Interactions near HAS')
            ax.set_title(f'HAS Proximity Enrichment\n(Fold: {has_enrichment_stats["fold_enrichment"]:.2f}x)')
            ax.set_ylim(0, max(proportions) * 1.2 if proportions else 1)
        
        # Plot 8: Real vs Null - X-cis
        ax = axes[2, 1]
        if x_enrichment_stats and 'x_cis' in x_enrichment_stats:
            labels = ['Real', 'Null']
            proportions = [
                x_enrichment_stats['x_cis']['real_proportion'] * 100,
                x_enrichment_stats['x_cis']['null_proportion'] * 100
            ]
            colors = ['steelblue', 'lightgray']
            ax.bar(labels, proportions, color=colors)
            
            if x_enrichment_stats['x_cis']['significant']:
                max_y = max(proportions) * 1.1
                pval = x_enrichment_stats['x_cis']['pvalue']
                stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*'
                ax.text(0.5, max_y, stars, ha='center', fontsize=16)
            
            ax.set_ylabel('% X-cis Interactions')
            ax.set_title(f'X-Cis Enrichment\n(OR: {x_enrichment_stats["x_cis"]["odds_ratio"]:.2f})')
            ax.set_ylim(0, max(proportions) * 1.2 if proportions else 1)
        
        # Plot 9: Real vs Null - X-trans
        ax = axes[2, 2]
        if x_enrichment_stats and 'x_trans' in x_enrichment_stats:
            labels = ['Real', 'Null']
            proportions = [
                x_enrichment_stats['x_trans']['real_proportion'] * 100,
                x_enrichment_stats['x_trans']['null_proportion'] * 100
            ]
            colors = ['steelblue', 'lightgray']
            ax.bar(labels, proportions, color=colors)
            
            if x_enrichment_stats['x_trans']['significant']:
                max_y = max(proportions) * 1.1
                pval = x_enrichment_stats['x_trans']['pvalue']
                stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*'
                ax.text(0.5, max_y, stars, ha='center', fontsize=16)
            
            ax.set_ylabel('% X-trans Interactions')
            ax.set_title(f'X-Trans Enrichment\n(OR: {x_enrichment_stats["x_trans"]["odds_ratio"]:.2f})')
            ax.set_ylim(0, max(proportions) * 1.2 if proportions else 1)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_x_regulation_analysis.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_prefix}_x_regulation_analysis.pdf")

def main():
    parser = argparse.ArgumentParser(description='Analyze X chromosome regulation')
    parser.add_argument('--has_sites', required=True, help='HAS/CES sites BED file')
    parser.add_argument('--mcool_files', nargs='+', required=True, help='Micro-C mcool files')
    parser.add_argument('--conditions', nargs='+', required=True, help='Condition names')
    parser.add_argument('--interaction_file', help='Differential interaction file')
    parser.add_argument('--null_interaction_file', help='Null model interaction file')
    parser.add_argument('--resolution', type=int, default=5000, help='Resolution for analysis')
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
    
    # Initialize results
    interaction_results = {}
    has_proximity_results = {}
    interaction_results_null = None
    has_proximity_results_null = None
    has_enrichment_stats = None
    x_enrichment_stats = None
    
    # Analyze real data
    if args.interaction_file:
        print("\n" + "="*60)
        print("ANALYZING REAL DATA")
        print("="*60)
        interaction_results = analyze_x_interactions(args.interaction_file)
        has_proximity_results = analyze_has_proximity(
            has_sites, args.interaction_file, 
            col_info=interaction_results.get('col_info')
        )
    
    # Analyze null model
    if args.null_interaction_file and args.interaction_file:
        print("\n" + "="*60)
        print("ANALYZING NULL MODEL")
        print("="*60)
        interaction_results_null = analyze_x_interactions(args.null_interaction_file)
        has_proximity_results_null = analyze_has_proximity(
            has_sites, args.null_interaction_file,
            col_info=interaction_results_null.get('col_info')
        )
        
        # Statistical comparisons
        print("\n" + "="*60)
        print("STATISTICAL COMPARISONS: REAL vs NULL")
        print("="*60)
        
        # Test HAS enrichment
        has_enrichment_stats = test_has_enrichment_vs_null(
            has_proximity_results, has_proximity_results_null,
            interaction_results['total_sig'], interaction_results_null['total_sig']
        )
        
        print("\n--- HAS Proximity Enrichment ---")
        print(f"Real data: {has_enrichment_stats['real_has_proximal']}/{has_enrichment_stats['real_total']}")
        print(f"  Proportion: {has_enrichment_stats['real_proportion']:.4f} ({has_enrichment_stats['real_proportion']*100:.2f}%)")
        print(f"Null model: {has_enrichment_stats['null_has_proximal']}/{has_enrichment_stats['null_total']}")
        print(f"  Proportion: {has_enrichment_stats['null_proportion']:.4f} ({has_enrichment_stats['null_proportion']*100:.2f}%)")
        print(f"Enrichment:")
        print(f"  Fold: {has_enrichment_stats['fold_enrichment']:.2f}x")
        print(f"  Odds ratio: {has_enrichment_stats['odds_ratio']:.2f}")
        print(f"  P-value: {has_enrichment_stats['pvalue']:.4e}")
        print(f"  Significant: {'YES ***' if has_enrichment_stats['significant'] else 'NO'}")
        
        # Test X enrichment
        x_enrichment_stats = test_x_enrichment_vs_null(
            interaction_results, interaction_results_null
        )
        
        print("\n--- X-Cis Enrichment ---")
        print(f"Real data: {x_enrichment_stats['x_cis']['real_count']}/{x_enrichment_stats['x_cis']['real_total']}")
        print(f"  Proportion: {x_enrichment_stats['x_cis']['real_proportion']:.4f}")
        print(f"Null model: {x_enrichment_stats['x_cis']['null_count']}/{x_enrichment_stats['x_cis']['null_total']}")
        print(f"  Proportion: {x_enrichment_stats['x_cis']['null_proportion']:.4f}")
        print(f"Enrichment:")
        print(f"  Odds ratio: {x_enrichment_stats['x_cis']['odds_ratio']:.2f}")
        print(f"  P-value: {x_enrichment_stats['x_cis']['pvalue']:.4e}")
        print(f"  Significant: {'YES ***' if x_enrichment_stats['x_cis']['significant'] else 'NO'}")
        
        print("\n--- X-Trans Enrichment ---")
        print(f"Real data: {x_enrichment_stats['x_trans']['real_count']}/{x_enrichment_stats['x_trans']['real_total']}")
        print(f"  Proportion: {x_enrichment_stats['x_trans']['real_proportion']:.4f}")
        print(f"Null model: {x_enrichment_stats['x_trans']['null_count']}/{x_enrichment_stats['x_trans']['null_total']}")
        print(f"  Proportion: {x_enrichment_stats['x_trans']['null_proportion']:.4f}")
        print(f"Enrichment:")
        print(f"  Odds ratio: {x_enrichment_stats['x_trans']['odds_ratio']:.2f}")
        print(f"  P-value: {x_enrichment_stats['x_trans']['pvalue']:.4e}")
        print(f"  Significant: {'YES ***' if x_enrichment_stats['x_trans']['significant'] else 'NO'}")
        
        # Save statistics
        comparison_summary = {
            'has_enrichment_fold': has_enrichment_stats['fold_enrichment'],
            'has_enrichment_odds_ratio': has_enrichment_stats['odds_ratio'],
            'has_enrichment_pvalue': has_enrichment_stats['pvalue'],
            'has_enrichment_significant': has_enrichment_stats['significant'],
            'x_cis_odds_ratio': x_enrichment_stats['x_cis']['odds_ratio'],
            'x_cis_pvalue': x_enrichment_stats['x_cis']['pvalue'],
            'x_cis_significant': x_enrichment_stats['x_cis']['significant'],
            'x_trans_odds_ratio': x_enrichment_stats['x_trans']['odds_ratio'],
            'x_trans_pvalue': x_enrichment_stats['x_trans']['pvalue'],
            'x_trans_significant': x_enrichment_stats['x_trans']['significant']
        }
        
        pd.DataFrame([comparison_summary]).to_csv(
            f"{args.output_prefix}_null_comparison_stats.tsv", 
            sep='\t', index=False
        )
        
        # Detailed comparison
        detailed_comparison = {
            'real_has_proximal': has_enrichment_stats['real_has_proximal'],
            'real_total_interactions': has_enrichment_stats['real_total'],
            'real_has_proportion': has_enrichment_stats['real_proportion'],
            'null_has_proximal': has_enrichment_stats['null_has_proximal'],
            'null_total_interactions': has_enrichment_stats['null_total'],
            'null_has_proportion': has_enrichment_stats['null_proportion'],
            'has_fold_enrichment': has_enrichment_stats['fold_enrichment'],
            'has_odds_ratio': has_enrichment_stats['odds_ratio'],
            'has_pvalue': has_enrichment_stats['pvalue'],
            'real_x_cis': x_enrichment_stats['x_cis']['real_count'],
            'real_x_cis_proportion': x_enrichment_stats['x_cis']['real_proportion'],
            'null_x_cis': x_enrichment_stats['x_cis']['null_count'],
            'null_x_cis_proportion': x_enrichment_stats['x_cis']['null_proportion'],
            'x_cis_odds_ratio': x_enrichment_stats['x_cis']['odds_ratio'],
            'x_cis_pvalue': x_enrichment_stats['x_cis']['pvalue'],
            'real_x_trans': x_enrichment_stats['x_trans']['real_count'],
            'real_x_trans_proportion': x_enrichment_stats['x_trans']['real_proportion'],
            'null_x_trans': x_enrichment_stats['x_trans']['null_count'],
            'null_x_trans_proportion': x_enrichment_stats['x_trans']['null_proportion'],
            'x_trans_odds_ratio': x_enrichment_stats['x_trans']['odds_ratio'],
            'x_trans_pvalue': x_enrichment_stats['x_trans']['pvalue']
        }
        
        pd.DataFrame([detailed_comparison]).to_csv(
            f"{args.output_prefix}_null_comparison_detailed.tsv",
            sep='\t', index=False
        )
    
    # Create visualizations
    create_x_regulation_plots(
        has_results, compartment_results, interaction_results, 
        has_proximity_results, args.output_prefix,
        interaction_results_null, has_proximity_results_null,
        has_enrichment_stats, x_enrichment_stats
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
    
    # Add interaction results
    if interaction_results:
        interaction_summary = {
            'x_cis_interactions': interaction_results['n_cis'],
            'x_trans_interactions': interaction_results['n_trans'],
            'total_sig_interactions': interaction_results['total_sig'],
            'x_cis_upregulated': interaction_results['upregulated_cis'],
            'x_cis_downregulated': interaction_results['downregulated_cis'],
            'x_trans_upregulated': interaction_results['upregulated_trans'],
            'x_trans_downregulated': interaction_results['downregulated_trans'],
            'mean_logfc_x_cis': interaction_results['mean_logfc_cis'],
            'mean_logfc_x_trans': interaction_results['mean_logfc_trans']
        }
        if summary:
            summary[0].update(interaction_summary)
    
    if has_proximity_results:
        proximity_summary = {
            'has_proximal_interactions': has_proximity_results['n_has_proximal'],
            'has_proximal_upregulated': has_proximity_results['has_proximal_upregulated'],
            'has_proximal_downregulated': has_proximity_results['has_proximal_downregulated'],
            'mean_logfc_has_proximal': has_proximity_results['mean_logfc_has_proximal']
        }
        if summary:
            summary[0].update(proximity_summary)
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(f"{args.output_prefix}_x_regulation_summary.tsv", sep='\t', index=False)
    
    # Save detailed results
    if has_results:
        all_has_results = []
        for condition, df in has_results.items():
            df_copy = df.copy()
            df_copy['condition'] = condition
            all_has_results.append(df_copy)
        
        if all_has_results:
            combined_has = pd.concat(all_has_results, ignore_index=True)
            combined_has.to_csv(f"{args.output_prefix}_has_contact_comparison.tsv", sep='\t', index=False)
    
    if compartment_results:
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
        comp_df.to_csv(f"{args.output_prefix}_x_compartment_changes.tsv", sep='\t', index=False)
    
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"\nResults saved to {args.output_prefix}_*")
    
    # Print summary
    print("\nSummary:")
    print(f"  HAS sites analyzed: {len(has_sites)}")
    print(f"  Conditions with HAS data: {len(has_results)}")
    print(f"  Conditions with compartment data: {len(compartment_results)}")
    
    if interaction_results:
        print(f"\nReal data interactions:")
        print(f"  Total significant: {interaction_results['total_sig']}")
        print(f"  X-chromosome cis: {interaction_results['n_cis']}")
        print(f"  X-chromosome trans: {interaction_results['n_trans']}")
    
    if has_proximity_results:
        print(f"  HAS-proximal: {has_proximity_results['n_has_proximal']}")
    
    if interaction_results_null:
        print(f"\nNull model interactions:")
        print(f"  Total significant: {interaction_results_null['total_sig']}")
        print(f"  X-chromosome cis: {interaction_results_null['n_cis']}")
        print(f"  X-chromosome trans: {interaction_results_null['n_trans']}")
        
    if has_proximity_results_null:
        print(f"  HAS-proximal: {has_proximity_results_null['n_has_proximal']}")

if __name__ == '__main__':
    main()
