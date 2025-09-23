#!/usr/bin/env python3
"""
Analyze enrichment of architectural proteins at chromatin features across conditions.
Compare DOX (uninfected) vs wMel, wRi, wWil with null model testing and FDR correction.
Uses existing null model and differential interactions from diffHic analysis.
"""
#Test with:
# python scripts/architectural_enrichment.py     --interactions /private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/structural_analysis_ms/results/diffhic_results/summary/all_results_combined.csv     --null_model /private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/structural_analysis_ms/results/diffhic_results/res_1000/null_model_results.csv     --chip_dir /private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/structural_analysis_ms/chip_peaks     --genome /private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/structural_analysis_ms/reference_files/dm6.genome     --output_prefix /private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/structural_analysis_ms/results/enrichment/architectural

import pandas as pd
import pybedtools
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict

# Define architectural proteins categories
ARCHITECTURAL_PROTEINS = {
    'DNA_binding': ['CTCF', 'Su(Hw)', 'BEAF-32', 'DREF', 'TFIIIC', 'Z4', 'Elba', 'ZIPIC', 'Ibf1', 'Ibf2'],
    'accessory': ['CP190', 'Mod(mdg4)', 'Rad21', 'Cap-H2', 'Fs(1)h-L', 'L3mbt', 'Chromator']
}

def scan_chip_directory(chip_dir):
    """Scan ChIP directory and suggest file mappings"""
    chip_dir_path = Path(chip_dir)
    # Look for all bed files including compressed ones
    bed_files = []
    bed_files.extend(chip_dir_path.glob("*.bed"))
    bed_files.extend(chip_dir_path.glob("*.bed.gz"))
    bed_files.extend(chip_dir_path.glob("*.narrowPeak"))
    bed_files.extend(chip_dir_path.glob("*.narrowPeak.gz"))
    
    print(f"\n=== CHIP DIRECTORY SCAN ===")

def load_chip_peaks(chip_dir):
    """Load ChIP-seq peaks for architectural proteins"""
    chip_data = {}
    
    # Create mapping from protein names to file patterns
    protein_patterns = {
        'CTCF': ['dCTCF', 'CTCF'],
        'Su(Hw)': ['SuHw', 'Su_Hw'],
        'BEAF-32': ['BEAF'],
        'DREF': ['DREF'],
        'TFIIIC': ['dTFIIIC', 'TFIIIC'],
        'Z4': ['Z4'],
        'Elba': ['Elba'],
        'ZIPIC': ['ZIPIC'],
        'Ibf1': ['Ibf1'],
        'Ibf2': ['Ibf2'],
        'CP190': ['CP190'],
        'Mod(mdg4)': ['Mod', 'mdg4'],
        'Rad21': ['Rad21'],
        'Cap-H2': ['CAPH2', 'Cap-H2'],
        'Fs(1)h-L': ['Fs1h'],
        'L3mbt': ['L3mbt'],
        'Chromator': ['Chromator']
    }
    
    # Get all bed files in the directory
    chip_dir_path = Path(chip_dir)
    bed_files = []
    bed_files.extend(chip_dir_path.glob("*.bed"))
    bed_files.extend(chip_dir_path.glob("*.bed.gz"))
    bed_files.extend(chip_dir_path.glob("*.narrowPeak"))
    bed_files.extend(chip_dir_path.glob("*.narrowPeak.gz"))
    
    print(f"Found {len(bed_files)} BED files in {chip_dir}")
    
    for protein, patterns in protein_patterns.items():
        found = False
        for pattern in patterns:
            # Look for files containing the pattern (case insensitive)
            matching_files = [f for f in bed_files if pattern.lower() in f.name.lower()]
            
            if matching_files:
                # Take the first matching file, prefer _peaks.bed files
                peaks_files = [f for f in matching_files if '_peaks.bed' in f.name.lower()]
                if peaks_files:
                    chip_file = peaks_files[0]
                else:
                    chip_file = matching_files[0]
                
                print(f"Loading {protein} peaks from {chip_file.name}")
                
                # Handle compressed files
                if chip_file.suffix == '.gz':
                    # For compressed files, we need to decompress first or use zcat
                    import gzip
                    try:
                        with gzip.open(chip_file, 'rt') as f:
                            # Read the file content and create a temporary BedTool
                            content = f.read()
                            chip_data[protein] = pybedtools.BedTool(content, from_string=True)
                    except Exception as e:
                        print(f"Warning: Could not read compressed file {chip_file}: {e}")
                        continue
                else:
                    chip_data[protein] = pybedtools.BedTool(str(chip_file))
                
                found = True
                break
        
        if not found:
            print(f"Warning: No ChIP data found for {protein}")
            # List some files for debugging
            if len(bed_files) > 0:
                print(f"  Available files (first 5): {[f.name for f in bed_files[:5]]}")
    
    return chip_data
    

def load_null_model(null_model_file):
    """Load null model results from diffHic analysis"""
    print(f"Loading null model from {null_model_file}")
    null_df = pd.read_csv(null_model_file)
    
    # Extract relevant statistics from null model
    null_stats = {
        'mean_logCPM': null_df['logCPM'].mean(),
        'std_logCPM': null_df['logCPM'].std(),
        'mean_pvalue': null_df['PValue'].mean(),
        'background_interactions': len(null_df)
    }
    
    print(f"Null model stats: {null_stats}")
    return null_stats

def load_differential_interactions(interactions_file, fdr_threshold=0.05):
    """Load differential interactions from combined diffHic results"""
    print(f"Loading differential interactions from {interactions_file}")
    
    # Read the combined results
    interactions_df = pd.read_csv(interactions_file)
    
    # Filter for significant interactions
    significant_interactions = interactions_df[
        (interactions_df['FDR'] < fdr_threshold) & 
        (abs(interactions_df['logFC']) > 1)
    ].copy()
    
    print(f"Found {len(significant_interactions)} significant differential interactions (FDR < {fdr_threshold})")
    
    # Convert to BED format for each condition
    features = {
        'interactions': {},
        'all_interactions': {}  # Include all interactions for background
    }
    
    # Group by infection condition
    for infection in significant_interactions['infection'].unique():
        infection_data = significant_interactions[significant_interactions['infection'] == infection]
        
        # Create BED entries for interaction anchors
        interaction_bed = []
        for _, row in infection_data.iterrows():
            # Add both anchors of each interaction
            interaction_bed.append([row['chr1'], row['start1'], row['end1']])
            interaction_bed.append([row['chr2'], row['start2'], row['end2']])
        
        if interaction_bed:
            features['interactions'][infection] = pybedtools.BedTool(interaction_bed)
            print(f"  {infection}: {len(interaction_bed)} interaction anchors")
    
    # Also create background set from all interactions (not just significant)
    all_interaction_bed = []
    for _, row in interactions_df.iterrows():
        all_interaction_bed.append([row['chr1'], row['start1'], row['end1']])
        all_interaction_bed.append([row['chr2'], row['start2'], row['end2']])
    
    features['all_interactions']['background'] = pybedtools.BedTool(all_interaction_bed)
    print(f"Background set: {len(all_interaction_bed)} total interaction anchors")
    
    return features, significant_interactions

def calculate_enrichment_with_existing_null(query_bed, reference_bed, genome_file, 
                                          null_stats, n_permutations=1000):
    """
    Calculate enrichment using existing null model statistics and permutation test.
    """
    # Calculate observed overlap
    observed = query_bed.intersect(reference_bed, u=True).count()
    
    # Use null model statistics to estimate expected background
    # This is a simplified approach - you might want to adjust based on your null model
    query_size = query_bed.count()
    reference_size = reference_bed.count()
    
    # Estimate expected based on null model background rate
    # Adjust this calculation based on your specific null model design
    background_rate = null_stats.get('background_interactions', 1000)
    expected_rate = min(0.1, reference_size / background_rate)  # Cap at 10% overlap
    expected = query_size * expected_rate
    
    # Permutation test for empirical p-value calculation
    null_overlaps = []
    for i in range(n_permutations):
        # Shuffle query regions maintaining chromosome
        try:
            shuffled = query_bed.shuffle(g=genome_file, chrom=True, seed=i)
            null_overlap = shuffled.intersect(reference_bed, u=True).count()
            null_overlaps.append(null_overlap)
        except:
            # If shuffling fails, use a random sample
            null_overlaps.append(np.random.poisson(expected))
    
    # Calculate statistics
    null_overlaps = np.array(null_overlaps)
    empirical_expected = np.median(null_overlaps)
    
    # Use empirical expected if available, otherwise use model-based
    final_expected = empirical_expected if empirical_expected > 0 else expected
    
    if final_expected > 0:
        enrichment = observed / final_expected
        log2_enrichment = np.log2(enrichment)
    else:
        enrichment = np.inf if observed > 0 else 1
        log2_enrichment = np.inf if observed > 0 else 0
    
    # Calculate p-value (two-tailed)
    if len(null_overlaps) > 0:
        if observed > final_expected:
            p_value = (np.sum(null_overlaps >= observed) + 1) / (len(null_overlaps) + 1)
        else:
            p_value = (np.sum(null_overlaps <= observed) + 1) / (len(null_overlaps) + 1)
        p_value = min(p_value * 2, 1)  # Two-tailed
    else:
        # Fallback to binomial test if permutation fails
        p_value = stats.binom_test(observed, query_size, expected_rate)
    
    return {
        'observed': observed,
        'expected': final_expected,
        'enrichment': enrichment,
        'log2_enrichment': log2_enrichment,
        'p_value': p_value,
        'null_mean': np.mean(null_overlaps) if len(null_overlaps) > 0 else final_expected,
        'null_std': np.std(null_overlaps) if len(null_overlaps) > 0 else 1
    }

def compare_enrichment_across_conditions(features, chip_data, genome_file, null_stats, window_size=5000):
    """
    Compare architectural protein enrichment between conditions using existing data.
    """
    results = defaultdict(list)
    
    print("\nAnalyzing differential interaction enrichment...")
    
    # Analyze enrichment at differential interaction sites
    for infected_condition in features['interactions'].keys():
        print(f"  Analyzing {infected_condition} interactions...")
        
        int_features = features['interactions'][infected_condition]
        int_windows = int_features.slop(b=window_size, g=genome_file)
        
        for protein, peaks in chip_data.items():
            print(f"    Analyzing {protein}...")
            
            enrichment = calculate_enrichment_with_existing_null(
                int_windows, peaks, genome_file, null_stats
            )
            
            results['differential_interactions'].append({
                'feature_type': 'differential_interactions',
                'protein': protein,
                'condition': infected_condition,
                'comparison': f'DOX_vs_{infected_condition}',
                'observed': enrichment['observed'],
                'expected': enrichment['expected'],
                'enrichment': enrichment['enrichment'],
                'log2_enrichment': enrichment['log2_enrichment'],
                'p_value': enrichment['p_value'],
                'null_mean': enrichment['null_mean'],
                'null_std': enrichment['null_std']
            })
    
    # Also analyze enrichment at background interaction sites for comparison
    print("\nAnalyzing background interaction enrichment...")
    if 'background' in features['all_interactions']:
        bg_features = features['all_interactions']['background']
        bg_windows = bg_features.slop(b=window_size, g=genome_file)
        
        for protein, peaks in chip_data.items():
            enrichment = calculate_enrichment_with_existing_null(
                bg_windows, peaks, genome_file, null_stats
            )
            
            results['background_interactions'].append({
                'feature_type': 'background_interactions',
                'protein': protein,
                'condition': 'background',
                'comparison': 'background',
                'observed': enrichment['observed'],
                'expected': enrichment['expected'],
                'enrichment': enrichment['enrichment'],
                'log2_enrichment': enrichment['log2_enrichment'],
                'p_value': enrichment['p_value'],
                'null_mean': enrichment['null_mean'],
                'null_std': enrichment['null_std']
            })
    
    return results

def apply_fdr_correction(results):
    """Apply FDR correction to all p-values"""
    corrected_results = {}
    
    for feature_type, data in results.items():
        df = pd.DataFrame(data)
        
        if 'p_value' in df.columns and len(df) > 0:
            # Apply FDR correction within each feature type
            _, fdr_values, _, _ = multipletests(df['p_value'], method='fdr_bh')
            df['fdr'] = fdr_values
            df['significant'] = df['fdr'] < 0.05
            
            # Also add a less stringent threshold
            df['suggestive'] = df['fdr'] < 0.1
        
        corrected_results[feature_type] = df
    
    return corrected_results

def analyze_protein_categories(corrected_results):
    """Analyze enrichment patterns by protein category"""
    category_analysis = []
    
    for feature_type, df in corrected_results.items():
        if len(df) == 0:
            continue
            
        # Add protein category
        df['category'] = df['protein'].apply(
            lambda x: 'DNA_binding' if x in ARCHITECTURAL_PROTEINS['DNA_binding'] 
            else 'accessory' if x in ARCHITECTURAL_PROTEINS['accessory'] 
            else 'other'
        )
        
        # Summarize by category and condition
        for condition in df['condition'].unique():
            cond_data = df[df['condition'] == condition]
            
            for category in ['DNA_binding', 'accessory']:
                cat_data = cond_data[cond_data['category'] == category]
                
                if len(cat_data) > 0:
                    category_analysis.append({
                        'feature_type': feature_type,
                        'condition': condition,
                        'category': category,
                        'n_proteins': len(cat_data),
                        'n_significant': cat_data['significant'].sum(),
                        'n_suggestive': cat_data['suggestive'].sum(),
                        'mean_log2_enrichment': cat_data['log2_enrichment'].mean(),
                        'median_log2_enrichment': cat_data['log2_enrichment'].median(),
                        'proteins_enriched': cat_data[
                            (cat_data['significant']) & 
                            (cat_data['log2_enrichment'] > 0)
                        ]['protein'].tolist(),
                        'proteins_depleted': cat_data[
                            (cat_data['significant']) & 
                            (cat_data['log2_enrichment'] < 0)
                        ]['protein'].tolist()
                    })
    
    return pd.DataFrame(category_analysis)

def create_enrichment_heatmaps(corrected_results, output_prefix):
    """Create comprehensive heatmaps for enrichment results"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    plot_idx = 0
    
    # Plot 1: Differential interactions enrichment heatmap
    if 'differential_interactions' in corrected_results:
        df = corrected_results['differential_interactions']
        
        if len(df) > 0:
            # Create pivot table for heatmap
            pivot_df = df.pivot_table(
                values='log2_enrichment',
                index='protein',
                columns='condition',
                fill_value=0
            )
            
            # Create significance annotation
            sig_pivot = df.pivot_table(
                values='significant',
                index='protein',
                columns='condition',
                fill_value=False
            )
            
            # Create annotation matrix
            annot = pivot_df.copy().astype(str)
            for i in range(len(annot)):
                for j in range(len(annot.columns)):
                    value = pivot_df.iloc[i, j]
                    if sig_pivot.iloc[i, j]:
                        annot.iloc[i, j] = f"{value:.2f}*"
                    else:
                        annot.iloc[i, j] = f"{value:.2f}"
            
            sns.heatmap(pivot_df, 
                       cmap='RdBu_r', center=0,
                       annot=annot, fmt='',
                       cbar_kws={'label': 'Log2 Enrichment'},
                       ax=axes[plot_idx])
            axes[plot_idx].set_title('Differential Interactions Enrichment')
            axes[plot_idx].set_ylabel('Architectural Protein')
            axes[plot_idx].set_xlabel('Wolbachia Strain')
    
    plot_idx += 1
    
    # Plot 2: Comparison with background
    if 'background_interactions' in corrected_results and 'differential_interactions' in corrected_results:
        bg_df = corrected_results['background_interactions']
        diff_df = corrected_results['differential_interactions']
        
        # Compare mean enrichments
        comparison_data = []
        for protein in set(bg_df['protein'].tolist() + diff_df['protein'].tolist()):
            bg_enrich = bg_df[bg_df['protein'] == protein]['log2_enrichment'].mean()
            
            # Average across strains for differential
            diff_enrich = diff_df[diff_df['protein'] == protein]['log2_enrichment'].mean()
            
            comparison_data.append({
                'protein': protein,
                'background': bg_enrich,
                'differential': diff_enrich,
                'difference': diff_enrich - bg_enrich
            })
        
        comp_df = pd.DataFrame(comparison_data)
        
        if len(comp_df) > 0:
            # Sort by difference
            comp_df = comp_df.sort_values('difference', ascending=True)
            
            # Create bar plot
            x = np.arange(len(comp_df))
            width = 0.35
            
            axes[plot_idx].bar(x - width/2, comp_df['background'], width, 
                             label='Background', alpha=0.7)
            axes[plot_idx].bar(x + width/2, comp_df['differential'], width, 
                             label='Differential', alpha=0.7)
            
            axes[plot_idx].set_xlabel('Architectural Protein')
            axes[plot_idx].set_ylabel('Mean Log2 Enrichment')
            axes[plot_idx].set_title('Background vs Differential Enrichment')
            axes[plot_idx].set_xticks(x)
            axes[plot_idx].set_xticklabels(comp_df['protein'], rotation=45, ha='right')
            axes[plot_idx].legend()
            axes[plot_idx].axhline(y=0, color='black', linestyle='--', alpha=0.3)
    
    plot_idx += 1
    
    # Plot 3: Protein category summary
    if 'differential_interactions' in corrected_results:
        df = corrected_results['differential_interactions']
        df['category'] = df['protein'].apply(
            lambda x: 'DNA_binding' if x in ARCHITECTURAL_PROTEINS['DNA_binding'] 
            else 'accessory'
        )
        
        category_stats = []
        for condition in df['condition'].unique():
            cond_data = df[df['condition'] == condition]
            
            for category in ['DNA_binding', 'accessory']:
                cat_data = cond_data[cond_data['category'] == category]
                if len(cat_data) > 0:
                    category_stats.append({
                        'condition': condition,
                        'category': category,
                        'mean_enrichment': cat_data['log2_enrichment'].mean(),
                        'n_significant': cat_data['significant'].sum(),
                        'n_total': len(cat_data)
                    })
        
        if category_stats:
            cat_df = pd.DataFrame(category_stats)
            
            # Create grouped bar plot
            pivot_enrich = cat_df.pivot(index='condition', columns='category', values='mean_enrichment')
            pivot_enrich.plot(kind='bar', ax=axes[plot_idx])
            axes[plot_idx].set_title('Mean Enrichment by Protein Category')
            axes[plot_idx].set_ylabel('Mean Log2 Enrichment')
            axes[plot_idx].set_xlabel('Wolbachia Strain')
            axes[plot_idx].legend(title='Protein Category')
            axes[plot_idx].axhline(y=0, color='black', linestyle='--', alpha=0.3)
    
    plot_idx += 1
    
    # Plot 4: Significance summary
    if 'differential_interactions' in corrected_results:
        df = corrected_results['differential_interactions']
        
        # Count significant proteins per condition
        sig_counts = df[df['significant']].groupby('condition').size()
        total_counts = df.groupby('condition').size()
        
        conditions = list(set(sig_counts.index.tolist() + total_counts.index.tolist()))
        sig_values = [sig_counts.get(c, 0) for c in conditions]
        total_values = [total_counts.get(c, 0) for c in conditions]
        
        x = np.arange(len(conditions))
        axes[plot_idx].bar(x, total_values, alpha=0.5, label='Total')
        axes[plot_idx].bar(x, sig_values, alpha=0.8, label='Significant')
        
        axes[plot_idx].set_xlabel('Wolbachia Strain')
        axes[plot_idx].set_ylabel('Number of Proteins')
        axes[plot_idx].set_title('Significant Enrichments per Condition')
        axes[plot_idx].set_xticks(x)
        axes[plot_idx].set_xticklabels(conditions)
        axes[plot_idx].legend()
        
        # Add percentages as text
        for i, (sig, total) in enumerate(zip(sig_values, total_values)):
            if total > 0:
                pct = sig / total * 100
                axes[plot_idx].text(i, max(sig_values) * 1.1, f'{pct:.1f}%', 
                                   ha='center', va='bottom')
    
    # Remove any unused subplots
    for i in range(plot_idx + 1, len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_enrichment_analysis.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def create_summary_plots(corrected_results, category_analysis, output_prefix):
    """Create additional summary plots"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Top enriched proteins across all conditions
    if 'differential_interactions' in corrected_results:
        df = corrected_results['differential_interactions']
        
        # Calculate mean enrichment per protein across conditions
        protein_means = df.groupby('protein')['log2_enrichment'].agg(['mean', 'std']).reset_index()
        protein_means = protein_means.sort_values('mean', ascending=False)
        
        top_proteins = protein_means.head(15)  # Top 15
        
        axes[0].barh(range(len(top_proteins)), top_proteins['mean'])
        axes[0].set_yticks(range(len(top_proteins)))
        axes[0].set_yticklabels(top_proteins['protein'])
        axes[0].set_xlabel('Mean Log2 Enrichment')
        axes[0].set_title('Top Enriched Architectural Proteins')
        axes[0].axvline(x=0, color='black', linestyle='--', alpha=0.3)
        
        # Add error bars if we have standard deviation
        if 'std' in top_proteins.columns:
            axes[0].errorbar(top_proteins['mean'], range(len(top_proteins)), 
                           xerr=top_proteins['std'], fmt='none', color='black', alpha=0.5)
    
    # Plot 2: Category comparison across strains
    if len(category_analysis) > 0:
        # Focus on differential interactions
        cat_diff = category_analysis[category_analysis['feature_type'] == 'differential_interactions']
        
        if len(cat_diff) > 0:
            # Plot mean enrichment by category and condition
            pivot_cat = cat_diff.pivot(index='condition', columns='category', values='mean_log2_enrichment')
            
            if not pivot_cat.empty:
                pivot_cat.plot(kind='bar', ax=axes[1])
                axes[1].set_title('Enrichment by Protein Category')
                axes[1].set_ylabel('Mean Log2 Enrichment')
                axes[1].set_xlabel('Wolbachia Strain')
                axes[1].legend(title='Protein Category')
                axes[1].axhline(y=0, color='black', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_summary_plots.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze architectural protein enrichment using existing diffHic results')
    parser.add_argument('--interactions', required=True, help='Combined differential interactions CSV file')
    parser.add_argument('--null_model', required=True, help='Null model results CSV file')
    parser.add_argument('--chip_dir', required=True, help='Directory with ChIP-seq peak files')
    parser.add_argument('--genome', required=True, help='Genome file for shuffling')
    parser.add_argument('--window_size', type=int, default=5000,
                       help='Window size around features')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significant interactions')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    parser.add_argument('--scan_only', action='store_true', 
                       help='Only scan ChIP directory and exit')
    
    args = parser.parse_args()
    
    # Scan ChIP directory first
    found_matches = scan_chip_directory(args.chip_dir)
    
    if args.scan_only:
        print("\nScan complete. Use the information above to verify file mappings.")
        return
    
    # Load null model
    null_stats = load_null_model(args.null_model)
    
    # Load ChIP-seq data
    chip_data = load_chip_peaks(args.chip_dir)
    
    if not chip_data:
        print("Error: No ChIP-seq data loaded")
        print("Run with --scan_only to see available files and their potential mappings")
        return
    
    print(f"\nSuccessfully loaded ChIP data for {len(chip_data)} proteins:")
    for protein in chip_data.keys():
        print(f"  - {protein}")
    
    # Load differential interactions
    features, interactions_df = load_differential_interactions(
        args.interactions, args.fdr_threshold
    )
    
    # Calculate enrichment using existing null model
    results = compare_enrichment_across_conditions(
        features, chip_data, args.genome, null_stats, args.window_size
    )
    
    # Apply FDR correction
    corrected_results = apply_fdr_correction(results)
    
    # Analyze by protein category
    category_analysis = analyze_protein_categories(corrected_results)
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create visualizations
    create_enrichment_heatmaps(corrected_results, args.output_prefix)
    create_summary_plots(corrected_results, category_analysis, args.output_prefix)
    
    # Save results
    for feature_type, df in corrected_results.items():
        df.to_csv(f"{args.output_prefix}_{feature_type}_enrichment.tsv", 
                 sep='\t', index=False)
    
    if len(category_analysis) > 0:
        category_analysis.to_csv(f"{args.output_prefix}_category_analysis.tsv",
                                sep='\t', index=False)
    
    # Create comprehensive summary
    summary = []
    total_significant = 0
    
    for feature_type, df in corrected_results.items():
        if 'significant' in df.columns and len(df) > 0:
            n_sig = df['significant'].sum()
            total_significant += n_sig
            
            summary.append({
                'feature_type': feature_type,
                'n_proteins_tested': len(df),
                'n_significant': n_sig,
                'percent_significant': n_sig / len(df) * 100 if len(df) > 0 else 0,
                'mean_log2_enrichment': df['log2_enrichment'].mean(),
                'median_log2_enrichment': df['log2_enrichment'].median()
            })
    
    if summary:
        summary_df = pd.DataFrame(summary)
        summary_df.to_csv(f"{args.output_prefix}_enrichment_summary.tsv",
                         sep='\t', index=False)
    
    # Print summary to console
    print(f"\n=== ARCHITECTURAL PROTEIN ENRICHMENT ANALYSIS SUMMARY ===")
    print(f"Total significant enrichments found: {total_significant}")
    print(f"ChIP-seq datasets analyzed: {len(chip_data)}")
    print(f"Differential interaction sites analyzed: {sum(len(f) for f in features['interactions'].values()) if features['interactions'] else 0}")
    
    if 'differential_interactions' in corrected_results:
        df = corrected_results['differential_interactions']
        top_proteins = df[df['significant']].nlargest(5, 'log2_enrichment')
        if len(top_proteins) > 0:
            print(f"\nTop 5 enriched proteins:")
            for _, row in top_proteins.iterrows():
                print(f"  {row['protein']} ({row['condition']}): {row['log2_enrichment']:.2f} log2 enrichment (FDR={row['fdr']:.3f})")
    
    print(f"\nResults saved to {args.output_prefix}_*")

if __name__ == '__main__':
    main()