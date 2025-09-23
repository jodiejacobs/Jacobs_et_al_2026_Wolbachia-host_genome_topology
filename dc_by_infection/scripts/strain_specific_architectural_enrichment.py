#!/usr/bin/env python3
"""
This is a working version as of 6/16/25
Direct strain-specific analysis script - save and run directly
"""

import pandas as pd
import pybedtools
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gzip

def load_bed_file_safe(file_path):
    """Safely load a BED file, skipping headers and invalid lines"""
    try:
        bed_lines = []
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments, headers, empty lines
                if not line or line.startswith('#') or line.startswith('track'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                
                try:
                    chrom = parts[0].replace('chr', '')
                    start = int(parts[1])
                    end = int(parts[2])
                    
                    # Basic validation
                    if start >= 0 and end > start and chrom in ['2L', '2R', '3L', '3R', '4', 'X', 'Y']:
                        bed_lines.append([chrom, start, end])
                except ValueError:
                    continue
        
        if bed_lines:
            bed_df = pd.DataFrame(bed_lines, columns=['chrom', 'start', 'end'])
            bed_df = bed_df.drop_duplicates()
            print(f"      Loaded {len(bed_df)} intervals from {file_path.name}")
            return pybedtools.BedTool.from_dataframe(bed_df)
        else:
            print(f"      No valid intervals in {file_path.name}")
            return None
            
    except Exception as e:
        print(f"      Error loading {file_path.name}: {e}")
        return None

def load_chip_peaks_working(chip_dir):
    """Load ChIP-seq peaks using only files known to work"""
    
    # Files that we know have actual data (not just headers)
    WORKING_FILES = {
        'CTCF': 'GSM762842_dCTCF_20HE_0hrs_ChIPSeq.bed.gz',
        'Su(Hw)': 'GSM762839_Su_Hw_20HE_0hrs_ChIPSeq.bed.gz', 
        'BEAF-32': 'GSM762845_BEAF-32_20HE_0hrs_ChIPSeq.bed.gz',
        'CP190': 'GSM762836_CP190_20HE_0hrs_ChIPSeq.bed.gz',
        'Ibf1': 'GSM2133766_Kc_Ibf1.bed.gz',
        'Ibf2': 'GSM2133767_Kc_Ibf2.bed.gz',
        'ZIPIC': 'GSM2133769_Kc_ZIPIC.bed.gz',
        'Pita': 'GSM2133768_Kc_Pita.bed.gz',
        'GAF': 'GSM2133762_Kc_GAF.bed.gz',
        'L3mbt': 'GSM892323_L3mbt.bed.gz',
        'Chromator': 'GSM1318357_Chromator.bed.gz',
        'Rad21': 'GSM1318352_Rad21.bed.gz'
    }
    
    chip_data = {}
    chip_dir = Path(chip_dir)
    
    print("Loading ChIP-seq data (reliable files only)...")
    
    for protein, filename in WORKING_FILES.items():
        print(f"  Loading {protein}...")
        chip_file = chip_dir / filename
        
        if chip_file.exists():
            bed_tool = load_bed_file_safe(chip_file)
            if bed_tool is not None and bed_tool.count() > 0:
                chip_data[protein] = bed_tool
            else:
                print(f"    Failed to load {protein}")
        else:
            print(f"    File not found: {filename}")
    
    return chip_data

def load_strain_interactions(interactions_file):
    """Load and separate interactions by strain"""
    print(f"Loading interactions from {interactions_file}")
    
    interactions = pd.read_csv(interactions_file)
    print(f"Loaded {len(interactions)} total interactions")
    
    # Find strain column
    strain_col = None
    for col in ['infection', 'filename']:
        if col in interactions.columns:
            strain_col = col
            break
    
    if strain_col is None:
        print("Error: No strain identifier column found")
        return None
    
    print(f"Using column '{strain_col}' for strain identification")
    
    # Separate by strain
    strain_data = {}
    
    for strain in ['wMel', 'wRi', 'wWil']:
        print(f"\nProcessing {strain}...")
        mask = interactions[strain_col].str.contains(strain, na=False)
        strain_df = interactions[mask].copy()
        
        print(f"  Found {len(strain_df)} interactions for {strain}")
        
        if len(strain_df) > 0:
            # Filter for significant interactions
            if 'FDR' in strain_df.columns:
                pre_fdr = len(strain_df)
                strain_df = strain_df[strain_df['FDR'] < 0.01].copy()
                print(f"  After FDR < 0.01: {len(strain_df)}/{pre_fdr}")
            
            if 'logFC' in strain_df.columns:
                pre_fc = len(strain_df)
                strain_df = strain_df[abs(strain_df['logFC']) > 1].copy()
                print(f"  After |logFC| > 1: {len(strain_df)}/{pre_fc}")
            
            if len(strain_df) > 0:
                strain_data[strain] = strain_df
                print(f"  Final: {len(strain_df)} significant interactions for {strain}")
            else:
                print(f"  No significant interactions remaining for {strain}")
        else:
            print(f"  No interactions found for {strain}")
    
    return strain_data

def create_bed_from_interactions(df):
    """Create BED regions from interaction endpoints"""
    regions = []
    for _, row in df.iterrows():
        regions.append([row['chr1'], row['start1'], row['end1']])
        regions.append([row['chr2'], row['start2'], row['end2']])
    
    bed_df = pd.DataFrame(regions, columns=['chrom', 'start', 'end'])
    bed_df = bed_df.drop_duplicates().reset_index(drop=True)
    
    print(f"    Created {len(bed_df)} unique regions from {len(df)} interactions")
    return pybedtools.BedTool.from_dataframe(bed_df)

def calculate_enrichment(query_bed, reference_bed, genome_file, n_perm=30):
    """Calculate enrichment with permutation test"""
    try:
        # Observed overlaps
        observed = query_bed.intersect(reference_bed, u=True).count()
        
        # Permutation test
        null_counts = []
        for i in range(n_perm):
            try:
                shuffled = query_bed.shuffle(g=genome_file, chrom=True, seed=i)
                null_count = shuffled.intersect(reference_bed, u=True).count()
                null_counts.append(null_count)
            except:
                continue
        
        if len(null_counts) < 5:
            return {
                'observed': observed,
                'expected': float(observed),
                'enrichment': 1.0,
                'log2_enrichment': 0.0,
                'p_value': 1.0
            }
        
        # Calculate statistics
        null_counts = np.array(null_counts)
        expected = np.median(null_counts)
        
        if expected > 0:
            enrichment = observed / expected
            log2_enrichment = np.log2(enrichment)
        else:
            enrichment = 10.0 if observed > 0 else 1.0
            log2_enrichment = 3.32 if observed > 0 else 0.0
        
        # P-value calculation
        if observed >= expected:
            p_value = (np.sum(null_counts >= observed) + 1) / (len(null_counts) + 1)
        else:
            p_value = (np.sum(null_counts <= observed) + 1) / (len(null_counts) + 1)
        
        p_value = min(p_value * 2, 1.0)  # Two-tailed
        
        return {
            'observed': observed,
            'expected': float(expected),
            'enrichment': float(enrichment),
            'log2_enrichment': float(log2_enrichment),
            'p_value': float(p_value)
        }
        
    except Exception as e:
        print(f"      Error: {e}")
        return {
            'observed': 0,
            'expected': 0.0,
            'enrichment': 1.0,
            'log2_enrichment': 0.0,
            'p_value': 1.0
        }

def analyze_strain_enrichment(strain_interactions, chip_data, genome_file, window_size=5000):
    """Analyze enrichment for each strain separately"""
    results = []
    
    print("\nAnalyzing strain-specific enrichment...")
    
    for strain, interactions_df in strain_interactions.items():
        print(f"\n=== Analyzing {strain} ===")
        
        # Create BED regions
        int_features = create_bed_from_interactions(interactions_df)
        
        # Add windows
        try:
            int_windows = int_features.slop(b=window_size, g=genome_file)
        except Exception as e:
            print(f"    Error adding windows: {e}")
            continue
        
        # Test each protein
        for protein, peaks in chip_data.items():
            print(f"    Testing {protein}...")
            
            enrichment = calculate_enrichment(int_windows, peaks, genome_file)
            
            result = {
                'strain': strain,
                'protein': protein,
                'n_interactions': len(interactions_df),
                'observed': enrichment['observed'],
                'expected': enrichment['expected'],
                'enrichment': enrichment['enrichment'],
                'log2_enrichment': enrichment['log2_enrichment'],
                'p_value': enrichment['p_value']
            }
            
            results.append(result)
            
            print(f"      obs={enrichment['observed']}, exp={enrichment['expected']:.1f}, "
                  f"fold={enrichment['enrichment']:.2f}, p={enrichment['p_value']:.3f}")
    
    return pd.DataFrame(results)

def apply_fdr_and_create_summary(results_df, output_prefix):
    """Apply FDR correction and create summary"""
    if len(results_df) == 0:
        return results_df
    
    # Add categories
    dna_binding = ['CTCF', 'Su(Hw)', 'BEAF-32', 'ZIPIC', 'Ibf1', 'Ibf2', 'Pita']
    results_df['category'] = results_df['protein'].apply(
        lambda x: 'DNA_binding' if x in dna_binding else 'accessory'
    )
    
    # FDR correction
    _, global_fdr, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    results_df['fdr'] = global_fdr
    results_df['significant'] = results_df['fdr'] < 0.05
    
    # Save results
    results_df.to_csv(f"{output_prefix}_strain_results.tsv", sep='\t', index=False)
    
    # Create summary plot
    if len(results_df) > 0:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Heatmap of enrichments
        pivot_enrichment = results_df.pivot(index='protein', columns='strain', values='log2_enrichment')
        pivot_enrichment = pivot_enrichment.fillna(0)
        
        ax1 = axes[0, 0]
        sns.heatmap(pivot_enrichment, cmap='RdBu_r', center=0, annot=True, fmt='.2f', ax=ax1)
        ax1.set_title('Log2 Enrichment by Strain')
        
        # Significance heatmap
        pivot_sig = results_df.pivot(index='protein', columns='strain', values='significant')
        pivot_sig = pivot_sig.fillna(False).astype(int)
        
        ax2 = axes[0, 1]
        sns.heatmap(pivot_sig, cmap='Reds', annot=True, fmt='d', ax=ax2)
        ax2.set_title('Significance by Strain')
        
        # Bar plot
        ax3 = axes[1, 0]
        sig_counts = results_df.groupby('strain')['significant'].sum()
        sig_counts.plot(kind='bar', ax=ax3)
        ax3.set_title('Significant Proteins by Strain')
        ax3.tick_params(axis='x', rotation=45)
        
        # Category plot
        ax4 = axes[1, 1]
        category_sig = results_df.groupby(['strain', 'category'])['significant'].sum().unstack(fill_value=0)
        category_sig.plot(kind='bar', stacked=True, ax=ax4)
        ax4.set_title('Significant Proteins by Category')
        ax4.tick_params(axis='x', rotation=45)
        ax4.legend(title='Category')
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_strain_plots.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Plots saved to {output_prefix}_strain_plots.pdf")
    
    # Print summary
    total_tests = len(results_df)
    global_sig = results_df['significant'].sum()
    
    print(f"\n=== STRAIN-SPECIFIC SUMMARY ===")
    print(f"Total tests: {total_tests}")
    print(f"Significant: {global_sig} ({global_sig/total_tests*100:.1f}%)")
    
    print(f"\nBy strain:")
    for strain in results_df['strain'].unique():
        strain_data = results_df[results_df['strain'] == strain]
        n_sig = strain_data['significant'].sum()
        n_total = len(strain_data)
        print(f"  {strain}: {n_sig}/{n_total} ({n_sig/n_total*100:.1f}%)")
        
        if n_sig > 0:
            sig_proteins = strain_data[strain_data['significant']]['protein'].tolist()
            print(f"    Significant: {', '.join(sig_proteins)}")
    
    # Save text summary
    with open(f"{output_prefix}_strain_summary.txt", 'w') as f:
        f.write("STRAIN-SPECIFIC ARCHITECTURAL PROTEIN ENRICHMENT\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total tests performed: {total_tests}\n")
        f.write(f"Globally significant: {global_sig} ({global_sig/total_tests*100:.1f}%)\n\n")
        
        for strain in results_df['strain'].unique():
            strain_data = results_df[results_df['strain'] == strain]
            n_sig = strain_data['significant'].sum()
            n_total = len(strain_data)
            
            f.write(f"{strain}:\n")
            f.write(f"  Proteins tested: {n_total}\n")
            f.write(f"  Significant: {n_sig} ({n_sig/n_total*100:.1f}%)\n")
            
            if n_sig > 0:
                sig_proteins = strain_data[strain_data['significant']]['protein'].tolist()
                f.write(f"  Significant proteins: {', '.join(sig_proteins)}\n")
                
                # Show enrichment values
                for _, row in strain_data[strain_data['significant']].iterrows():
                    direction = "enriched" if row['log2_enrichment'] > 0 else "depleted"
                    f.write(f"    {row['protein']}: {row['enrichment']:.2f}x {direction} "
                           f"(p={row['p_value']:.3f}, FDR={row['fdr']:.3f})\n")
            f.write("\n")
    
    return results_df

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Strain-specific architectural protein enrichment analysis')
    parser.add_argument('--interactions', required=True, help='Combined differential interactions CSV file')
    parser.add_argument('--chip_dir', required=True, help='Directory with ChIP-seq peak files')
    parser.add_argument('--genome', required=True, help='Genome file for shuffling')
    parser.add_argument('--null_model', help='Null model file from diffHic (optional)')
    parser.add_argument('--window_size', type=int, default=5000, help='Window size around interactions')
    parser.add_argument('--fdr_threshold', type=float, default=0.05, help='FDR threshold for significance')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output_prefix).parent.mkdir(parents=True, exist_ok=True)
    
    # Load ChIP data
    chip_data = load_chip_peaks_working(args.chip_dir)
    
    if not chip_data:
        print("No ChIP data loaded")
        return 1
    
    print(f"\nSuccessfully loaded {len(chip_data)} proteins:")
    for protein, bed_tool in chip_data.items():
        print(f"  {protein}: {bed_tool.count()} peaks")
    
    # Load strain interactions
    strain_interactions = load_strain_interactions(args.interactions)
    
    if not strain_interactions:
        print("No strain interactions loaded")
        return 1
    
    print(f"\nLoaded interactions for {len(strain_interactions)} strains:")
    for strain, df in strain_interactions.items():
        print(f"  {strain}: {len(df)} interactions")
    
    # Analyze enrichment
    results = analyze_strain_enrichment(strain_interactions, chip_data, args.genome, args.window_size)
    
    if len(results) == 0:
        print("No results obtained")
        return 1
    
    # Apply FDR and create summary
    results = apply_fdr_and_create_summary(results, args.output_prefix)
    
    print(f"\nResults saved to {args.output_prefix}_*")
    return 0

if __name__ == '__main__':
    exit(main())
