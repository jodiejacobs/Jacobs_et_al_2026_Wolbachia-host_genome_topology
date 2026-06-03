#!/usr/bin/env python3
"""
Create mapping for ChIP-seq files and updated architectural enrichment script
"""

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
import os
import gzip

# Define architectural proteins with their actual file mappings
CHIP_FILE_MAPPING = {
    # DNA-binding proteins
    'CTCF': ['GSM1015410_Ct-CTCF_peaks.bed.gz', 'GSM1015411_Ab-CTCF_peaks.bed.gz', 
             'GSM762842_dCTCF_20HE_0hrs_ChIPSeq.bed.gz', 'GSM762843_dCTCF_20HE_3hrs_ChIPSeq.bed.gz'],
    'Su(Hw)': ['GSM1015406_CtSuHw_peaks.bed.gz', 'GSM1015407_AbSuHw_peaks.bed.gz',
               'GSM762839_Su_Hw_20HE_0hrs_ChIPSeq.bed.gz', 'GSM762840_Su_Hw_20HE_3hrs_ChIPSeq.bed.gz'],
    'BEAF-32': ['GSM762845_BEAF-32_20HE_0hrs_ChIPSeq.bed.gz', 'GSM762846_BEAF-32_20HE_3hrs_ChIPSeq.bed.gz'],
    'DREF': ['GSE39664_DREF_TOTAL.bed.gz', 'GSE39664_DREF_MITOTIC.bed.gz'],
    'TFIIIC': ['GSM1318349_dTFIIIC220_1.bed.gz', 'GSM1318350_dTFIIIC220_2.bed.gz', 'GSM1318351_dTFIIIC220_3.bed.gz'],
    'Z4': [],  # No Z4 files found
    'Elba': [],  # No Elba files found
    'ZIPIC': ['GSM1313421_ZIPIC_coverage_dm3.bed.gz', 'GSM2133769_Kc_ZIPIC.bed.gz'],
    'Ibf1': ['GSM2133766_Kc_Ibf1.bed.gz'],
    'Ibf2': ['GSM2133767_Kc_Ibf2.bed.gz'],
    'Pita': ['GSM1313420_Pita_coverage_dm3.bed.gz', 'GSM2133768_Kc_Pita.bed.gz'],
    
    # Accessory proteins
    'CP190': ['GSM1015404_CtCP190_peaks.bed.gz', 'GSM1015405_AbCP190_peaks.bed.gz',
              'GSM762836_CP190_20HE_0hrs_ChIPSeq.bed.gz', 'GSM762837_CP190_20HE_3hrs_ChIPSeq.bed.gz',
              'GSM1318359_CP190.bed.gz'],
    'Mod(mdg4)': ['GSM1015408_CtMod_peaks.bed.gz', 'GSM1015409_AbMod_peaks.bed.gz',
                  'GSM892321_Modmdg42.2.bed.gz', 'GSM892322_Modmdg4BTB.bed.gz'],
    'Rad21': ['GSM1318352_Rad21.bed.gz', 'GSM1363353_Rad21_dCTCF-RNAi.bed.gz'],
    'Cap-H2': ['GSM1318355_CAPH2_int.bed.gz', 'GSM1318356_CAPH2.bed.gz', 'GSM1363354_CAPH2_dCTCF-RNAi.bed.gz'],
    'Fs(1)h-L': ['GSM1032228_Fs1h_L_Kc_peaks.bed.gz', 'GSM1032229_Fs1h_LS_Kc_peaks.bed.gz'],
    'L3mbt': ['GSM892323_L3mbt.bed.gz'],
    'Chromator': ['GSM1318357_Chromator.bed.gz'],
    'GAF': ['GSM1318358_GAF.bed.gz', 'GSM2133762_Kc_GAF.bed.gz']
}

# Define protein categories
ARCHITECTURAL_PROTEINS = {
    'DNA_binding': ['CTCF', 'Su(Hw)', 'BEAF-32', 'DREF', 'TFIIIC', 'ZIPIC', 'Ibf1', 'Ibf2', 'Pita'],
    'accessory': ['CP190', 'Mod(mdg4)', 'Rad21', 'Cap-H2', 'Fs(1)h-L', 'L3mbt', 'Chromator', 'GAF']
}

def load_chip_peaks(chip_dir):
    """Load ChIP-seq peaks for architectural proteins using file mapping"""
    chip_data = {}
    chip_peak_counts = {}  # NEW: track peak counts per protein
    chip_dir = Path(chip_dir)
    
    for protein, file_list in CHIP_FILE_MAPPING.items():
        if not file_list:  # Skip proteins with no files
            print(f"No files available for {protein}")
            continue
            
        protein_beds = []
        
        for filename in file_list:
            chip_file = chip_dir / filename
            if chip_file.exists():
                print(f"Loading {protein} peaks from {filename}")
                try:
                    # Handle compressed files
                    if filename.endswith('.gz'):
                        # Read compressed BED file
                        with gzip.open(chip_file, 'rt') as f:
                            # Check if it's a proper BED file by reading first line
                            first_line = f.readline().strip()
                            if first_line and not first_line.startswith('#'):
                                # Reset file pointer
                                f.seek(0)
                                
                                # Read all lines and create BedTool
                                bed_lines = []
                                for line in f:
                                    line = line.strip()
                                    if line and not line.startswith('#') and not line.startswith('track'):
                                        parts = line.split('\t')
                                        if len(parts) >= 3:
                                            # Ensure we have at least chrom, start, end
                                            bed_lines.append([parts[0], parts[1], parts[2]])
                                
                                if bed_lines:
                                    bed_df = pd.DataFrame(bed_lines, columns=['chrom', 'start', 'end'])
                                    bed_df['start'] = pd.to_numeric(bed_df['start'])
                                    bed_df['end'] = pd.to_numeric(bed_df['end'])
                                    protein_beds.append(pybedtools.BedTool.from_dataframe(bed_df))
                    else:
                        # Uncompressed file
                        protein_beds.append(pybedtools.BedTool(str(chip_file)))
                        
                except Exception as e:
                    print(f"Warning: Could not load {filename} for {protein}: {e}")
                    continue
            else:
                print(f"File not found: {filename}")
        
        # Combine all BED intervals for this protein
        if protein_beds:
            if len(protein_beds) == 1:
                chip_data[protein] = protein_beds[0]
            else:
                # Merge multiple files for the same protein
                try:
                    combined = protein_beds[0]
                    for bed in protein_beds[1:]:
                        combined = combined.cat(bed, postmerge=False)
                    # Sort and merge overlapping intervals
                    chip_data[protein] = combined.sort().merge()
                except Exception as e:
                    print(f"Warning: Could not combine files for {protein}: {e}")
                    # Use the first file as fallback
                    chip_data[protein] = protein_beds[0]
            chip_peak_counts[protein] = chip_data[protein].count()  # NEW
        else:
            print(f"Warning: No valid ChIP data loaded for {protein}")
    
    return chip_data, chip_peak_counts  # MODIFIED: return counts too

def load_differential_interactions(interactions_file):
    """Load and parse differential interactions file"""
    print(f"Loading interactions from {interactions_file}")
    
    try:
        # Try to read the CSV file
        interactions = pd.read_csv(interactions_file)
        print(f"Loaded {len(interactions)} interactions")
        print(f"Columns: {interactions.columns.tolist()}")
        
        # Check what columns we have and standardize them
        required_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
        missing_cols = [col for col in required_cols if col not in interactions.columns]
        
        if missing_cols:
            print(f"Error: Missing required columns: {missing_cols}")
            print(f"Available columns: {interactions.columns.tolist()}")
            return None
        
        # Filter for significant interactions if FDR column exists
        if 'FDR' in interactions.columns:
            sig_interactions = interactions[interactions['FDR'] < 0.01].copy()
            print(f"Found {len(sig_interactions)} significant interactions (FDR < 0.01)")
        else:
            print("Warning: No FDR column found, using all interactions")
            sig_interactions = interactions.copy()
        
        # Further filter by logFC if available
        if 'logFC' in sig_interactions.columns:
            sig_interactions = sig_interactions[abs(sig_interactions['logFC']) > 1].copy()
            print(f"After logFC filter: {len(sig_interactions)} interactions")
        
        return sig_interactions
        
    except Exception as e:
        print(f"Error loading interactions file: {e}")
        return None

def create_bed_from_interactions(interactions_df):
    """Convert interaction endpoints to BED format"""
    bed_regions = []
    
    # Add anchor1 regions
    for _, row in interactions_df.iterrows():
        bed_regions.append([row['chr1'], row['start1'], row['end1']])
    
    # Add anchor2 regions  
    for _, row in interactions_df.iterrows():
        bed_regions.append([row['chr2'], row['start2'], row['end2']])
    
    # Remove duplicates
    bed_df = pd.DataFrame(bed_regions, columns=['chrom', 'start', 'end'])
    bed_df = bed_df.drop_duplicates().reset_index(drop=True)
    
    print(f"Created {len(bed_df)} unique genomic regions from interactions")
    return pybedtools.BedTool.from_dataframe(bed_df)

def calculate_enrichment_with_null(query_bed, reference_bed, genome_file, 
                                  n_permutations=1000, null_model=None):
    """Calculate enrichment with permutation test"""
    try:
        # Calculate observed overlap
        observed = query_bed.intersect(reference_bed, u=True).count()
        
        # Permutation test
        null_overlaps = []
        for i in range(min(n_permutations, 100)):  # Reduce permutations for speed
            try:
                shuffled = query_bed.shuffle(g=genome_file, chrom=True, seed=i)
                null_overlap = shuffled.intersect(reference_bed, u=True).count()
                null_overlaps.append(null_overlap)
            except Exception as e:
                if i < 10:  # Only print first few errors
                    print(f"Warning: Permutation {i} failed: {e}")
                continue
        
        if len(null_overlaps) < 10:
            print(f"Warning: Only {len(null_overlaps)} successful permutations")
            return {
                'observed': observed,
                'expected': observed,
                'enrichment': 1.0,
                'log2_enrichment': 0.0,
                'p_value': 1.0,
                'null_mean': observed,
                'null_std': 0.0,
                'null_overlaps': null_overlaps,  # NEW
            }
        
        # Calculate statistics
        null_overlaps_arr = np.array(null_overlaps)
        expected = np.median(null_overlaps_arr)
        
        if expected > 0:
            enrichment = observed / expected
            log2_enrichment = np.log2(enrichment)
        else:
            enrichment = np.inf if observed > 0 else 1
            log2_enrichment = np.inf if observed > 0 else 0
        
        # Calculate p-value (two-tailed)
        if observed > expected:
            p_value = (np.sum(null_overlaps_arr >= observed) + 1) / (len(null_overlaps_arr) + 1)
        else:
            p_value = (np.sum(null_overlaps_arr <= observed) + 1) / (len(null_overlaps_arr) + 1)
        
        p_value = min(p_value * 2, 1)  # Two-tailed
        
        return {
            'observed': observed,
            'expected': expected,
            'enrichment': enrichment,
            'log2_enrichment': log2_enrichment,
            'p_value': p_value,
            'null_mean': np.mean(null_overlaps_arr),
            'null_std': np.std(null_overlaps_arr),
            'null_overlaps': null_overlaps,  # NEW: preserve raw null values
        }
        
    except Exception as e:
        print(f"Error in enrichment calculation: {e}")
        return {
            'observed': 0,
            'expected': 0,
            'enrichment': 1.0,
            'log2_enrichment': 0.0,
            'p_value': 1.0,
            'null_mean': 0.0,
            'null_std': 0.0,
            'null_overlaps': [],  # NEW
        }

def compute_per_interval_overlap(int_windows, chip_data):
    """For each query interval, flag which insulator proteins overlap it.

    Returns a DataFrame with chrom/start/end, a 0/1 column per protein,
    plus 'n_insulators' and 'insulators_overlapping' (comma-separated names).
    """
    # Sort once so every per-protein intersect returns rows in the same order
    query_sorted = int_windows.sort()
    overlap_df = query_sorted.to_dataframe(names=['chrom', 'start', 'end'],
                                           dtype={'chrom': str})
    proteins = list(chip_data.keys())

    for protein, peaks in chip_data.items():
        # -c appends, per query interval, the number of overlapping peaks
        counted = query_sorted.intersect(peaks, c=True)
        cdf = counted.to_dataframe(names=['chrom', 'start', 'end', 'cnt'],
                                   dtype={'chrom': str})
        overlap_df[protein] = (cdf['cnt'].values > 0).astype(int)

    overlap_df['n_insulators'] = overlap_df[proteins].sum(axis=1)
    overlap_df['insulators_overlapping'] = overlap_df[proteins].apply(
        lambda row: ','.join([p for p in proteins if row[p] == 1]) or 'none',
        axis=1
    )
    return overlap_df, proteins


def analyze_differential_interactions(interactions_df, chip_data, genome_file, window_size=5000):
    """Analyze architectural protein enrichment at differential interaction sites"""
    print("\nAnalyzing differential interaction enrichment...")
    
    if interactions_df is None or len(interactions_df) == 0:
        print("No interactions to analyze")
        return pd.DataFrame(), {}, 0, pd.DataFrame()  # MODIFIED
    
    # Create BED regions from interactions
    int_features = create_bed_from_interactions(interactions_df)
    
    # Add windows around features
    try:
        int_windows = int_features.slop(b=window_size, g=genome_file)
    except Exception as e:
        print(f"Error adding windows to interactions: {e}")
        return pd.DataFrame(), {}, 0, pd.DataFrame()  # MODIFIED
    
    n_query_intervals = int_windows.count()  # NEW
    
    # NEW: which insulator proteins overlap each interval
    print("  Computing per-interval insulator overlaps...")
    try:
        per_interval_df, _ = compute_per_interval_overlap(int_windows, chip_data)
    except Exception as e:
        print(f"  Warning: per-interval overlap failed: {e}")
        per_interval_df = pd.DataFrame()
    
    results = []
    null_distributions = {}  # NEW: collect null values per protein
    
    # For each protein
    for protein, peaks in chip_data.items():
        print(f"  Analyzing {protein}...")
        
        try:
            enrichment = calculate_enrichment_with_null(
                int_windows, peaks, genome_file
            )
            
            null_distributions[protein] = enrichment.get('null_overlaps', [])  # NEW
            
            results.append({
                'feature_type': 'differential_interactions',
                'protein': protein,
                'observed': enrichment['observed'],
                'expected': enrichment['expected'],
                'enrichment': enrichment['enrichment'],
                'log2_enrichment': enrichment['log2_enrichment'],
                'p_value': enrichment['p_value'],
                'null_mean': enrichment['null_mean'],  # NEW
                'null_std': enrichment['null_std'],   # NEW
                'n_permutations': len(enrichment.get('null_overlaps', [])),  # NEW
            })
            
        except Exception as e:
            print(f"  Error analyzing {protein}: {e}")
            continue
    
    return pd.DataFrame(results), null_distributions, n_query_intervals, per_interval_df  # MODIFIED

def apply_fdr_correction(results_df):
    """Apply FDR correction to p-values"""
    if len(results_df) == 0 or 'p_value' not in results_df.columns:
        return results_df
    
    _, fdr_values, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    results_df['fdr'] = fdr_values
    results_df['significant'] = results_df['fdr'] < 0.05
    
    return results_df

def create_enrichment_plots(results_df, output_prefix):
    """Create visualizations for enrichment results"""
    if len(results_df) == 0:
        print("No results to plot")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Enrichment by protein
    ax = axes[0, 0]
    if 'log2_enrichment' in results_df.columns:
        # Sort by enrichment for better visualization
        sorted_df = results_df.sort_values('log2_enrichment', ascending=True)
        y_pos = np.arange(len(sorted_df))
        
        colors = ['red' if sig else 'gray' for sig in sorted_df.get('significant', [False]*len(sorted_df))]
        ax.barh(y_pos, sorted_df['log2_enrichment'], color=colors, alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(sorted_df['protein'])
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel('Log2 Enrichment')
        ax.set_title('Protein Enrichment at Differential Interactions')
    
    # Plot 2: Observed vs Expected
    ax = axes[0, 1]
    if 'observed' in results_df.columns and 'expected' in results_df.columns:
        ax.scatter(results_df['expected'], results_df['observed'], alpha=0.7)
        
        # Add diagonal line
        max_val = max(results_df['expected'].max(), results_df['observed'].max())
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
        
        ax.set_xlabel('Expected Overlaps')
        ax.set_ylabel('Observed Overlaps')
        ax.set_title('Observed vs Expected Overlaps')
        
        # Add protein labels for significant points
        if 'significant' in results_df.columns:
            sig_results = results_df[results_df['significant']]
            for _, row in sig_results.iterrows():
                ax.annotate(row['protein'], 
                          (row['expected'], row['observed']),
                          fontsize=8, alpha=0.8)
    
    # Plot 3: Category analysis
    ax = axes[1, 0]
    results_df['category'] = results_df['protein'].apply(
        lambda x: 'DNA_binding' if x in ARCHITECTURAL_PROTEINS['DNA_binding'] 
        else 'accessory' if x in ARCHITECTURAL_PROTEINS['accessory'] 
        else 'other'
    )
    
    if 'log2_enrichment' in results_df.columns:
        category_data = results_df.groupby('category').agg({
            'log2_enrichment': ['mean', 'std'],
            'significant': 'sum'
        }).round(3)
        
        categories = category_data.index
        means = category_data[('log2_enrichment', 'mean')]
        stds = category_data[('log2_enrichment', 'std')]
        
        ax.bar(categories, means, yerr=stds, alpha=0.7, capsize=5)
        ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        ax.set_ylabel('Mean Log2 Enrichment')
        ax.set_title('Enrichment by Protein Category')
        ax.tick_params(axis='x', rotation=45)
    
    # Plot 4: Volcano plot
    ax = axes[1, 1]
    if 'p_value' in results_df.columns and 'log2_enrichment' in results_df.columns:
        x = results_df['log2_enrichment']
        y = -np.log10(results_df['p_value'])
        
        # Color by significance
        colors = ['red' if sig else 'gray' for sig in results_df.get('significant', [False]*len(results_df))]
        ax.scatter(x, y, c=colors, alpha=0.7)
        
        # Add protein names for significant points
        if 'significant' in results_df.columns:
            sig_results = results_df[results_df['significant']]
            for _, row in sig_results.iterrows():
                ax.annotate(row['protein'], 
                          (row['log2_enrichment'], -np.log10(row['p_value'])),
                          fontsize=8, alpha=0.8)
        
        ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel('Log2 Enrichment')
        ax.set_ylabel('-Log10 P-value')
        ax.set_title('Enrichment Volcano Plot')
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_enrichment_analysis.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_prefix}_enrichment_analysis.pdf")


# ============================================================
# NEW: text-file output functions for raw values
# ============================================================

def write_overlap_values(results_df, null_distributions, chip_peak_counts,
                         n_query_intervals, output_prefix):
    """Write the raw values used in enrichment calculations to text files."""

    # File 1: per-protein observed/expected and null summary stats
    values_file = f"{output_prefix}_overlap_values.tsv"
    cols = ['protein', 'observed', 'expected', 'enrichment', 'log2_enrichment',
            'p_value', 'fdr', 'significant', 'null_mean', 'null_std',
            'n_permutations', 'n_chip_peaks', 'n_query_intervals']
    df_out = results_df.copy()
    df_out['n_chip_peaks'] = df_out['protein'].map(chip_peak_counts)
    df_out['n_query_intervals'] = n_query_intervals
    # keep only columns that exist
    cols = [c for c in cols if c in df_out.columns]
    df_out[cols].to_csv(values_file, sep='\t', index=False)
    print(f"Wrote per-protein overlap values to {values_file}")

    # File 2: full null distributions (one row per permutation per protein)
    null_file = f"{output_prefix}_null_distributions.tsv"
    with open(null_file, 'w') as f:
        f.write("protein\tpermutation\tnull_overlap\n")
        for protein, nulls in null_distributions.items():
            for i, n in enumerate(nulls):
                f.write(f"{protein}\t{i}\t{n}\n")
    print(f"Wrote full null distributions to {null_file}")


def write_chip_peak_counts(chip_peak_counts, output_prefix):
    """Write the number of merged peaks loaded per protein."""
    chip_file = f"{output_prefix}_chip_peak_counts.tsv"
    with open(chip_file, 'w') as f:
        f.write("protein\tcategory\tn_merged_peaks\n")
        for protein, n in chip_peak_counts.items():
            if protein in ARCHITECTURAL_PROTEINS['DNA_binding']:
                cat = 'DNA_binding'
            elif protein in ARCHITECTURAL_PROTEINS['accessory']:
                cat = 'accessory'
            else:
                cat = 'other'
            f.write(f"{protein}\t{cat}\t{n}\n")
    print(f"Wrote ChIP peak counts to {chip_file}")


def write_interactions_bed(interactions_df, output_prefix):
    """Write the unique anchor BED regions used as the query."""
    bed_file = f"{output_prefix}_interaction_anchors.bed"
    rows = []
    for _, r in interactions_df.iterrows():
        rows.append([r['chr1'], int(r['start1']), int(r['end1'])])
        rows.append([r['chr2'], int(r['start2']), int(r['end2'])])
    bed_df = pd.DataFrame(rows, columns=['chrom', 'start', 'end']).drop_duplicates()
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)
    print(f"Wrote {len(bed_df)} unique anchor regions to {bed_file}")


def write_insulator_overlap(per_interval_df, output_prefix):
    """Write, for each query interval, which insulator proteins overlap it."""
    if per_interval_df is None or len(per_interval_df) == 0:
        print("No per-interval overlap data to write")
        return
    overlap_file = f"{output_prefix}_insulator_overlap.tsv"
    per_interval_df.to_csv(overlap_file, sep='\t', index=False)
    print(f"Wrote per-interval insulator overlaps to {overlap_file}")

    # Also write a compact version: only intervals that overlap >=1 insulator,
    # with just coords + the list of proteins
    hits = per_interval_df[per_interval_df['n_insulators'] > 0]
    compact_file = f"{output_prefix}_insulator_overlap_hits.tsv"
    hits[['chrom', 'start', 'end', 'n_insulators', 'insulators_overlapping']].to_csv(
        compact_file, sep='\t', index=False)
    print(f"Wrote {len(hits)} overlapping intervals to {compact_file}")


# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Analyze architectural protein enrichment at differential interactions')
    parser.add_argument('--interactions', required=True, help='Differential interactions CSV file')
    parser.add_argument('--chip_dir', required=True, help='Directory with ChIP-seq peak files')
    parser.add_argument('--genome', required=True, help='Genome file for shuffling')
    parser.add_argument('--null_model', help='Null model file from diffHic (optional)')
    parser.add_argument('--window_size', type=int, default=5000,
                       help='Window size around interactions')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold for significance')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load ChIP-seq data
    print("Loading ChIP-seq data...")
    chip_data, chip_peak_counts = load_chip_peaks(args.chip_dir)  # MODIFIED
    
    if not chip_data:
        print("Error: No ChIP-seq data loaded")
        return 1
    
    print(f"Loaded ChIP data for {len(chip_data)} proteins: {list(chip_data.keys())}")
    
    # NEW: dump ChIP peak counts
    write_chip_peak_counts(chip_peak_counts, args.output_prefix)
    
    # Load differential interactions
    interactions_df = load_differential_interactions(args.interactions)
    
    if interactions_df is None:
        print("Error: Could not load interactions")
        return 1
    
    # NEW: dump query anchor BED
    write_interactions_bed(interactions_df, args.output_prefix)
    
    # Analyze enrichment
    print("\nAnalyzing enrichment...")
    results, null_distributions, n_query_intervals, per_interval_df = analyze_differential_interactions(  # MODIFIED
        interactions_df, chip_data, args.genome, args.window_size
    )
    
    if len(results) == 0:
        print("No results obtained")
        return 1
    
    # Apply FDR correction
    print("Applying FDR correction...")
    corrected_results = apply_fdr_correction(results)
    
    # Create visualizations
    print("Creating visualizations...")
    create_enrichment_plots(corrected_results, args.output_prefix)
    
    # Save results
    print("Saving results...")
    corrected_results.to_csv(f"{args.output_prefix}_enrichment.tsv", 
                           sep='\t', index=False)
    
    # NEW: dump raw overlap values and null distributions
    write_overlap_values(corrected_results, null_distributions, chip_peak_counts,
                         n_query_intervals, args.output_prefix)
    
    # NEW: dump which insulator proteins overlap each interval
    write_insulator_overlap(per_interval_df, args.output_prefix)
    
    # Summary statistics
    n_significant = corrected_results['significant'].sum() if 'significant' in corrected_results.columns else 0
    summary = {
        'n_proteins_tested': len(corrected_results),
        'n_significant': n_significant,
        'percent_significant': n_significant / len(corrected_results) * 100 if len(corrected_results) > 0 else 0,
        'proteins_tested': list(corrected_results['protein']),
        'significant_proteins': list(corrected_results[corrected_results['significant']]['protein']) if 'significant' in corrected_results.columns else []
    }
    
    # Save summary
    with open(f"{args.output_prefix}_summary.txt", 'w') as f:
        f.write("ARCHITECTURAL PROTEIN ENRICHMENT ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Proteins tested: {summary['n_proteins_tested']}\n")
        f.write(f"Significant enrichments: {summary['n_significant']}\n")
        f.write(f"Percentage significant: {summary['percent_significant']:.1f}%\n\n")
        f.write(f"Proteins analyzed: {', '.join(summary['proteins_tested'])}\n\n")
        if summary['significant_proteins']:
            f.write(f"Significantly enriched proteins: {', '.join(summary['significant_proteins'])}\n")
        else:
            f.write("No significantly enriched proteins found.\n")
    
    print(f"\nAnalysis complete. Results saved to {args.output_prefix}_*")
    print(f"Tested {len(corrected_results)} proteins")
    print(f"Found {n_significant} significant enrichments")
    
    return 0

if __name__ == '__main__':
    exit(main())