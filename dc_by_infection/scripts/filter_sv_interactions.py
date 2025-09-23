#!/usr/bin/env python3
"""
This one works as of 6/17/25
Optimized version: Filter differential chromatin interactions that overlap with structural variants.
Following methods from Cubenas-Potts et al. 2017 and Ghavi-Helm et al. 2019.
"""

import pandas as pd
import numpy as np
import pybedtools
import argparse
import os
import tempfile
from pathlib import Path
import sys

def read_sv_vcf(vcf_file, quality_threshold=10):
    """Read and filter structural variants from VCF file"""
    print(f"Reading SVs from {vcf_file}")
    
    svs = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
                
            chrom = fields[0]
            pos = int(fields[1])
            
            # Handle quality score
            try:
                qual = float(fields[5]) if fields[5] != '.' else 0
            except ValueError:
                qual = 0
            
            # Skip low quality variants
            if qual < quality_threshold:
                continue
                
            # Extract SV info from INFO field
            info = fields[7]
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
            
            # Get SV type and end position
            sv_type = info_dict.get('SVTYPE', 'UNK')
            try:
                end = int(info_dict.get('END', pos + 1))
            except ValueError:
                end = pos + 1
            
            # Add padding around breakpoints (200kb as in methods)
            padding = 200000
            svs.append({
                'chrom': chrom,
                'start': max(0, pos - padding),
                'end': end + padding,
                'sv_type': sv_type,
                'quality': qual
            })
    
    sv_df = pd.DataFrame(svs)
    print(f"Found {len(sv_df)} high-quality SVs")
    return sv_df

def create_temp_bed_file(df, columns):
    """Create a temporary BED file from DataFrame"""
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
    
    for idx, row in df.iterrows():
        temp_file.write(f"{row[columns[0]]}\t{row[columns[1]]}\t{row[columns[2]]}\t{idx}\n")
    
    temp_file.close()
    return temp_file.name

def efficient_overlap_detection(interactions_df, sv_df, chunk_size=50000):
    """
    Use chunked processing for efficient overlap detection
    """
    print("Using efficient chunked overlap detection...")
    
    # Create temporary BED file for SVs
    sv_bed_file = create_temp_bed_file(sv_df, ['chrom', 'start', 'end'])
    sv_bed = pybedtools.BedTool(sv_bed_file)
    
    overlapping_indices = set()
    
    # Process interactions in chunks
    n_chunks = (len(interactions_df) // chunk_size) + 1
    print(f"Processing {len(interactions_df)} interactions in {n_chunks} chunks...")
    
    for chunk_idx in range(n_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min((chunk_idx + 1) * chunk_size, len(interactions_df))
        
        if start_idx >= len(interactions_df):
            break
            
        chunk = interactions_df.iloc[start_idx:end_idx]
        print(f"  Processing chunk {chunk_idx + 1}/{n_chunks} ({len(chunk)} interactions)...")
        
        # Process anchor1
        anchor1_bed_file = create_temp_bed_file(chunk, ['chr1', 'start1', 'end1'])
        try:
            anchor1_bed = pybedtools.BedTool(anchor1_bed_file)
            anchor1_overlaps = anchor1_bed.intersect(sv_bed, wa=True)
            
            for overlap in anchor1_overlaps:
                original_idx = int(overlap.name)
                actual_idx = chunk.index[original_idx]  # Map back to original DataFrame index
                overlapping_indices.add(actual_idx)
                
        except Exception as e:
            print(f"Error processing anchor1 in chunk {chunk_idx}: {e}")
        finally:
            os.unlink(anchor1_bed_file)
        
        # Process anchor2
        anchor2_bed_file = create_temp_bed_file(chunk, ['chr2', 'start2', 'end2'])
        try:
            anchor2_bed = pybedtools.BedTool(anchor2_bed_file)
            anchor2_overlaps = anchor2_bed.intersect(sv_bed, wa=True)
            
            for overlap in anchor2_overlaps:
                original_idx = int(overlap.name)
                actual_idx = chunk.index[original_idx]  # Map back to original DataFrame index
                overlapping_indices.add(actual_idx)
                
        except Exception as e:
            print(f"Error processing anchor2 in chunk {chunk_idx}: {e}")
        finally:
            os.unlink(anchor2_bed_file)
    
    # Clean up
    os.unlink(sv_bed_file)
    
    print(f"Found {len(overlapping_indices)} interactions overlapping with SVs")
    return overlapping_indices

def pandas_overlap_detection_optimized(interactions_df, sv_df):
    """
    Optimized pandas-based overlap detection using vectorized operations
    """
    print("Using optimized pandas-based overlap detection...")
    
    overlapping_indices = set()
    
    # Create chromosome-indexed dictionaries for faster lookup
    sv_by_chr = {}
    for chrom in sv_df['chrom'].unique():
        chr_svs = sv_df[sv_df['chrom'] == chrom].copy()
        # Sort by start position for faster searching
        chr_svs = chr_svs.sort_values('start').reset_index(drop=True)
        sv_by_chr[chrom] = chr_svs
    
    print(f"Created SV index for {len(sv_by_chr)} chromosomes")
    
    # Process interactions by chromosome for efficiency
    for chrom in interactions_df['chr1'].unique():
        if chrom not in sv_by_chr:
            continue
            
        chr_interactions = interactions_df[
            (interactions_df['chr1'] == chrom) | (interactions_df['chr2'] == chrom)
        ]
        
        if len(chr_interactions) == 0:
            continue
            
        print(f"  Processing chromosome {chrom}: {len(chr_interactions)} interactions")
        chr_svs = sv_by_chr[chrom]
        
        # Check anchor1 overlaps for this chromosome
        anchor1_chr = chr_interactions[chr_interactions['chr1'] == chrom]
        if len(anchor1_chr) > 0:
            for sv_idx, sv in chr_svs.iterrows():
                # Vectorized overlap check
                overlaps = (
                    (anchor1_chr['start1'] < sv['end']) & 
                    (anchor1_chr['end1'] > sv['start'])
                )
                overlapping_indices.update(anchor1_chr[overlaps].index)
        
        # Check anchor2 overlaps for this chromosome  
        anchor2_chr = chr_interactions[chr_interactions['chr2'] == chrom]
        if len(anchor2_chr) > 0:
            for sv_idx, sv in chr_svs.iterrows():
                # Vectorized overlap check
                overlaps = (
                    (anchor2_chr['start2'] < sv['end']) & 
                    (anchor2_chr['end2'] > sv['start'])
                )
                overlapping_indices.update(anchor2_chr[overlaps].index)
    
    print(f"Found {len(overlapping_indices)} interactions overlapping with SVs")
    return overlapping_indices

def filter_interactions_by_sv(interactions_file, sv_df, min_distance=4000):
    """
    Filter out interactions that might be artifacts due to SVs.
    """
    print(f"Reading interactions from {interactions_file}")
    
    # Read interactions in chunks if file is large
    file_size = os.path.getsize(interactions_file)
    print(f"Interactions file size: {file_size / (1024**2):.1f} MB")
    
    if file_size > 500 * 1024 * 1024:  # If larger than 500MB
        print("Large file detected, using chunked reading...")
        chunk_iter = pd.read_csv(interactions_file, chunksize=50000)
        interactions = pd.concat(chunk_iter, ignore_index=True)
    else:
        interactions = pd.read_csv(interactions_file)
    
    print(f"Loaded {len(interactions)} interactions")
    
    # Ensure coordinate columns are numeric
    coord_cols = ['start1', 'end1', 'start2', 'end2']
    for col in coord_cols:
        if col in interactions.columns:
            interactions[col] = pd.to_numeric(interactions[col], errors='coerce')
    
    # Remove rows with invalid coordinates
    before_clean = len(interactions)
    interactions = interactions.dropna(subset=coord_cols)
    print(f"Removed {before_clean - len(interactions)} interactions with invalid coordinates")
    
    # Try bedtools first, fall back to optimized pandas if it fails
    try:
        overlapping_indices = efficient_overlap_detection(interactions, sv_df)
    except Exception as e:
        print(f"Bedtools approach failed: {e}")
        print("Falling back to optimized pandas approach...")
        overlapping_indices = pandas_overlap_detection_optimized(interactions, sv_df)
    
    # Mark overlapping interactions
    interactions['overlaps_sv'] = False
    interactions.loc[interactions.index.isin(overlapping_indices), 'overlaps_sv'] = True
    
    # Calculate distances for cis interactions
    print("Calculating interaction distances...")
    cis_mask = interactions['chr1'] == interactions['chr2']
    interactions['distance'] = 0
    interactions.loc[cis_mask, 'distance'] = abs(
        interactions.loc[cis_mask, 'start2'] - interactions.loc[cis_mask, 'start1']
    )
    
    # Apply filters
    print("Applying filters...")
    filtered = interactions[
        (~interactions['overlaps_sv']) & 
        ((~cis_mask) | (interactions['distance'] >= min_distance))
    ].copy()
    
    # Print filtering statistics
    print(f"\nFiltering results:")
    print(f"  Original interactions: {len(interactions)}")
    print(f"  Overlapping SVs: {sum(interactions['overlaps_sv'])}")
    print(f"  Short-range cis interactions: {sum(cis_mask & (interactions['distance'] < min_distance))}")
    print(f"  Final filtered interactions: {len(filtered)}")
    print(f"  Retention rate: {len(filtered)/len(interactions)*100:.1f}%")
    
    return filtered

def main():
    parser = argparse.ArgumentParser(description='Filter interactions by structural variants')
    parser.add_argument('--interactions', required=True, help='Differential interactions CSV file')
    parser.add_argument('--vcf', required=True, help='SVIM VCF file with breakpoints')
    parser.add_argument('--output', required=True, help='Output filtered interactions CSV')
    parser.add_argument('--quality', type=int, default=10, help='Minimum SV quality score')
    parser.add_argument('--min_distance', type=int, default=4000, help='Minimum cis interaction distance')
    parser.add_argument('--chunk_size', type=int, default=50000, help='Chunk size for processing')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.interactions):
        print(f"Error: Interactions file not found: {args.interactions}")
        sys.exit(1)
    
    if not os.path.exists(args.vcf):
        print(f"Error: VCF file not found: {args.vcf}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        # Read structural variants
        sv_df = read_sv_vcf(args.vcf, args.quality)
        
        if len(sv_df) == 0:
            print("Warning: No structural variants found with specified quality threshold")
            # Just copy the input file if no SVs to filter
            interactions = pd.read_csv(args.interactions)
            interactions.to_csv(args.output, index=False)
            print(f"No filtering applied. Output saved to {args.output}")
            return
        
        # Filter interactions
        filtered = filter_interactions_by_sv(args.interactions, sv_df, args.min_distance)
        
        # Save results
        filtered.to_csv(args.output, index=False)
        print(f"\nFiltered interactions saved to {args.output}")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
