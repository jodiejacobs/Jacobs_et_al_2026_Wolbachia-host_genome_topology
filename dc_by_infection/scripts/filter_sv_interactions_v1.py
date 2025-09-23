#!/usr/bin/env python3
"""
Filter differential chromatin interactions that overlap with structural variants.
Following methods from Cubenas-Potts et al. 2017 and Ghavi-Helm et al. 2019.
"""

import pandas as pd
import pybedtools
import argparse
import os
from pathlib import Path

def read_sv_vcf(vcf_file, quality_threshold=10):
    """Read and filter structural variants from VCF file"""
    print(f"Reading SVs from {vcf_file}")
    
    svs = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            qual = float(fields[5]) if fields[5] != '.' else 0
            
            # Skip low quality variants
            if qual < quality_threshold:
                continue
                
            # Extract SV info from INFO field
            info = fields[7]
            info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)
            
            # Get SV type and end position
            sv_type = info_dict.get('SVTYPE', 'UNK')
            end = int(info_dict.get('END', pos + 1))
            
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

def filter_interactions_by_sv(interactions_file, sv_df, min_distance=4000):
    """
    Filter out interactions that might be artifacts due to SVs.
    Following Cubenas-Potts et al. 2017 method.
    """
    print(f"Reading interactions from {interactions_file}")
    interactions = pd.read_csv(interactions_file)
    
    # Add index to track original positions
    interactions['original_index'] = interactions.index
    
    # Convert to BED format for both anchors with proper names
    anchor1_bed = interactions[['chr1', 'start1', 'end1', 'original_index']].copy()
    anchor1_bed.columns = ['chrom', 'start', 'end', 'name']
    
    anchor2_bed = interactions[['chr2', 'start2', 'end2', 'original_index']].copy()
    anchor2_bed.columns = ['chrom', 'start', 'end', 'name']
    
    # Convert SVs to proper BED format
    sv_bed_df = sv_df[['chrom', 'start', 'end']].copy()
    
    # Convert DataFrames to BedTool objects
    sv_bed = pybedtools.BedTool.from_dataframe(sv_bed_df)
    anchor1_bt = pybedtools.BedTool.from_dataframe(anchor1_bed)
    anchor2_bt = pybedtools.BedTool.from_dataframe(anchor2_bed)
    
    # Find overlaps with SVs - use a simpler approach
    print("Finding overlaps between anchor1 and SVs...")
    anchor1_overlaps = anchor1_bt.intersect(sv_bed, wa=True)
    anchor1_overlap_indices = set()
    for interval in anchor1_overlaps:
        try:
            idx = int(interval.name)
            anchor1_overlap_indices.add(idx)
        except (ValueError, AttributeError):
            # Skip if we can't parse the index
            continue
    
    print("Finding overlaps between anchor2 and SVs...")
    anchor2_overlaps = anchor2_bt.intersect(sv_bed, wa=True)
    anchor2_overlap_indices = set()
    for interval in anchor2_overlaps:
        try:
            idx = int(interval.name)
            anchor2_overlap_indices.add(idx)
        except (ValueError, AttributeError):
            # Skip if we can't parse the index
            continue
    
    # Alternative approach using pandas for overlap detection
    print("Using pandas-based overlap detection as backup...")
    all_overlapping_indices = set()
    
    # For each interaction, check if either anchor overlaps with any SV
    for idx, row in interactions.iterrows():
        # Check anchor1
        anchor1_overlaps_sv = False
        for _, sv in sv_df.iterrows():
            if (row['chr1'] == sv['chrom'] and 
                not (row['end1'] < sv['start'] or row['start1'] > sv['end'])):
                anchor1_overlaps_sv = True
                break
        
        # Check anchor2
        anchor2_overlaps_sv = False
        for _, sv in sv_df.iterrows():
            if (row['chr2'] == sv['chrom'] and 
                not (row['end2'] < sv['start'] or row['start2'] > sv['end'])):
                anchor2_overlaps_sv = True
                break
        
        if anchor1_overlaps_sv or anchor2_overlaps_sv:
            all_overlapping_indices.add(idx)
    
    # Combine results from both methods
    all_overlapping = anchor1_overlap_indices | anchor2_overlap_indices | all_overlapping_indices
    
    # Mark interactions overlapping SVs
    interactions['overlaps_sv'] = interactions['original_index'].isin(all_overlapping)
    
    # Also filter by minimum distance for cis interactions
    cis_mask = interactions['chr1'] == interactions['chr2']
    interactions['distance'] = 0
    interactions.loc[cis_mask, 'distance'] = abs(
        interactions.loc[cis_mask, 'start2'] - interactions.loc[cis_mask, 'start1']
    )
    
    # Apply filters
    filtered = interactions[
        (~interactions['overlaps_sv']) & 
        ((~cis_mask) | (interactions['distance'] >= min_distance))
    ].copy()
    
    # Remove the temporary index column
    if 'original_index' in filtered.columns:
        filtered = filtered.drop('original_index', axis=1)
    
    print(f"Filtered {len(interactions)} to {len(filtered)} interactions")
    print(f"Removed {sum(interactions['overlaps_sv'])} interactions overlapping SVs")
    print(f"Removed {sum(cis_mask & (interactions['distance'] < min_distance))} short-range cis interactions")
    
    return filtered

def main():
    parser = argparse.ArgumentParser(description='Filter interactions by structural variants')
    parser.add_argument('--interactions', required=True, help='Differential interactions CSV file')
    parser.add_argument('--vcf', required=True, help='SVIM VCF file with breakpoints')
    parser.add_argument('--output', required=True, help='Output filtered interactions CSV')
    parser.add_argument('--quality', type=int, default=10, help='Minimum SV quality score')
    parser.add_argument('--min_distance', type=int, default=4000, help='Minimum cis interaction distance')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Read structural variants
    sv_df = read_sv_vcf(args.vcf, args.quality)
    
    # Filter interactions
    filtered = filter_interactions_by_sv(args.interactions, sv_df, args.min_distance)
    
    # Save results
    filtered.to_csv(args.output, index=False)
    print(f"Saved filtered interactions to {args.output}")
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"- Total SVs considered: {len(sv_df)}")
    print(f"- Final filtered interactions: {len(filtered)}")
    if 'overlaps_sv' in filtered.columns:
        print(f"- SV overlap filtering removed interactions")
    if 'distance' in filtered.columns:
        cis_filtered = filtered[filtered['chr1'] == filtered['chr2']]
        if len(cis_filtered) > 0:
            print(f"- Minimum cis distance: {cis_filtered['distance'].min()}")
            print(f"- Maximum cis distance: {cis_filtered['distance'].max()}")

if __name__ == '__main__':
    main()