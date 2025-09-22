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
    
    # Convert to BED format for both anchors
    anchor1_bed = interactions[['chr1', 'start1', 'end1']].copy()
    anchor2_bed = interactions[['chr2', 'start2', 'end2']].copy()
    
    # Convert DataFrames to BedTool objects
    sv_bed = pybedtools.BedTool.from_dataframe(sv_df[['chrom', 'start', 'end']])
    
    # Find overlaps with SVs
    anchor1_bt = pybedtools.BedTool.from_dataframe(anchor1_bed)
    anchor2_bt = pybedtools.BedTool.from_dataframe(anchor2_bed)
    
    # Get overlapping indices
    anchor1_overlaps = set()
    for overlap in anchor1_bt.intersect(sv_bed, wa=True):
        # Extract original index from the interaction
        idx = int(overlap.name) if hasattr(overlap, 'name') else len(anchor1_overlaps)
        anchor1_overlaps.add(idx)
    
    anchor2_overlaps = set()
    for overlap in anchor2_bt.intersect(sv_bed, wa=True):
        idx = int(overlap.name) if hasattr(overlap, 'name') else len(anchor2_overlaps)
        anchor2_overlaps.add(idx)
    
    # Mark interactions overlapping SVs
    interactions['overlaps_sv'] = False
    interactions.loc[interactions.index.isin(anchor1_overlaps | anchor2_overlaps), 'overlaps_sv'] = True
    
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
    
    # Read structural variants
    sv_df = read_sv_vcf(args.vcf, args.quality)
    
    # Filter interactions
    filtered = filter_interactions_by_sv(args.interactions, sv_df, args.min_distance)
    
    # Save results
    filtered.to_csv(args.output, index=False)
    print(f"Saved filtered interactions to {args.output}")

if __name__ == '__main__':
    main()