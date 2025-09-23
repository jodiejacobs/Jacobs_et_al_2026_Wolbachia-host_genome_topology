#!/usr/bin/env python3
import sys
import re

def svim_to_bed(vcf_file, bed_file):
    """Convert SVIM VCF to BED format"""
    
    with open(vcf_file, 'r') as vcf, open(bed_file, 'w') as bed:
        # Write BED header
        bed.write("#chrom\tstart\tend\tsv_type\tsv_id\tquality\tsupport\tsvlen\tstrand\n")
        
        for line in vcf:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1]) - 1  # Convert to 0-based coordinates
            sv_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            info = fields[7]
            
            # Parse INFO field for SVIM-specific tags
            sv_type = "UNKNOWN"
            end_pos = pos + len(ref)  # Default end position
            support = "0"
            svlen = "0"
            strand = "."
            
            for item in info.split(';'):
                if item.startswith('SVTYPE='):
                    sv_type = item.split('=')[1]
                elif item.startswith('END='):
                    end_pos = int(item.split('=')[1])
                elif item.startswith('SUPPORT='):
                    support = item.split('=')[1]
                elif item.startswith('SVLEN='):
                    svlen = item.split('=')[1]
                elif item.startswith('STRAND='):
                    strand = item.split('=')[1]
            
            # Handle different SV types
            if sv_type == "INS":
                end_pos = pos + 1  # Insertions are point events
            elif sv_type == "BND":
                # For breakends, use a small window around the breakpoint
                end_pos = pos + 10
            
            # Write BED entry
            bed.write(f"{chrom}\t{pos}\t{end_pos}\t{sv_type}\t{sv_id}\t{qual}\t{support}\t{svlen}\t{strand}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 svim_to_bed.py input.vcf output.bed")
        sys.exit(1)
    
    svim_to_bed(sys.argv[1], sys.argv[2])
    print(f"Converted {sys.argv[1]} to {sys.argv[2]}")
