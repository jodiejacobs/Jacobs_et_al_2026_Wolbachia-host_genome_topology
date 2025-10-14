#!/bin/bash

mkdir -p resources/chipseq

# Genome and annotations
echo "Downloading genome..."
wget -O resources/dmel-all-chromosome-r6.55.fasta.gz \
    https://ftp.flybase.net/releases/current/dmel_r6.55/fasta/dmel-all-chromosome-r6.55.fasta.gz
gunzip resources/dmel-all-chromosome-r6.55.fasta.gz
samtools faidx resources/dmel-all-chromosome-r6.55.fasta

echo "Downloading gene annotations..."
wget -O resources/genes.gtf.gz \
    https://ftp.flybase.net/releases/current/dmel_r6.55/gtf/dmel-all-r6.55.gtf.gz

# REDfly
echo "Downloading TFBS..."
wget -O resources/redfly_tfbs.gff \
    http://redfly.ccr.buffalo.edu/datadumps/redfly_tfbs.gff
wget -O resources/redfly_crm.gff \
    http://redfly.ccr.buffalo.edu/datadumps/redfly_crm.gff

# ChIP-seq data - ADD ACTUAL URLS HERE
echo "Downloading ChIP-seq data..."
# Replace ENCFF### with actual accessions you find

# H3K36me3 - CRITICAL FOR NDF
wget -O resources/chipseq/H3K36me3.bed.gz \
    "https://www.encodeproject.org/files/ENCFF###/@@download/ENCFF###.bed.gz"

# H4K16ac - CRITICAL FOR MSL
wget -O resources/chipseq/H4K16ac.bed.gz \
    "https://www.encodeproject.org/files/ENCFF###/@@download/ENCFF###.bed.gz"

# MSL2 - CRITICAL FOR DOSAGE COMP
wget -O resources/chipseq/MSL2.bed.gz \
    "https://www.encodeproject.org/files/ENCFF###/@@download/ENCFF###.bed.gz"

# Add other marks...

echo "Done!"
