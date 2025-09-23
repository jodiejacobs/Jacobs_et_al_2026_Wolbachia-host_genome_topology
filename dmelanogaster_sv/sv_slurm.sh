#!/bin/bash
#SBATCH --job-name=svim_detection
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=svim_%j.out
#SBATCH --error=svim_%j.err
#SBATCH --partition=medium

#source ~/.bashrc
#mamba activate sv_masking

# Set variables
REF="dmel-all-chromosome-r6.63.fasta.gz"
READS="Dmel_all_reads.rmdup.fastq.gz"
THREADS=${SLURM_CPUS_PER_TASK}
OUTPUT_DIR="svim_results_$(date +%Y%m%d)"

mkdir -p ${OUTPUT_DIR}

echo "Starting SVIM analysis at $(date)"
echo "Reference: ${REF}"
echo "Reads: ${READS}"
echo "Threads: ${THREADS}"

# Step 1: Map reads
echo "Mapping reads..."
minimap2 -ax map-ont ${REF} ${READS} | \
samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/aligned_reads.bam

samtools index ${OUTPUT_DIR}/aligned_reads.bam

# Step 2: Call SVs with SVIM
echo "Calling structural variants..."
svim alignment ./svim_output \
     aligned_reads.bam \
     dmel-all-chromosome-r6.63.fasta.gz \
     --minimum_length 50 \
     --minimum_depth 3

# Step 3: Filter and convert
echo "Filtering and converting to BED..."
bcftools filter -i 'QUAL>=10 && INFO/SUPPORT>=3' \
    ${OUTPUT_DIR}/svim_output/variants.vcf > ${OUTPUT_DIR}/filtered_svs.vcf

python3 svim_to_bed.py ${OUTPUT_DIR}/filtered_svs.vcf ${OUTPUT_DIR}/sv_regions.bed

# Step 4: Generate summary
echo "Generating summary..."
echo "SVIM Analysis Complete - $(date)" > ${OUTPUT_DIR}/summary.txt
echo "Total SVs: $(grep -v '^#' ${OUTPUT_DIR}/filtered_svs.vcf | wc -l)" >> ${OUTPUT_DIR}/summary.txt
echo "BED file: ${OUTPUT_DIR}/sv_regions.bed" >> ${OUTPUT_DIR}/summary.txt

echo "Analysis complete! Results in: ${OUTPUT_DIR}"
