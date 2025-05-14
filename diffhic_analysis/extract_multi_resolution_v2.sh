#!/bin/bash
# Multi-Resolution mcool Extraction Script works with snakemake pipeline
# This script extracts both cis and trans interactions from two replicates
# at specified resolution for differential analysis
#
# Usage: ./extract_multi_resolution.sh [cooler_file1] [cooler_file2] [output_file]

# Handle command line arguments
COOLER1="$1"
COOLER2="$2"
OUTPUT_FILE="$3"

# Validate inputs
if [ -z "$COOLER1" ] || [ -z "$COOLER2" ] || [ -z "$OUTPUT_FILE" ]; then
    echo "Usage: $0 <cooler_file1> <cooler_file2> <output_file>"
    exit 1
fi

# Check if files exist
for file in "$COOLER1" "$COOLER2"; do
  if [ ! -f "$file" ]; then
    echo "Error: File $file not found!"
    exit 1
  fi
done

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Extract sample name from cooler file path
SAMPLE_NAME=$(basename "$COOLER1" | cut -d'-' -f1-2)
echo "Processing sample: $SAMPLE_NAME"

# Define resolution - default to 1000
RESOLUTION=1000

# Manually define chromosomes for Drosophila melanogaster
CHROMS=("2L" "2R" "3L" "3R" "4" "X" "Y")
echo "Will process chromosomes: ${CHROMS[*]}"

# Create temp directories for intermediate files
TEMP_DIR="${OUTPUT_DIR}/temp_${SAMPLE_NAME}"
mkdir -p "$TEMP_DIR"

#
# Function to process cis interactions for a replicate
#
process_cis_interactions() {
  local mcool_file=$1
  local replicate_num=$2
  local output_file="${TEMP_DIR}/${SAMPLE_NAME}_rep${replicate_num}_cis_contacts_${RESOLUTION}bp.tsv"
  
  echo "Processing ${SAMPLE_NAME} replicate ${replicate_num} cis interactions at ${RESOLUTION}bp"
  
  # Create header
  echo -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcount" > "$output_file"
  
  # Process each chromosome for cis interactions
  for chr in "${CHROMS[@]}"; do
    echo "  Extracting cis interactions for $chr"
    
    # Extract interactions and filter out NaN values
    cooler dump -t pixels --join "${mcool_file}::resolutions/${RESOLUTION}" -r "$chr" -r2 "$chr" | \
      awk 'NR>1 && !($8 == "NaN" || $8 == "" || $8 == 0)' >> "$output_file"
    
    # Check if any data was extracted
    if [ $(wc -l < "$output_file") -le 1 ]; then
      echo "  Warning: No valid balanced interactions found for $chr. Trying without filtering..."
      
      # Try without filtering NaN values
      cooler dump -t pixels --join "${mcool_file}::resolutions/${RESOLUTION}" -r "$chr" -r2 "$chr" | \
        awk 'NR>1' >> "$output_file"
    fi
  done
  
  # Check if we got any data
  if [ $(wc -l < "$output_file") -le 1 ]; then
    echo "  Warning: No cis interactions extracted for ${SAMPLE_NAME} replicate ${replicate_num} at ${RESOLUTION}bp"
  else
    echo "  Successfully extracted $(( $(wc -l < "$output_file") - 1 )) cis interactions"
  fi
}

#
# Function to process trans interactions for a replicate
#
process_trans_interactions() {
  local mcool_file=$1
  local replicate_num=$2
  local output_file="${TEMP_DIR}/${SAMPLE_NAME}_rep${replicate_num}_trans_contacts_${RESOLUTION}bp.tsv"
  
  echo "Processing ${SAMPLE_NAME} replicate ${replicate_num} trans interactions at ${RESOLUTION}bp"
  
  # Create header
  echo -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcount" > "$output_file"
  
  # Process chromosome pairs for trans interactions
  for ((i=0; i<${#CHROMS[@]}; i++)); do
    for ((j=i+1; j<${#CHROMS[@]}; j++)); do
      chrA="${CHROMS[i]}"
      chrB="${CHROMS[j]}"
      echo "  Extracting trans interactions for $chrA vs $chrB"
      
      # Extract trans interactions and filter out NaN values
      cooler dump -t pixels --join "${mcool_file}::resolutions/${RESOLUTION}" -r "$chrA" -r2 "$chrB" | \
        awk 'NR>1 && !($8 == "NaN" || $8 == "" || $8 == 0)' >> "$output_file"
    done
  done
  
  # Check if we got any data
  if [ $(wc -l < "$output_file") -le 1 ]; then
    echo "  Warning: No trans interactions extracted for ${SAMPLE_NAME} replicate ${replicate_num} at ${RESOLUTION}bp"
    
    # Try again without filtering
    echo "  Trying without filtering NaN values..."
    rm "$output_file"
    echo -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcount" > "$output_file"
    
    for ((i=0; i<${#CHROMS[@]}; i++)); do
      for ((j=i+1; j<${#CHROMS[@]}; j++)); do
        chrA="${CHROMS[i]}"
        chrB="${CHROMS[j]}"
        
        cooler dump -t pixels --join "${mcool_file}::resolutions/${RESOLUTION}" -r "$chrA" -r2 "$chrB" | \
          awk 'NR>1' >> "$output_file"
      done
    done
    
    if [ $(wc -l < "$output_file") -le 1 ]; then
      echo "  Warning: Still no trans interactions found even without filtering."
    else
      echo "  Successfully extracted $(( $(wc -l < "$output_file") - 1 )) unfiltered trans interactions"
    fi
  else
    echo "  Successfully extracted $(( $(wc -l < "$output_file") - 1 )) trans interactions"
  fi
}

#
# Main processing
#

# Process replicate 1
process_cis_interactions "$COOLER1" "1"
process_trans_interactions "$COOLER1" "1"

# Process replicate 2
process_cis_interactions "$COOLER2" "2"
process_trans_interactions "$COOLER2" "2"

#
# Create merged files for final output
#
echo "Creating merged interaction files..."

# Merge cis interactions from both replicates
echo -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcount\treplicate\tcondition" > "$OUTPUT_FILE"

# Add replicate 1 cis interactions
cis_file1="${TEMP_DIR}/${SAMPLE_NAME}_rep1_cis_contacts_${RESOLUTION}bp.tsv"
if [ -f "$cis_file1" ] && [ $(wc -l < "$cis_file1") -gt 1 ]; then
  awk -v sample="$SAMPLE_NAME" 'NR>1 {print $0 "\t1\t" sample}' "$cis_file1" >> "$OUTPUT_FILE"
fi

# Add replicate 2 cis interactions
cis_file2="${TEMP_DIR}/${SAMPLE_NAME}_rep2_cis_contacts_${RESOLUTION}bp.tsv"
if [ -f "$cis_file2" ] && [ $(wc -l < "$cis_file2") -gt 1 ]; then
  awk -v sample="$SAMPLE_NAME" 'NR>1 {print $0 "\t2\t" sample}' "$cis_file2" >> "$OUTPUT_FILE"
fi

# Add replicate 1 trans interactions
trans_file1="${TEMP_DIR}/${SAMPLE_NAME}_rep1_trans_contacts_${RESOLUTION}bp.tsv"
if [ -f "$trans_file1" ] && [ $(wc -l < "$trans_file1") -gt 1 ]; then
  awk -v sample="$SAMPLE_NAME" 'NR>1 {print $0 "\t1\t" sample}' "$trans_file1" >> "$OUTPUT_FILE"
fi

# Add replicate 2 trans interactions
trans_file2="${TEMP_DIR}/${SAMPLE_NAME}_rep2_trans_contacts_${RESOLUTION}bp.tsv"
if [ -f "$trans_file2" ] && [ $(wc -l < "$trans_file2") -gt 1 ]; then
  awk -v sample="$SAMPLE_NAME" 'NR>1 {print $0 "\t2\t" sample}' "$trans_file2" >> "$OUTPUT_FILE"
fi

# Clean up temporary files (comment out if you want to keep them)
rm -rf "$TEMP_DIR"

# Report completion
echo "All interactions extracted and merged successfully!"
echo "Results saved to $OUTPUT_FILE"