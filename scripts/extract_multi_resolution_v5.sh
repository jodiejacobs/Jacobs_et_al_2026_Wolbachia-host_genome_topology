#!/bin/bash
# Multi-Replicate mcool Extraction Script for Snakemake
# This script extracts both cis and trans interactions from multiple replicates
# at different resolutions and creates a consolidated output file
#
# Usage: ./extract_multi_resolution_v5.sh <input_mcool_rep1> <input_mcool_rep2> <output_file>

# Handle command line arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_mcool_rep1> <input_mcool_rep2> <output_file>"
    exit 1
fi

INPUT_MCOOL_REP1=$1
INPUT_MCOOL_REP2=$2
OUTPUT_FILE=$3

# Extract sample name from output file
SAMPLE_NAME=$(basename "$OUTPUT_FILE" | cut -d'_' -f1)
echo "Processing sample: $SAMPLE_NAME"

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Create temporary directory for intermediate files
TEMP_DIR="${OUTPUT_DIR}/temp_${SAMPLE_NAME}"
mkdir -p "$TEMP_DIR"

# Define resolutions to process
resolutions=(1000 8000 32000 128000)
echo "Will process resolutions: ${resolutions[*]}"

# Define mcool files with replicates
SAMPLE_REPS=(
  $INPUT_MCOOL_REP1
  $INPUT_MCOOL_REP2
)

# Manually define chromosomes for Drosophila melanogaster
CHROMS=("2L" "2R" "3L" "3R" "4" "X" "Y")

echo "Will process chromosomes: ${CHROMS[*]}"
echo "Using resolutions: ${resolutions[*]}"

# Check if files exist
for file in "${SAMPLE_REPS[@]}"; do
  if [ ! -f "$file" ]; then
    echo "Warning: File $file not found!"
    echo "Make sure the file paths are correct."
    exit 1
  fi
done

# Create directory for each resolution
for res in "${resolutions[@]}"; do
  mkdir -p "${TEMP_DIR}/res_${res}"
done

#
# Function to process cis interactions for a replicate
#
process_cis_interactions() {
  local condition=$1
  local replicate_num=$2
  local mcool_file=$3
  local resolution=$4
  local output_file="${TEMP_DIR}/res_${resolution}/${condition}-rep${replicate_num}_cis_contacts_${resolution}bp.tsv"
  
  echo "Processing ${condition} replicate ${replicate_num} cis interactions at ${resolution}bp"
  
  # Create header
  echo -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcount" > "$output_file"
  
  # Process each chromosome for cis interactions
  for chr in "${CHROMS[@]}"; do
    echo "  Extracting cis interactions for $chr"
    
    # Extract interactions and filter out NaN values
    cooler dump -t pixels --join "${mcool_file}::resolutions/${resolution}" -r "$chr" -r2 "$chr" | \
      awk 'NR>1 && !($8 == "NaN" || $8 == "" || $8 == 0)' >> "$output_file"
    
    # Check if any data was extracted
    if [ $(wc -l < "$output_file") -le 1 ]; then
      echo "  Warning: No valid balanced interactions found for $chr. Trying without filtering..."
      
      # Try without filtering NaN values
      cooler dump -t pixels --join "${mcool_file}::resolutions/${resolution}" -r "$chr" -r2 "$chr" | \
        awk 'NR>1' >> "$output_file"
    fi
  done
  
  # Check if we got any data
  if [ $(wc -l < "$output_file") -le 1 ]; then
    echo "  Warning: No cis interactions extracted for ${condition} replicate ${replicate_num} at ${resolution}bp"
  else
    echo "  Successfully extracted $(( $(wc -l < "$output_file") - 1 )) cis interactions"
  fi
}

#
# Function to process trans interactions for a replicate
#
process_trans_interactions() {
  local condition=$1
  local replicate_num=$2
  local mcool_file=$3
  local resolution=$4
  local output_file="${TEMP_DIR}/res_${resolution}/${condition}-rep${replicate_num}_trans_contacts_${resolution}bp.tsv"
  
  echo "Processing ${condition} replicate ${replicate_num} trans interactions at ${resolution}bp"
  
  # Create header
  echo -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcount" > "$output_file"
  
  # Process chromosome pairs for trans interactions
  for ((i=0; i<${#CHROMS[@]}; i++)); do
    for ((j=i+1; j<${#CHROMS[@]}; j++)); do
      chrA="${CHROMS[i]}"
      chrB="${CHROMS[j]}"
      echo "  Extracting trans interactions for $chrA vs $chrB"
      
      # Extract trans interactions and filter out NaN values
      cooler dump -t pixels --join "${mcool_file}::resolutions/${resolution}" -r "$chrA" -r2 "$chrB" | \
        awk 'NR>1 && !($8 == "NaN" || $8 == "" || $8 == 0)' >> "$output_file"
    done
  done
  
  # Check if we got any data
  if [ $(wc -l < "$output_file") -le 1 ]; then
    echo "  Warning: No trans interactions extracted for ${condition} replicate ${replicate_num} at ${resolution}bp"
    
    # Try again without filtering
    echo "  Trying without filtering NaN values..."
    rm "$output_file"
    echo -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcount" > "$output_file"
    
    for ((i=0; i<${#CHROMS[@]}; i++)); do
      for ((j=i+1; j<${#CHROMS[@]}; j++)); do
        chrA="${CHROMS[i]}"
        chrB="${CHROMS[j]}"
        
        cooler dump -t pixels --join "${mcool_file}::resolutions/${resolution}" -r "$chrA" -r2 "$chrB" | \
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
# Main processing loop
#

# Process each resolution
for res in "${resolutions[@]}"; do
  echo "==== Processing resolution: ${res}bp ===="
  
  # Process sample replicates
  echo "Processing ${SAMPLE_NAME} replicates at ${res}bp resolution..."
  for i in "${!SAMPLE_REPS[@]}"; do
    rep_num=$((i+1))
    process_cis_interactions "${SAMPLE_NAME}" "$rep_num" "${SAMPLE_REPS[$i]}" "$res"
    process_trans_interactions "${SAMPLE_NAME}" "$rep_num" "${SAMPLE_REPS[$i]}" "$res"
  done

  echo "Completed processing at ${res}bp resolution"
done

#
# Create consolidated output file
#
echo "Creating consolidated output file..."

# Create the final output file with headers
echo -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcount\treplicate\tcondition\tinteraction_type\tresolution" > "$OUTPUT_FILE"

# Add all data to the consolidated file
for res in "${resolutions[@]}"; do
  echo "Adding resolution ${res}bp data to consolidated file..."
  
  # Add cis interactions
  for i in "${!SAMPLE_REPS[@]}"; do
    rep_num=$((i+1))
    cis_file="${TEMP_DIR}/res_${res}/${SAMPLE_NAME}-rep${rep_num}_cis_contacts_${res}bp.tsv"
    if [ -f "$cis_file" ] && [ $(wc -l < "$cis_file") -gt 1 ]; then
      echo "  Adding cis interactions from replicate ${rep_num}..."
      awk -v rep="$rep_num" -v res="$res" -v cond="$SAMPLE_NAME" 'NR>1 {print $0 "\t" rep "\t" cond "\tcis\t" res}' "$cis_file" >> "$OUTPUT_FILE"
    fi
  done
  
  # Add trans interactions
  for i in "${!SAMPLE_REPS[@]}"; do
    rep_num=$((i+1))
    trans_file="${TEMP_DIR}/res_${res}/${SAMPLE_NAME}-rep${rep_num}_trans_contacts_${res}bp.tsv"
    if [ -f "$trans_file" ] && [ $(wc -l < "$trans_file") -gt 1 ]; then
      echo "  Adding trans interactions from replicate ${rep_num}..."
      awk -v rep="$rep_num" -v res="$res" -v cond="$SAMPLE_NAME" 'NR>1 {print $0 "\t" rep "\t" cond "\ttrans\t" res}' "$trans_file" >> "$OUTPUT_FILE"
    fi
  done
done

# Count the number of rows in the final output file
TOTAL_ROWS=$(wc -l < "$OUTPUT_FILE")
TOTAL_INTERACTIONS=$((TOTAL_ROWS - 1))

echo "All interactions extracted and consolidated successfully!"
echo "Total interactions: $TOTAL_INTERACTIONS"
echo "Results saved to $OUTPUT_FILE"

# Keep temporary files for debugging
echo "Keeping temporary files for debugging in: $TEMP_DIR"
# rm -rf "$TEMP_DIR"

echo "Processing complete."