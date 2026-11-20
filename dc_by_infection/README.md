# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Hi-C/Micro-C analysis pipeline studying differential chromatin interactions between Wolbachia-infected (wMel) and uninfected (DOX-treated) Drosophila melanogaster samples. The analysis focuses on:
- Differential chromatin interactions using diffHic
- X chromosome dosage compensation changes
- Enhancer-promoter interaction changes
- TAD and loop boundary alterations
- Architectural protein (CTCF, BEAF-32, etc.) enrichment at differential contact sites

## Running the Pipeline

### Main Command

Run the complete Snakemake workflow using SLURM:
```bash
mamba activate snakemake
snakemake --executor slurm --default-resources slurm_partition=medium runtime=720 mem_mb=1000000 -j 10 -s Snakefile
```

Add `-n` for a dry run to preview what will execute.

### Individual Rule Execution

Run specific analysis steps:
```bash
# Run only diffHic analysis
snakemake --executor slurm -j 10 diffhic_analysis -n

# Run only chromatin structure analysis
snakemake --executor slurm -j 10 chromatin_structure_analysis -n

# Run specific sample contact extraction
snakemake --executor slurm -j 10 results/contacts/JW18wMel_contacts.tsv -n
```

### Conda Environments

The pipeline uses three main conda environments (paths are hardcoded in Snakefile:31-34):
- `micro-c`: For cooler/cooltools operations (contact extraction, TAD calling)
- `diffhic`: For R-based differential interaction analysis
- `micro-c`: For Python-based analyses (filtering, enrichment, enhancer analysis)

## Architecture

### Data Flow

1. **Contact Extraction** (`identify_contacts`): Extracts multi-resolution contacts from `.mcool` files for both replicates
   - Input: Two `.mcool` replicates per sample
   - Output: `results/contacts/{sample}_contacts.tsv`
   - Script: `scripts/extract_multi_resolution_v5.sh`
   - Resolutions: 128000, 32000, 8000, 1000 bp

2. **SV Filtering** (`filter_sv_interactions`): Removes interactions overlapping structural variants
   - Input: Contact TSV + VCF file of structural variants
   - Output: `results/sv_filtered/{sample}_filtered_contacts.csv`
   - Script: `scripts/filter_sv_interactions.py`
   - Method: Adds 200bp padding around SV breakpoints (following Cubenas-Potts 2017)

3. **DiffHic Analysis** (`diffhic_analysis`): Identifies differential interactions using edgeR
   - Input: Filtered contacts from all samples
   - Output: Per-resolution results, combined tables, volcano plots, null model
   - Script: `scripts/analyze_multi_resolution_4.R`
   - Creates comprehensive summary and null model for downstream analyses

4. **Downstream Analyses**: All use the combined diffHic results and null model
   - Chromatin structure (TADs, loops, compartments) via cooltools
   - Architectural protein enrichment at differential contacts
   - Enhancer-enhancer and enhancer-promoter interaction changes
   - X chromosome-specific dosage compensation analysis

### Key Design Patterns

**Null Model Usage**: The diffHic analysis generates a null model (`null_model_results.csv`) representing expected variation between technical replicates. All downstream scripts use this to:
- Compare observed differential interactions against technical noise
- Calculate z-scores and enrichment significance
- Filter out interactions within expected technical variation

**Multi-Resolution Analysis**: Contacts are extracted at 4 resolutions (1kb-128kb) to capture both:
- Fine-scale enhancer-promoter interactions (1-8kb)
- Larger chromatin domain changes (32-128kb)

**Chromosome Handling**: Drosophila chromosomes are referenced without "chr" prefix (2L, 2R, 3L, 3R, 4, X, Y). Scripts handle both formats and standardize internally.

## Reference Files

Located in `/private/groups/russelllab/jodie/wolbachia_induced_DE/Jacobs_et_al_2026_Wolbachia-host_genome_topology/reference_files/`:
- `dm6.genome`: Chromosome sizes
- `JW18_structural_variants.vcf`: Strain-specific SVs to filter
- `chip_peaks/`: Directory of ChIP-seq BED files for architectural proteins
- `dm6_enhancers.bed`: Enhancer annotations
- `enhancer_classification.tsv`: Enhancer functional classes
- `dm6_tss.bed`: Transcription start sites
- `dm6_has_ces_sites.bed`: High-affinity sites (HAS) for dosage compensation complex

## Important Implementation Details

### Script Versioning
The `scripts/` directory contains versioned files (v0, v1, etc.). The Snakefile always uses the latest non-versioned filename (e.g., `filter_sv_interactions.py`), which represents the current working version.

### Memory and Threading
- DiffHic analysis: 128GB RAM, 125 threads (Snakefile:129-131)
- Chromatin structure: 200GB RAM, 16 threads
- Most Python scripts: 100-128GB RAM, 4-8 threads

### Cooler Resolution Handling
Scripts that load `.mcool` files use `get_closest_resolution()` to handle cases where exact resolutions may not be available, selecting the nearest available resolution.

### ChIP-seq Peak Loading
`strain_specific_architectural_enrichment.py` only loads validated ChIP-seq files (WORKING_FILES dict at line 61) to avoid parsing errors from header-only or malformed files.

## Common Development Workflows

### Adding a New Analysis

1. Add new rule to Snakefile with appropriate dependencies
2. Use the null model for statistical comparison: `input.null_model = f"{output_dir}/diffhic_results/res_8000/null_model_results.csv"`
3. Use combined interactions: `input.interactions = f"{output_dir}/diffhic_results/comprehensive_tables/all_sig_interactions_combined.csv"`
4. Activate appropriate conda environment in shell command
5. Add output to `rule all` inputs to include in default workflow

### Modifying DiffHic Parameters

Key parameters in Snakefile (lines 8-14):
- `samples`: Sample names matching input .mcool prefixes
- `conditions`: Condition labels for comparison
- `resolutions`: List of resolutions to analyze (reversed order for processing)

DiffHic script parameters (passed via CLI, defined in Snakefile:137-142):
- `--fdr`: False discovery rate threshold (default 0.05)
- `--resolution`: Comma-separated list of resolutions
- `--threads`: Parallelization level

### Testing Individual Scripts

All Python scripts accept `--help` to show parameters:
```bash
conda activate micro-c
python scripts/filter_sv_interactions.py --help
python scripts/enhancer_class_analysis.py --help
```

R scripts use optparse:
```bash
conda activate diffhic
Rscript scripts/analyze_multi_resolution_4.R --help
```
