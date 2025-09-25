#!/usr/bin/env Rscript
# Single Resolution diffHic Analysis for Differential Chromatin Interactions
# 
# This script performs differential analysis of chromatin interactions between
# multiple Wolbachia infections at a single resolution using the diffHic package.
#
# Usage: Rscript analyze_single_resolution.R --data_dir interactions_data --resolution 1000 
#                                          --min_count 10 --fdr 0.01 --output_dir diffhic_results

#=============================================================================
# 1. Setup Environment and Load Libraries
#=============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicInteractions)
  library(InteractionSet)
  library(diffHic)
  library(edgeR)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(pheatmap)
  library(VennDiagram)
  library(gridExtra)
  library(data.table)
  library(BiocParallel)
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("--data_dir"), type="character", default="interactions_data",
              help="Directory containing the contact data [default=%default]"),
  make_option(c("--resolution"), type="integer", default=1000,
              help="Resolution to analyze [default=%default]"),
  make_option(c("--min_count"), type="integer", default=10,
              help="Minimum read count to keep an interaction [default=%default]"),
  make_option(c("--min_samples"), type="integer", default=2,
              help="Minimum number of samples with min_count [default=%default]"),
  make_option(c("--fdr"), type="double", default=0.01,
              help="FDR threshold for significance [default=%default]"),
  make_option(c("--output_dir"), type="character", default="diffhic_results",
              help="Output directory [default=%default]"),
  make_option(c("--reference"), type="character", default="DOX",
              help="Reference level for comparisons [default=%default]"),
  make_option(c("--threads"), type="integer", default=4,
              help="Number of CPU threads to use [default=%default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list, 
                          description="Single-resolution differential Micro-C analysis")
opt <- parse_args(opt_parser)

# Set parameters from command line arguments
data_dir <- opt$data_dir
resolution <- opt$resolution
min_count <- opt$min_count
min_samples <- opt$min_samples
fdr_threshold <- opt$fdr
output_dir <- opt$output_dir
reference_level <- opt$reference
threads <- opt$threads

# Register parallel processing
register(MulticoreParam(threads))

cat("Using parameters:\n")
cat("  Data directory:", data_dir, "\n")
cat("  Resolution:", resolution, "bp\n")
cat("  Min count:", min_count, "\n")
cat("  Min samples:", min_samples, "\n")
cat("  FDR threshold:", fdr_threshold, "\n")
cat("  Output directory:", output_dir, "\n")
cat("  Reference level:", reference_level, "\n")
cat("  Threads:", threads, "\n")

# Create output directories
dir.create(output_dir, showWarnings = FALSE)
output_res_dir <- file.path(output_dir, paste0("res_", resolution))
dir.create(output_res_dir, showWarnings = FALSE)

#=============================================================================
# 2. Define Functions for Processing
#=============================================================================

# Function to create a unique identifier for interactions
create_interaction_id <- function(chr1, start1, end1, chr2, start2, end2) {
  # Sort chromosome pairs to ensure consistency
  is_swap <- as.character(chr2) < as.character(chr1) | 
            (as.character(chr2) == as.character(chr1) & start2 < start1)
  
  # Swap if needed to ensure consistent ordering
  result <- character(length(chr1))
  
  for (i in 1:length(chr1)) {
    if (is_swap[i]) {
      # Swap chromosome information for consistent IDs
      id <- paste(
        chr2[i], start2[i], end2[i],
        chr1[i], start1[i], end1[i],
        sep = "_"
      )
    } else {
      id <- paste(
        chr1[i], start1[i], end1[i],
        chr2[i], start2[i], end2[i],
        sep = "_"
      )
    }
    result[i] <- id
  }
  
  return(result)
}

# Function to read data from the Snakemake pipeline output with flexible sample handling
read_all_data <- function(resolution, data_dir) {
  # List all contact files in the data directory
  contact_files <- list.files(data_dir, pattern = "_contacts\\.tsv$", full.names = TRUE)
  
  cat("Found contact files:", paste(basename(contact_files), collapse=", "), "\n")
  
  # Initialize a list to store data for each sample and replicate
  sample_data <- list()
  
  # Process each file
  for (file in contact_files) {
    cat("Processing file:", basename(file), "for resolution", resolution, "bp\n")
    
    # Extract sample name from filename (adjust this pattern as needed)
    filename <- basename(file)
    sample_name <- gsub("_contacts\\.tsv$", "", filename)
    
    # Read the data
    data <- fread(file, header = TRUE)
    
    # Check if the resolution column exists
    if (!"resolution" %in% colnames(data)) {
      # If resolution not in columns, check if the file is already filtered by resolution
      # or try to infer resolution from filename patterns
      if (grepl(paste0("_", resolution, "bp"), filename)) {
        # File is already filtered for this resolution
        resolution_data <- data
        resolution_data <- resolution_data[resolution_data$count >= min_count, ]
      } else {
        warning("Resolution column not found in ", filename, 
                " and cannot infer resolution from filename. Assuming all data is for requested resolution.")
        resolution_data <- data
        resolution_data <- resolution_data[resolution_data$count >= min_count, ]

      }
    } else {
      # Filter for the specified resolution3.
      resolution_data <- data[data$resolution == resolution, ]
      resolution_data <- resolution_data[resolution_data$count >= min_count, ]

    }
    
    if (nrow(resolution_data) == 0) {
      warning("No data found for resolution ", resolution, " in file ", filename)
      next
    }
    
    # Ensure the data has all required columns
    required_cols <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "count")
    missing_cols <- setdiff(required_cols, colnames(resolution_data))
    if (length(missing_cols) > 0) {
      stop("Missing columns in data file: ", paste(missing_cols, collapse=", "))
    }
    
    # Convert columns to numeric if they're not already
    cols_to_convert <- c("start1", "end1", "start2", "end2", "count")
    resolution_data[, (cols_to_convert) := lapply(.SD, function(x) {
      if (is.character(x)) as.numeric(gsub(",", "", x)) else x
    }), .SDcols = cols_to_convert]
    
    # Extract condition from sample name
    # Pattern for JW18-wMel_contacts.tsv, JW18-DOX_contacts.tsv, etc.
    condition <- NA
    if (grepl("JW18-wMel", sample_name)) {
      condition <- "wMel"
    } else if (grepl("JW18-wRi", sample_name)) {
      condition <- "wRi"
    } else if (grepl("JW18-wWil", sample_name)) {
      condition <- "wWil"
    } else if (grepl("JW18-DOX", sample_name)) {
      condition <- "DOX"
    } else {
      # Default to using the filename if we can't match
      condition <- sample_name
    }
    
    # Add condition column if it doesn't exist
    if (!"condition" %in% colnames(resolution_data)) {
      resolution_data$condition <- condition
    }
    
    # Add replicate column if it doesn't exist
    if (!"replicate" %in% colnames(resolution_data)) {
      # Check if each row already has a replicate column we can use
      # This is relevant if data has one row per replicate
      
      # If not, check the filename for replicate info
      if (grepl("_rep1", sample_name, ignore.case = TRUE) || grepl("_r1", sample_name, ignore.case = TRUE)) {
        resolution_data$replicate <- 1
      } else if (grepl("_rep2", sample_name, ignore.case = TRUE) || grepl("_r2", sample_name, ignore.case = TRUE)) {
        resolution_data$replicate <- 2
      } else {
        # If replicate information is not in the filename,
        # check if there's a column that might contain it
        possible_rep_cols <- grep("rep|replicate|batch", colnames(resolution_data), ignore.case = TRUE, value = TRUE)
        
        if (length(possible_rep_cols) > 0) {
          cat("  - Found possible replicate column:", possible_rep_cols[1], "\n")
          resolution_data$replicate <- resolution_data[[possible_rep_cols[1]]]
        } else {
          # As a last resort, check if we can parse replicate from another pattern
          # Assuming we need to split the samples by replicate
          if (nrow(resolution_data) > 0) {
            # Try to figure out which values might be replicates
            unique_values <- unique(resolution_data$sample)
            if (length(unique_values) <= 2) {
              cat("  - Found possible replicate values in 'sample' column:", 
                 paste(unique_values, collapse=", "), "\n")
              # Map these values to replicates 1 and 2
              rep_map <- setNames(1:length(unique_values), unique_values)
              resolution_data$replicate <- rep_map[resolution_data$sample]
            } else {
              # Default to replicate 1 if we can't determine
              warning("Could not determine replicate information. Splitting by sample column.")
              sample_groups <- split(resolution_data, resolution_data$sample)
              
              # Process each sample group separately
              for (i in seq_along(sample_groups)) {
                sg <- sample_groups[[i]]
                sg$replicate <- i
                
                # Add to sample_data
                key <- paste0(condition, "_rep", i)
                sample_data[[key]] <- as.data.frame(sg)
                cat("  - Found", nrow(sg), "interactions for", key, "\n")
              }
              
              # Skip the rest of the loop since we've already added to sample_data
              next
            }
          } else {
            # Default to replicate 1
            resolution_data$replicate <- 1
            warning("Could not determine replicate. Defaulting to replicate 1.")
          }
        }
      }
    }
    
    # Split data by replicate
    for (rep in unique(resolution_data$replicate)) {
      subset_data <- as.data.frame(resolution_data[resolution_data$replicate == rep, ])
      
      if (nrow(subset_data) > 0) {
        key <- paste0(condition, "_rep", rep)
        if (is.null(sample_data[[key]])) {
          sample_data[[key]] <- subset_data
        } else {
          sample_data[[key]] <- rbind(sample_data[[key]], subset_data)
        }
        
        cat("  - Found", nrow(subset_data), "interactions for", key, "\n")
      }
    }
  }
  
  # Return all sample data
  return(sample_data)
}

# Function to create GenomicInteractions object
create_gi <- function(data) {
  # Create GRanges for anchors
  anchor1 <- GRanges(
    seqnames = data$chrom1,
    ranges = IRanges(start = data$start1, end = data$end1)
  )
  
  anchor2 <- GRanges(
    seqnames = data$chrom2,
    ranges = IRanges(start = data$start2, end = data$end2)
  )
  
  # Create GenomicInteractions object with raw counts
  gi <- GenomicInteractions(anchor1, anchor2, as.numeric(data$count))
  
  return(gi)
}

# Function to find all unique interactions across samples
find_all_interactions <- function(iset_list) {
  cat("Finding all unique interactions across samples...\n")
  
  # First approach: convert everything to data frames and work with those
  interactions_data <- list()
  
  # Extract all interactions as data frames
  for (i in seq_along(iset_list)) {
    sample_name <- names(iset_list)[i]
    cat("Processing", sample_name, "\n")
    
    current_iset <- iset_list[[i]]
    current_gi <- interactions(current_iset)
    
    # Extract regions
    regions_df <- as.data.frame(regions(current_gi))
    
    # Extract interactions
    a1_idx <- anchors(current_gi, type="first", id=TRUE)
    a2_idx <- anchors(current_gi, type="second", id=TRUE)
    
    # Get the actual GRanges for each anchor
    a1_ranges <- regions_df[a1_idx, ]
    a2_ranges <- regions_df[a2_idx, ]
    
    # Get counts
    counts <- assay(current_iset)[,1]
    
    # Combine into a single data frame
    int_df <- data.frame(
      chrom1 = a1_ranges$seqnames,
      start1 = a1_ranges$start,
      end1 = a1_ranges$end,
      chrom2 = a2_ranges$seqnames,
      start2 = a2_ranges$start,
      end2 = a2_ranges$end,
      count = counts,
      sample = i
    )
    
    interactions_data[[sample_name]] <- int_df
  }
  
  # Combine all interactions
  all_ints <- do.call(rbind, interactions_data)
  
  # Create a unique identifier for each interaction
  all_ints$int_id <- paste(
    all_ints$chrom1, all_ints$start1, all_ints$end1,
    all_ints$chrom2, all_ints$start2, all_ints$end2,
    sep = "_"
  )
  
  # Find unique interactions
  unique_ints <- unique(all_ints[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "int_id")])
  cat("Found", nrow(unique_ints), "unique interactions\n")
  
  # Create a map from int_id to row index in unique_ints
  int_id_to_index <- setNames(1:nrow(unique_ints), unique_ints$int_id)
  
  # Create count matrix
  n_unique <- nrow(unique_ints)
  n_samples <- length(iset_list)
  count_matrix <- matrix(0, nrow=n_unique, ncol=n_samples)
  colnames(count_matrix) <- names(iset_list)
  
  cat("Filling count matrix using a faster approach...\n")
  
  # Fill the count matrix using a vectorized approach
  for (j in 1:n_samples) {
    cat("Processing sample", j, "of", n_samples, "\n")
    
    # Get the subset of data for this sample
    sample_data <- all_ints[all_ints$sample == j, ]
    
    # For each interaction in this sample, add its count to the right position in the matrix
    if (nrow(sample_data) > 0) {
      # Look up the row indices in the count matrix
      row_indices <- int_id_to_index[sample_data$int_id]
      
      # Process in chunks to avoid memory issues
      chunk_size <- 100000
      n_chunks <- ceiling(length(row_indices) / chunk_size)
      
      for (chunk in 1:n_chunks) {
        cat("  Processing chunk", chunk, "of", n_chunks, "\n")
        start_idx <- (chunk - 1) * chunk_size + 1
        end_idx <- min(chunk * chunk_size, length(row_indices))
        
        # Get the subset for this chunk
        chunk_indices <- row_indices[start_idx:end_idx]
        chunk_counts <- sample_data$count[start_idx:end_idx]
        
        # Aggregate by row index (some interactions might appear multiple times)
        for (k in seq_along(chunk_indices)) {
          count_matrix[chunk_indices[k], j] <- count_matrix[chunk_indices[k], j] + chunk_counts[k]
        }
      }
    }
  }
  
  # Create GInteractions object
  cat("Creating GInteractions object...\n")
  anchor1 <- GRanges(
    seqnames = unique_ints$chrom1,
    ranges = IRanges(start = unique_ints$start1, end = unique_ints$end1)
  )
  
  anchor2 <- GRanges(
    seqnames = unique_ints$chrom2,
    ranges = IRanges(start = unique_ints$start2, end = unique_ints$end2)
  )
  
  all_gi <- GInteractions(anchor1, anchor2)
  
  # Create combined InteractionSet
  cat("Creating final InteractionSet...\n")
  combined_iset <- InteractionSet(list(counts=count_matrix), interactions=all_gi)
  
  # Final verification
  cat("Count matrix summary:", summary(as.vector(count_matrix)), "\n")
  cat("Non-zero entries:", sum(count_matrix > 0), "\n")
  
  return(combined_iset)
}

# Optimized function to add annotations to results
add_annotations <- function(results_table, iset) {
  # Get the interactions
  gi <- interactions(iset)
  
  # Get anchors
  anchor1 <- anchors(gi, type="first")
  anchor2 <- anchors(gi, type="second")
  
  # Determine if each interaction is cis or trans - vectorized
  chr1 <- as.character(seqnames(anchor1))
  chr2 <- as.character(seqnames(anchor2))
  is_cis <- chr1 == chr2
  
  # Calculate interaction distance for cis interactions - vectorized
  interaction_distance <- rep(NA, length(gi))
  if (any(is_cis)) {
    # For cis interactions, calculate distance as end2 - start1
    # This is more efficient than a loop
    interaction_distance[is_cis] <- end(anchor2)[is_cis] - start(anchor1)[is_cis]
  }
  
  # Create a data frame with interaction info - all at once
  annotations <- data.frame(
    chr1 = chr1,
    start1 = start(anchor1),
    end1 = end(anchor1),
    chr2 = chr2,
    start2 = start(anchor2),
    end2 = end(anchor2),
    interaction_type = ifelse(is_cis, "cis", "trans"),
    interaction_distance = interaction_distance
  )
  
  # Make sure row names match
  row.names(annotations) <- row.names(results_table)
  
  # Combine with results
  combined_results <- cbind(results_table, annotations)
  return(combined_results)
}

# Function to convert GenomicInteractions to GInteractions
convertToGInteractions <- function(gi) {
  # Extract the anchors from the GenomicInteractions object
  anchor1 <- anchors(gi, type="first")
  anchor2 <- anchors(gi, type="second")
  
  # Create GInteractions object
  gint <- GInteractions(anchor1, anchor2)
  
  # Copy the metadata columns if needed
  mcols(gint) <- mcols(gi)
  
  return(gint)
}

#=============================================================================
# 3. Process Single Resolution
#=============================================================================

cat("\n=======================================================\n")
cat("Processing resolution:", resolution, "bp\n")
cat("=======================================================\n")

# Load data from the combined files
all_data_list <- read_all_data(resolution, data_dir)

# Check if we have any data
if (length(all_data_list) == 0) {
  stop("No data found for resolution ", resolution, "bp")
}

# Convert columns to numeric if not already
cols_to_numeric <- c("start1", "end1", "start2", "end2", "count")

# Process each dataset separately
for (dataset_name in names(all_data_list)) {
  # Convert each data frame to data.table first if needed
  if (!is.data.table(all_data_list[[dataset_name]])) {
    all_data_list[[dataset_name]] <- as.data.table(all_data_list[[dataset_name]])
  }
  
  # Now convert the numeric columns
  all_data_list[[dataset_name]][, (cols_to_numeric) := lapply(.SD, as.numeric), 
                                .SDcols = cols_to_numeric]
}

# Create GInteractions objects for each sample
gi_list <- list()
iset_list <- list()

for (sample_name in names(all_data_list)) {
  data <- all_data_list[[sample_name]]
  if (nrow(data) > 0) {
    # Create GRanges for anchors
    anchor1 <- GRanges(
      seqnames = data$chrom1,
      ranges = IRanges(start = data$start1, end = data$end1)
    )
    
    anchor2 <- GRanges(
      seqnames = data$chrom2,
      ranges = IRanges(start = data$start2, end = data$end2)
    )
    
    # Create GInteractions object
    gint <- GInteractions(anchor1, anchor2)
    gi_list[[sample_name]] <- gint
    
    # Create InteractionSet object
    iset_list[[sample_name]] <- InteractionSet(list(counts=matrix(data$count, ncol=1)), gint)
  }
}

# Create unified InteractionSet with all samples
cat("Creating combined interaction set...\n")
combined_iset <- find_all_interactions(iset_list)

# Calculate library sizes
lib_sizes <- colSums(assay(combined_iset))
combined_iset$totals <- lib_sizes

# Check the dimensions of the combined object
cat("Combined interaction set has", nrow(combined_iset), "interactions across", 
    ncol(assay(combined_iset)), "samples\n")

# Filter out low-abundance interactions
cat("Filtering low-abundance interactions...\n")
keep <- rowSums(assay(combined_iset) >= min_count) >= min_samples
filtered_iset <- combined_iset[keep,]

# Check filtering results
cat("Filtered from", nrow(combined_iset), "to", nrow(filtered_iset), "interactions\n")

# Create DGEList object
cat("Creating DGEList object for analysis...\n")
y <- asDGEList(filtered_iset)

# Extract sample information from column names
samples <- colnames(y)

# Parse condition and replicate from sample names (assumes format like condition_rep#)
conditions <- sapply(strsplit(samples, "_rep"), function(x) x[1])
replicates <- sapply(strsplit(samples, "_rep"), function(x) x[2])

# Create condition factor - ensure reference level is first
all_conditions <- unique(conditions)
if (reference_level %in% all_conditions) {
  # Put reference level first
  condition_levels <- c(reference_level, setdiff(all_conditions, reference_level))
} else {
  # Reference level not found, use alphabetical order
  warning("Reference level ", reference_level, " not found in data. Using alphabetical order for conditions.")
  condition_levels <- sort(all_conditions)
}

infection <- factor(conditions, levels=condition_levels)
replicate <- factor(replicates)

# Create a sample table
sample_table <- data.frame(
  sample = samples,
  infection = infection,
  replicate = replicate
)
rownames(sample_table) <- samples
print(sample_table)

# Define model matrix with infection as the only factor
cat("Creating model matrix for infection effect...\n")
design <- model.matrix(~infection, data=sample_table)
print(design)

# Calculate normalization factors
cat("Calculating normalization factors...\n")
y <- calcNormFactors(y)

# Log-transform counts for QC
logcounts <- cpm(y, log=TRUE)

# Save PCA plot
pdf(file.path(output_res_dir, "sample_pca_plot.pdf"), width=8, height=6)
pca <- prcomp(t(logcounts))
pca_data <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], 
                       Sample=colnames(logcounts),
                       Group=conditions)
plot(pca_data$PC1, pca_data$PC2, main=paste0("PCA of Samples (", resolution, "bp)"), 
     xlab=paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab=paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     col=as.factor(pca_data$Group), pch=19)
text(pca_data$PC1, pca_data$PC2, labels=pca_data$Sample, pos=3)
legend("topright", legend=levels(as.factor(pca_data$Group)), 
       col=1:length(levels(as.factor(pca_data$Group))), pch=19)
dev.off()

# Sample correlation heatmap
pdf(file.path(output_res_dir, "sample_correlation_heatmap.pdf"), width=7, height=6)
cor_matrix <- cor(logcounts)
pheatmap(cor_matrix, main=paste0("Sample Correlation Heatmap (", resolution, "bp)"),
         color=colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers=TRUE, number_format="%.3f")
dev.off()

# Estimate dispersion
cat("Estimating dispersion...\n")
y <- estimateDisp(y, design)

# Plot BCV to assess dispersion estimates
pdf(file.path(output_res_dir, "dispersion_BCV_plot.pdf"), width=7, height=6)
plotBCV(y)
title(paste0("Biological Coefficient of Variation (", resolution, "bp)"))
dev.off()

# Fit the full model (~infection)
cat("Fitting GLM model with infection effects...\n")
fit <- glmQLFit(y, design, BPPARAM=bpparam())

# ADD THIS SECTION TO GENERATE NULL MODEL
cat("Generating null model data...\n")

# Extract fitted values (expected counts under the model)
fitted_counts <- fit$fitted.values

# Get the null model parameters
null_model_data <- data.frame(
  interaction_id = rownames(y),
  
  # Expected counts under null hypothesis for each sample
  fitted_counts,
  
  # Dispersion estimates
  common_dispersion = y$common.dispersion,
  trended_dispersion = y$trended.dispersion,
  tagwise_dispersion = y$tagwise.dispersion,
  
  # Library size and normalization factors
  lib_sizes = paste(y$samples$lib.size, collapse = ";"),
  norm_factors = paste(y$samples$norm.factors, collapse = ";"),
  
  # Average log2 CPM (baseline expression level)
  avg_logCPM = aveLogCPM(y),
  
  # Add resolution info
  resolution = resolution
)

# Add coordinate information
if (!is.null(filtered_iset)) {
  tryCatch({
    annotations <- add_annotations(data.frame(row.names = rownames(y)), filtered_iset)
    null_model_data <- cbind(null_model_data, annotations[, c("chr1", "start1", "end1", "chr2", "start2", "end2", "interaction_type")])
  }, error = function(e) {
    cat("Could not add coordinate annotations to null model\n")
  })
}

# Save null model for this resolution
write.csv(null_model_data, 
          file.path(output_res_dir, "null_model_results.csv"), 
          row.names = TRUE)

cat("Saved null model with", nrow(null_model_data), "interactions\n")
# Test for differential interactions - get coefficients other than intercept
cat("Testing for differential interactions...\n")
coef_names <- colnames(design)[-1]  # Skip intercept
qlf_list <- list()

for (i in seq_along(coef_names)) {
  coef_name <- coef_names[i]
  cat("Testing coefficient:", coef_name, "\n")
  qlf_list[[coef_name]] <- glmQLFTest(fit, coef=i+1)
}
  
# Create the null model (intercept only)
cat("Creating and testing null model...\n")
null_design <- model.matrix(~1, data=sample_table)
null_fit <- glmFit(y, null_design)

# For null model, use deviance directly
null_dev <- sum(null_fit$deviance)
cat("Null model deviance:", null_dev, "\n")

# Perform model comparison
cat("Performing model comparison (full vs null)...\n")

# Calculate AIC and BIC for both models
full_model_aic <- -2 * sum(fit$loglik) + 2 * fit$df.prior
null_model_aic <- -2 * sum(null_fit$loglik) + 2
full_model_bic <- -2 * sum(fit$loglik) + log(ncol(y)) * fit$df.prior
null_model_bic <- -2 * sum(null_fit$loglik) + log(ncol(y))

# LRT test (using deviance difference)
lr_stat <- sum(null_fit$deviance) - sum(fit$deviance)
lr_df <- ncol(design) - 1
lr_pval <- 1 - pchisq(lr_stat, lr_df)

# Create a model comparison data frame
model_comparison <- data.frame(
  Model = c("Full Model (~infection)", "Null Model (Intercept only)"),
  LogLikelihood = c(sum(fit$loglik), sum(null_fit$loglik)),
  AIC = c(full_model_aic, null_model_aic),
  BIC = c(full_model_bic, null_model_bic),
  Deviance = c(sum(fit$deviance), sum(null_fit$deviance)),
  LRT_statistic = c(lr_stat, NA),
  LRT_df = c(lr_df, NA),
  LRT_pvalue = c(lr_pval, NA)
)

# Save model comparison
write.csv(model_comparison, file.path(output_res_dir, "model_comparison.csv"), row.names=FALSE)

# Create null results object
null_results <- list(
  table = data.frame(
    logFC = rep(0, nrow(y)),
    logCPM = aveLogCPM(y),
    PValue = rep(1, nrow(y)),
    FDR = rep(1, nrow(y)),
    row.names = rownames(y)
  )
)

# Write null results to file
write.csv(null_results$table, file.path(output_res_dir, "null_model_results.csv"), row.names=TRUE)

# Get the top differential interactions from each comparison
cat("Extracting results from coefficients...\n")
top_hits_list <- list()
for (coef_name in names(qlf_list)) {
  top_hits_list[[coef_name]] <- topTags(qlf_list[[coef_name]], n=Inf)
  
  # Write the results to file
  write.csv(top_hits_list[[coef_name]]$table, 
            file.path(output_res_dir, paste0("differential_interactions_", coef_name, ".csv")), 
            row.names=TRUE)
  
  # Also write a filtered version with only significant hits
  sig_hits <- top_hits_list[[coef_name]]$table[top_hits_list[[coef_name]]$table$FDR < fdr_threshold, ]
  if (nrow(sig_hits) > 0) {
    write.csv(sig_hits, 
              file.path(output_res_dir, paste0("significant_interactions_", coef_name, ".csv")), 
              row.names=TRUE)
  }
}

# Generate visualizations for each comparison
for (coef_name in names(qlf_list)) {
  # Get the topTags results for this coefficient
  qlf <- qlf_list[[coef_name]]
  top_hits <- top_hits_list[[coef_name]]
  
  # Filter by FDR
  significant <- top_hits$table[top_hits$table$FDR < fdr_threshold,]
  cat("Found", nrow(significant), "significant differential interactions for", 
      coef_name, "at FDR <", fdr_threshold, "\n")
  
  # Visualize the results
  pdf_file <- file.path(output_res_dir, paste0("differential_interactions_", coef_name, "_plots.pdf"))
  pdf(pdf_file, width=10, height=8)
  
  # MA plot
  plotMD(qlf)
  abline(h=c(-1, 0, 1), col=c("blue", "red", "blue"), lty=c(2,1,2))
  title(paste0("MA plot for ", coef_name, " (", resolution, "bp)"))
  
  # Volcano plot
  plot(top_hits$table$logFC, -log10(top_hits$table$PValue), 
       pch=20, col=ifelse(top_hits$table$FDR < fdr_threshold, "red", "black"),
       xlab="log2 Fold Change", ylab="-log10(P-value)",
       main=paste0("Volcano plot for ", coef_name, " (", resolution, "bp)"))
  abline(v=c(-1, 1), h=-log10(0.05), lty=2, col="blue")
  
  dev.off()
  
  # Add annotations to results
  cat("Adding annotations to results for", coef_name, "...\n")
  annotated_results <- add_annotations(top_hits$table, filtered_iset)
  annotated_file <- file.path(output_res_dir, paste0("annotated_differential_interactions_", coef_name, ".csv"))
  write.csv(annotated_results, annotated_file, row.names=TRUE)
  
  # Summarize differential interactions by type
  sig_results <- annotated_results[annotated_results$FDR < fdr_threshold,]
  
  if (nrow(sig_results) > 0) {
    # Count by interaction type
    type_summary <- table(sig_results$interaction_type)
    cat("\nSummary by interaction type for", coef_name, ":\n")
    print(type_summary)
    
    # Create a summary file
    summary_file <- file.path(output_res_dir, paste0("results_summary_", coef_name, ".txt"))
    sink(summary_file)
    cat("DIFFERENTIAL INTERACTION ANALYSIS SUMMARY FOR", coef_name, "\n")
    cat("=========================================\n")
    cat("Resolution:", resolution, "bp\n")
    cat("Total interactions analyzed:", nrow(filtered_iset), "\n")
    cat("Significant differential interactions (FDR <", fdr_threshold, "):", nrow(sig_results), "\n\n")
    
    cat("Summary by interaction type:\n")
    print(type_summary)
    cat("\n")
    
    # For cis interactions, summarize by chromosome
    if(sum(sig_results$interaction_type == "cis") > 0) {
      cis_summary <- table(sig_results$chr1[sig_results$interaction_type == "cis"])
      cat("Summary of cis interactions by chromosome:\n")
      print(cis_summary)
      cat("\n")
      
      # Summarize by distance
      if (any(!is.na(sig_results$interaction_distance[sig_results$interaction_type == "cis"]))) {
        dist_breaks <- c(0, 10000, 50000, 100000, 500000, 1000000, Inf)
        dist_labels <- c("<10kb", "10-50kb", "50-100kb", "100-500kb", "500kb-1Mb", ">1Mb")
        distance_summary <- table(cut(sig_results$interaction_distance[sig_results$interaction_type == "cis"], 
                                     breaks=dist_breaks, labels=dist_labels))
        cat("Summary of cis interactions by distance:\n")
        print(distance_summary)
        cat("\n")
      }
    }
    
    # For trans interactions, summarize by chromosome pairs
    if(sum(sig_results$interaction_type == "trans") > 0) {
      trans_pairs <- paste(sig_results$chr1, sig_results$chr2, sep="-")[sig_results$interaction_type == "trans"]
      trans_summary <- sort(table(trans_pairs), decreasing=TRUE)
      cat("Summary of trans interactions by chromosome pairs:\n")
      print(trans_summary)
    }
    
    sink()
    
    # Create supplementary files for significant interactions
    cat("Creating supplementary files for significant interactions for", coef_name, "...\n")
    
    # Extract significant cis and trans interactions
    sig_cis <- sig_results[sig_results$interaction_type == "cis", ]
    sig_trans <- sig_results[sig_results$interaction_type == "trans", ]
    
    # Save to separate files
    if (nrow(sig_cis) > 0) {
      write.csv(sig_cis, file.path(output_res_dir, paste0("significant_cis_interactions_", coef_name, ".csv")), 
                row.names=TRUE)
    }
    
    if (nrow(sig_trans) > 0) {
      write.csv(sig_trans, file.path(output_res_dir, paste0("significant_trans_interactions_", coef_name, ".csv")), 
                row.names=TRUE)
    }
    
    # Create a BED file for visualization in genome browsers
    bed_file <- file.path(output_res_dir, paste0("significant_interactions_", coef_name, ".bedpe"))
    write.table(
      data.frame(
        chr1 = sig_results$chr1,
        start1 = sig_results$start1,
        end1 = sig_results$end1,
        chr2 = sig_results$chr2,
        start2 = sig_results$start2,
        end2 = sig_results$end2,
        name = paste0("DI_", 1:nrow(sig_results)),
        score = -log10(sig_results$FDR),
        strand1 = ".",
        strand2 = ".",
        type = sig_results$interaction_type,
        logFC = sig_results$logFC
      ),
      bed_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE
    )
  }
}

#=============================================================================
# 4. Session Information and Completion
#=============================================================================

cat("\nDifferential analysis at", resolution, "bp resolution complete!\n")
cat("Results saved to:", output_res_dir, "\n")

# Session info for reproducibility
cat("\n=======================================================\n")
cat("Session Information\n")
cat("=======================================================\n")
print(sessionInfo())