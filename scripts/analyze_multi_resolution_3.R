#!/usr/bin/env Rscript
# As of 5/19/24 this script runs all the way with test dataset. 
# Multi-Resolution diffHic Analysis for Differential Chromatin Interactions
# 
# This script performs differential analysis of chromatin interactions between
# multiple Wolbachia infections at multiple resolutions using the diffHic package.
# It processes each resolution separately and compares the results across resolutions.
#
# Usage: Rscript analyze_multi_resolution_3.R --data_dir interactions_data --resolutions 1000,8000,32000,128000 
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
  make_option(c("--resolutions"), type="character", default="1000,8000,32000,128000",
              help="Comma-separated list of resolutions to analyze [default=%default]"),
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
                          description="Multi-resolution differential Micro-C analysis")
opt <- parse_args(opt_parser)

# Set parameters from command line arguments
data_dir <- opt$data_dir
resolutions <- as.numeric(unlist(strsplit(opt$resolutions, ",")))
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
cat("  Resolutions:", paste(resolutions, collapse=", "), "bp\n")
cat("  Min count:", min_count, "\n")
cat("  Min samples:", min_samples, "\n")
cat("  FDR threshold:", fdr_threshold, "\n")
cat("  Output directory:", output_dir, "\n")
cat("  Reference level:", reference_level, "\n")
cat("  Threads:", threads, "\n")

# Create output directories
dir.create(output_dir, showWarnings = FALSE)
for (res in resolutions) {
  dir.create(file.path(output_dir, paste0("res_", res)), showWarnings = FALSE)
}
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, showWarnings = FALSE)

#=============================================================================
# 2. Define Parameters
#=============================================================================

# # Define all resolutions to process
# resolutions <- c(1000, 8000, 32000, 128000) # This is hard coded, change this to the command line argument
# # resolutions <- c(128000)  #Test set

# # Set filtering parameters
# min_count <- 10  # Minimum read count to keep an interaction # This is hard coded, change this to the command line argument
# min_samples <- 2  # Minimum number of samples with min_count # This is hard coded, change this to the command line argument
# fdr_threshold <- 0.01  # FDR threshold for significant interactions # This is hard coded, change this to the command line argument

# # Set output directory
# output_dir <- "diffhic_results" # This is hard coded, change this to the command line argument
# dir.create(output_dir, showWarnings = FALSE)

# # Create directory for each resolution
# for (res in resolutions) {
#   dir.create(file.path(output_dir, paste0("res_", res)), showWarnings = FALSE)
# }

# # Create summary directory
# summary_dir <- file.path(output_dir, "summary")
# dir.create(summary_dir, showWarnings = FALSE)

#=============================================================================
# 3. Define Functions for Processing
#=============================================================================

# # Create directories for output
# dir.create(output_dir, showWarnings = FALSE)
# for (res in resolutions) {
#   dir.create(file.path(output_dir, paste0("res_", res)), showWarnings = FALSE)
# }
# summary_dir <- file.path(output_dir, "summary")
# dir.create(summary_dir, showWarnings = FALSE)

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

# Modified function to process a single resolution with DESeq2-like design
process_resolution <- function(resolution) {
  output_res_dir <- file.path(output_dir, paste0("res_", resolution))
  
  cat("\n=======================================================\n")
  cat("Processing resolution:", resolution, "bp\n")
  cat("=======================================================\n")
  
  # Load data from the combined files
  all_data_list <- read_all_data(resolution, data_dir)
  
  # Check if we have any data
  if (length(all_data_list) == 0) {
    warning("No data found for resolution ", resolution, "bp")
    return(NULL)
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
  
  # Test for differential interactions - get coefficients other than intercept
  cat("Testing for differential interactions...\n")
  coef_names <- colnames(design)[-1]  # Skip intercept
  qlf_list <- list()
  
  for (i in seq_along(coef_names)) {
    coef_name <- coef_names[i]
    cat("Testing coefficient:", coef_name, "\n")
    qlf_list[[coef_name]] <- glmQLFTest(fit, coef=i+1)  # Remove BPPARAM parameter
  }
    
  # Create the null model (intercept only)
  cat("Creating and testing null model...\n")
  null_design <- model.matrix(~1, data=sample_table)
  null_fit <- glmFit(y, null_design)  # Use glmFit instead of glmQLFit

  # For null model, use deviance directly
  null_dev <- sum(null_fit$deviance)
  cat("Null model deviance:", null_dev, "\n")

  # Perform model comparison
  cat("Performing model comparison (full vs null)...\n")

  # Calculate AIC and BIC for both models
  full_model_aic <- -2 * sum(fit$loglik) + 2 * fit$df.prior
  null_model_aic <- -2 * sum(null_fit$loglik) + 2  # Adjust for null model 
  full_model_bic <- -2 * sum(fit$loglik) + log(ncol(y)) * fit$df.prior
  null_model_bic <- -2 * sum(null_fit$loglik) + log(ncol(y))  # Adjust for null model

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

  # Create null results object (used for topTags in the original code)
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
  
  # Save null model results
  # null_results <- topTags(null_qlf, n=Inf)
  write.csv(null_results$table, file.path(output_res_dir, "null_model_results.csv"), row.names=TRUE)
  
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
  
  # Create a results object to return
  result <- list(
    resolution = resolution,
    significant = if (length(top_hits_list) > 0) {
      # Get the first coefficient's significant results
      top_hits_list[[1]]$table[top_hits_list[[1]]$table$FDR < fdr_threshold, ]
    } else NULL,
    all_interactions = if (length(top_hits_list) > 0) {
      add_annotations(top_hits_list[[1]]$table, filtered_iset)
    } else NULL,
    filtered_iset = filtered_iset,
    logcounts = logcounts,
    pca = pca,
    y = y,
    model_comparison = model_comparison,
    null_results = null_results$table,  # Updated to match your new structure
    all_results = top_hits_list
  )
  cat("\nDifferential analysis at", resolution, "bp resolution complete!\n")
  
  return(result)
}

#=============================================================================
# 4. Process Each Resolution
#=============================================================================

# Run the analysis for each resolution
results_list <- list()

for (res in resolutions) {
  results_list[[as.character(res)]] <- process_resolution(res)
}

#============================================================================= 
# 5. Compare Results Across Resolutions
#=============================================================================

cat("\n=======================================================\n")
cat("Comparing results across resolutions\n")
cat("=======================================================\n")

# Extract significant interactions from each resolution
sig_interactions <- list()
for (res in resolutions) {
  res_str <- as.character(res)
  
  if (!is.null(results_list[[res_str]]) && 
      !is.null(results_list[[res_str]]$significant) && 
      nrow(results_list[[res_str]]$significant) > 0) {
    # Create interaction IDs
    sig_results <- results_list[[res_str]]$significant
    sig_results$interaction_id <- create_interaction_id(
      sig_results$chr1, sig_results$start1, sig_results$end1,
      sig_results$chr2, sig_results$start2, sig_results$end2
    )
    sig_interactions[[res_str]] <- sig_results$interaction_id
  } else {
    sig_interactions[[res_str]] <- character(0)
  }
  
  cat("Resolution", res_str, "bp:", length(sig_interactions[[res_str]]), "significant interactions\n")
}

# Create Venn diagram to compare overlaps between resolutions
if (sum(sapply(sig_interactions, length)) > 0) {
  cat("Creating Venn diagrams to compare overlaps between resolutions...\n")
  
  # Create directory for Venn diagrams
  venn_dir <- file.path(summary_dir, "venn_diagrams")
  dir.create(venn_dir, showWarnings = FALSE)
  
  # Function to compare two resolutions
  compare_resolutions <- function(res1, res2) {
    res1_str <- as.character(res1)
    res2_str <- as.character(res2)
    
    if (length(sig_interactions[[res1_str]]) > 0 && length(sig_interactions[[res2_str]]) > 0) {
      overlap <- length(intersect(sig_interactions[[res1_str]], sig_interactions[[res2_str]]))
      
      cat("Overlap between", res1, "bp and", res2, "bp:", overlap, "interactions\n")
      
      # Create Venn diagram
      venn_file <- file.path(venn_dir, paste0("venn_", res1, "_vs_", res2, ".pdf"))
      pdf(venn_file, width=7, height=7)
      draw.pairwise.venn(
        area1 = length(sig_interactions[[res1_str]]),
        area2 = length(sig_interactions[[res2_str]]),
        cross.area = overlap,
        category = c(paste0(res1, "bp"), paste0(res2, "bp")),
        fill = c("blue", "red"),
        alpha = 0.5,
        cat.pos = c(0, 0),
        cat.dist = c(0.025, 0.025)
      )
      dev.off()
    }
  }
  
  # Compare all pairs of resolutions
  for (i in 1:(length(resolutions)-1)) {
    for (j in (i+1):length(resolutions)) {
      compare_resolutions(resolutions[i], resolutions[j])
    }
  }
}

#=============================================================================
# 6. Create Summary Reports and Visualizations
#=============================================================================

cat("\n=======================================================\n")
cat("Creating summary reports and visualizations\n")
cat("=======================================================\n")

# Create a summary table of results across all resolutions
summary_table <- data.frame(
  Resolution = resolutions,
  Total_Interactions = sapply(results_list, function(x) if(!is.null(x$filtered_iset)) nrow(x$filtered_iset) else 0),
  Significant_Interactions = sapply(results_list, function(x) if(!is.null(x$significant)) nrow(x$significant) else 0),
  Cis_Interactions = sapply(results_list, function(x) if(!is.null(x$significant)) sum(x$significant$interaction_type == "cis") else 0),
  Trans_Interactions = sapply(results_list, function(x) if(!is.null(x$significant)) sum(x$significant$interaction_type == "trans") else 0)
)

# Save summary table
write.csv(summary_table, file.path(summary_dir, "resolution_summary.csv"), row.names = FALSE)

# Create bar plot of significant interactions
pdf(file.path(summary_dir, "significant_interactions_barplot.pdf"), width=10, height=6)
par(mar=c(5, 5, 4, 5))
bp <- barplot(t(summary_table[,c("Cis_Interactions", "Trans_Interactions")]), 
              beside=TRUE, 
              names.arg=paste0(summary_table$Resolution, "bp"),
              col=c("darkblue", "darkred"),
              border="white",
              main="Significant Differential Interactions by Resolution",
              ylab="Number of Interactions",
              xlab="Resolution (bp)")
legend("topright", legend=c("Cis", "Trans"), fill=c("darkblue", "darkred"), border="white")

# Add total counts on top of bars
for (i in 1:nrow(summary_table)) {
  cis_x <- bp[1,i]
  trans_x <- bp[2,i]
  cis_y <- summary_table$Cis_Interactions[i]
  trans_y <- summary_table$Trans_Interactions[i]
  
  if (cis_y > 0) text(cis_x, cis_y + max(summary_table$Significant_Interactions) * 0.05, cis_y)
  if (trans_y > 0) text(trans_x, trans_y + max(summary_table$Significant_Interactions) * 0.05, trans_y)
}
dev.off()

# Create a comprehensive summary report
cat("Creating comprehensive summary report...\n")
sink(file.path(summary_dir, "comprehensive_summary.txt"))

cat("=====================================================\n")
cat("MULTI-RESOLUTION DIFFERENTIAL INTERACTION ANALYSIS\n")
cat("=====================================================\n\n")

cat("Analysis Parameters:\n")
cat("- Minimum count threshold:", min_count, "\n")
cat("- Minimum samples with minimum count:", min_samples, "\n")
cat("- FDR threshold for significance:", fdr_threshold, "\n\n")

cat("Summary of Results by Resolution:\n")
print(summary_table)
cat("\n")

# Report comparison to null:
cat("Model Comparison (Full vs Null Model):\n")
for (res in resolutions) {
  res_str <- as.character(res)
  if (!is.null(results_list[[res_str]]) && !is.null(results_list[[res_str]]$model_comparison)) {
    mc <- results_list[[res_str]]$model_comparison
    
    # Calculate LRT test result
    lr_stat <- 2 * (mc$LogLikelihood[1] - mc$LogLikelihood[2])
    lr_df <- mc$LRT_df[1]
    lr_pval <- 1 - pchisq(lr_stat, lr_df)
    
    cat(sprintf("- %sbp resolution:\n", res))
    cat(sprintf("  * Full model log-likelihood: %.2f\n", mc$LogLikelihood[1]))
    cat(sprintf("  * Null model log-likelihood: %.2f\n", mc$LogLikelihood[2]))
    cat(sprintf("  * LRT statistic: %.2f (p-value: %.2e)\n", lr_stat, lr_pval))
    cat(sprintf("  * Full model AIC: %.2f vs Null model AIC: %.2f\n", mc$AIC[1], mc$AIC[2]))
    cat(sprintf("  * Full model BIC: %.2f vs Null model BIC: %.2f\n", mc$BIC[1], mc$BIC[2]))
    
    # Determine if infection effect is significant
    if (lr_pval < 0.05) {
      cat("  * The infection effect is statistically significant (p < 0.05)\n")
    } else {
      cat("  * The infection effect is NOT statistically significant (p >= 0.05)\n")
    }
    cat("\n")
  }
}

# Report overlaps between resolutions
cat("Overlap of Significant Interactions Between Resolutions:\n")
for (i in 1:(length(resolutions)-1)) {
  for (j in (i+1):length(resolutions)) {
    res1 <- resolutions[i]
    res2 <- resolutions[j]
    res1_str <- as.character(res1)
    res2_str <- as.character(res2)
    
    if (length(sig_interactions[[res1_str]]) > 0 && length(sig_interactions[[res2_str]]) > 0) {
      overlap <- length(intersect(sig_interactions[[res1_str]], sig_interactions[[res2_str]]))
      
      cat(sprintf("- %sbp and %sbp: %d interactions (%.1f%% of %sbp, %.1f%% of %sbp)\n", 
                 res1, res2, overlap, 
                 100 * overlap / length(sig_interactions[[res1_str]]),
                 res1,
                 100 * overlap / length(sig_interactions[[res2_str]]),
                 res2))
    }
  }
}
cat("\n")

# Report interaction types by resolution
cat("Interaction Types by Resolution:\n")
for (res in resolutions) {
  res_str <- as.character(res)
  if (!is.null(results_list[[res_str]]) && !is.null(results_list[[res_str]]$significant)) {
    sig_results <- results_list[[res_str]]$significant
    
    if (nrow(sig_results) > 0) {
      type_counts <- table(sig_results$interaction_type)
      cat(sprintf("- %sbp: %d cis (%.1f%%), %d trans (%.1f%%)\n", 
                 res, 
                 type_counts["cis"], 100 * type_counts["cis"] / sum(type_counts),
                 type_counts["trans"], 100 * type_counts["trans"] / sum(type_counts)))
    }
  }
}
cat("\n")

# Report distance distribution for cis interactions
cat("Distance Distribution of Significant Cis Interactions:\n")
for (res in resolutions) {
  res_str <- as.character(res)
  if (!is.null(results_list[[res_str]]) && !is.null(results_list[[res_str]]$significant)) {
    sig_results <- results_list[[res_str]]$significant
    
    if (nrow(sig_results) > 0 && sum(sig_results$interaction_type == "cis") > 0) {
      cis_data <- sig_results[sig_results$interaction_type == "cis", ]
      
      # Create distance bins
      dist_breaks <- c(0, 10000, 50000, 100000, 500000, 1000000, Inf)
      dist_labels <- c("<10kb", "10-50kb", "50-100kb", "100-500kb", "500kb-1Mb", ">1Mb")
      
      distance_summary <- table(cut(cis_data$interaction_distance, 
                                   breaks=dist_breaks, labels=dist_labels))
      
      cat(sprintf("- %sbp resolution:\n", res))
      for (bin in names(distance_summary)) {
        cat(sprintf("  * %s: %d interactions (%.1f%%)\n", 
                   bin, distance_summary[bin], 
                   100 * distance_summary[bin] / sum(distance_summary)))
      }
      cat("\n")
    }
  }
}

cat("Conclusions and Recommendations:\n")
# Determine which resolution had the most interactions
most_sig_resolution <- resolutions[which.max(summary_table$Significant_Interactions)]
most_cis_resolution <- resolutions[which.max(summary_table$Cis_Interactions)]
most_trans_resolution <- resolutions[which.max(summary_table$Trans_Interactions)]

cat("- The highest number of significant interactions was found at", most_sig_resolution, "bp resolution.\n")
cat("- For cis interactions,", most_cis_resolution, "bp resolution showed the best detection.\n")
cat("- For trans interactions,", most_trans_resolution, "bp resolution was most effective.\n")

# Determine which resolution had the best representation of different distances
distance_representation <- character(0)
for (res in resolutions) {
  res_str <- as.character(res)
  if (!is.null(results_list[[res_str]]) && !is.null(results_list[[res_str]]$significant)) {
    sig_results <- results_list[[res_str]]$significant
    
    if (nrow(sig_results) > 0 && sum(sig_results$interaction_type == "cis") > 0) {
      cis_data <- sig_results[sig_results$interaction_type == "cis", ]
      
      # Create distance bins
      dist_breaks <- c(0, 10000, 50000, 100000, 500000, 1000000, Inf)
      dist_labels <- c("<10kb", "10-50kb", "50-100kb", "100-500kb", "500kb-1Mb", ">1Mb")
      
      distance_summary <- table(cut(cis_data$interaction_distance, 
                                   breaks=dist_breaks, labels=dist_labels))
      
      if (length(distance_summary) > 0 && sum(distance_summary > 0) >= 3) {
        distance_representation <- c(distance_representation, res_str)
      }
    }
  }
}

#####################################################################################
################### VOLCANO PLOTS ###################################################
######################################################################################

# Function to generate comprehensive interaction tables
generate_interaction_tables <- function(results_list, output_dir, fdr_threshold = 0.01) {
  cat("\nGenerating comprehensive interaction tables...\n")
  
  # Create directory for comprehensive tables
  tables_dir <- file.path(output_dir, "comprehensive_tables")
  dir.create(tables_dir, showWarnings = FALSE)
  
  # Combined table for all significant interactions across all resolutions
  all_sig_interactions <- NULL
  
  # Process each resolution
  for (res in names(results_list)) {
    if (!is.null(results_list[[res]]) && 
        !is.null(results_list[[res]]$significant) && 
        nrow(results_list[[res]]$significant) > 0) {
      
      # Get significant results for this resolution
      sig_results <- results_list[[res]]$significant
      
      # Add resolution column
      sig_results$resolution <- as.numeric(res)
      
      # Add to combined table
      all_sig_interactions <- rbind(all_sig_interactions, sig_results)
      
      # Save detailed table for this resolution
      write.csv(sig_results, 
                file.path(tables_dir, paste0("sig_interactions_", res, "bp.csv")), 
                row.names = FALSE)
      
      # Generate separate tables for cis and trans
      cis_results <- sig_results[sig_results$interaction_type == "cis", ]
      trans_results <- sig_results[sig_results$interaction_type == "trans", ]
      
      if (nrow(cis_results) > 0) {
        write.csv(cis_results, 
                  file.path(tables_dir, paste0("sig_cis_interactions_", res, "bp.csv")), 
                  row.names = FALSE)
      }
      
      if (nrow(trans_results) > 0) {
        write.csv(trans_results, 
                  file.path(tables_dir, paste0("sig_trans_interactions_", res, "bp.csv")), 
                  row.names = FALSE)
      }
    }
  }
  
  # Save combined table if there are any significant interactions
  if (!is.null(all_sig_interactions) && nrow(all_sig_interactions) > 0) {
    write.csv(all_sig_interactions, 
              file.path(tables_dir, "all_sig_interactions_combined.csv"), 
              row.names = FALSE)
    
    # Create separate tables for all cis and trans
    all_cis <- all_sig_interactions[all_sig_interactions$interaction_type == "cis", ]
    all_trans <- all_sig_interactions[all_sig_interactions$interaction_type == "trans", ]
    
    if (nrow(all_cis) > 0) {
      write.csv(all_cis, file.path(tables_dir, "all_sig_cis_interactions.csv"), row.names = FALSE)
    }
    
    if (nrow(all_trans) > 0) {
      write.csv(all_trans, file.path(tables_dir, "all_sig_trans_interactions.csv"), row.names = FALSE)
    }
  }
  
  cat("Completed generating comprehensive interaction tables.\n")
  cat("Tables saved to:", tables_dir, "\n")
}

# Function to create custom volcano plots using ggplot2 with fewer labels
create_custom_volcano_plot <- function(data, title, output_file, interaction_type = NULL, 
                                      logFC_col = "logFC", pval_col = "PValue", fdr_col = "FDR", 
                                      fdr_threshold = 0.01, logFC_threshold = 1) {
  
  # Ensure ggplot2 and ggrepel are loaded
  library(ggplot2)
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    install.packages("ggrepel")
  }
  library(ggrepel)
  
  # Filter data by interaction type if specified
  if (!is.null(interaction_type)) {
    plot_data <- data[data$interaction_type == interaction_type, ]
    subtitle <- paste0(" (", interaction_type, " interactions)")
  } else {
    plot_data <- data
    subtitle <- " (all interactions)"
  }
  
  # Skip if no data to plot
  if (nrow(plot_data) == 0) {
    warning(paste0("No data to plot for ", title, " ", interaction_type, " interactions"))
    return(NULL)
  }
  
  # Create custom labels for the points
  plot_data$label <- paste0(
    plot_data$chr1, ":", 
    format(plot_data$start1, scientific = FALSE), "-", 
    format(plot_data$end1, scientific = FALSE), "_",
    plot_data$chr2, ":", 
    format(plot_data$start2, scientific = FALSE), "-",
    format(plot_data$end2, scientific = FALSE)
  )
  
  # Calculate -log10 p-values
  plot_data$neg_log10_pval <- -log10(plot_data[[pval_col]])
  
  # Add a column for significance category
  plot_data$sig_cat <- "Not significant"
  plot_data$sig_cat[abs(plot_data[[logFC_col]]) > logFC_threshold] <- "Log2 fold change"
  plot_data$sig_cat[plot_data[[fdr_col]] < fdr_threshold] <- "P-value"
  plot_data$sig_cat[plot_data[[fdr_col]] < fdr_threshold & abs(plot_data[[logFC_col]]) > logFC_threshold] <- "Significant"
  
  # Add a column for color
  plot_data$color <- "grey80"  # Default (not significant)
  plot_data$color[abs(plot_data[[logFC_col]]) > logFC_threshold] <- "grey60"  # Only FC significant
  plot_data$color[plot_data[[fdr_col]] < fdr_threshold] <- "grey40"  # Only p-value significant
  
  # Color for significant points based on fold change direction
  sig_idx <- plot_data[[fdr_col]] < fdr_threshold & abs(plot_data[[logFC_col]]) > logFC_threshold
  plot_data$color[sig_idx & plot_data[[logFC_col]] > 0] <- "#8ecc85"  # UP in uninf - light green
  plot_data$color[sig_idx & plot_data[[logFC_col]] < 0] <- "#1bab4b"  # UP in wMel - dark green
  
  # MODIFIED: Use more stringent criteria for labeling
  # Label only the top N most significant points with the highest fold changes
  # This reduces the number of labels to prevent overcrowding
  
  # Filter to get only significant points
  sig_points <- plot_data[plot_data[[fdr_col]] < fdr_threshold & abs(plot_data[[logFC_col]]) > logFC_threshold, ]
  
  # By default, don't label any points
  plot_data$to_label <- FALSE
  
  if (nrow(sig_points) > 0) {
    # Sort by significance (most significant first) and then by fold change magnitude
    sig_points <- sig_points[order(sig_points[[fdr_col]], -abs(sig_points[[logFC_col]])), ]
    
    # Take only the top 15 points (adjust this number based on your preference)
    top_n <- min(15, nrow(sig_points))
    top_sig_idx <- rownames(sig_points)[1:top_n]
    
    # Set to_label flag for these points
    plot_data$to_label[rownames(plot_data) %in% top_sig_idx] <- TRUE
  }
  
  # Count significant points
  sig_uninf <- sum(sig_idx & plot_data[[logFC_col]] > 0)
  sig_wmel <- sum(sig_idx & plot_data[[logFC_col]] < 0)
  
  # Create the volcano plot
  p <- ggplot(plot_data, aes(x = .data[[logFC_col]], y = neg_log10_pval)) +
    # Add points
    geom_point(aes(color = color), size = 2, alpha = 0.75) +
    # Use color as is (don't map to a scale)
    scale_color_identity() +
    # Add labels for significant points
    ggrepel::geom_text_repel(
      data = subset(plot_data, to_label),
      aes(label = label),
      size = 2.5,
      box.padding = 0.5,
      point.padding = 0.2,
      segment.color = "grey50",
      max.overlaps = 15,  # Reduced from previous value
      min.segment.length = 0
    ) +
    # Add horizontal and vertical lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "grey") +
    # Add theme and labels
    theme_bw(base_size = 12) +
    labs(
      title = paste0(title, subtitle),
      subtitle = paste0("uninf. vs wMel, FDR < ", fdr_threshold),
      x = "log2 Fold Change",
      y = "-log10(P-value)",
      caption = paste0('Total = ', nrow(plot_data), ' interactions | Labeled top ', sum(plot_data$to_label), ' interactions')
    ) +
    # Add annotation for counts
    annotate(
      "text",
      x = max(plot_data[[logFC_col]], na.rm = TRUE) * 0.85,
      y = max(plot_data$neg_log10_pval, na.rm = TRUE) * 0.9,
      label = paste0("UP in uninf: ", sig_uninf, "\nUP in wMel: ", sig_wmel),
      hjust = 1,
      size = 3.5
    ) 
  
  # Save the plot to file
  pdf(output_file, width = 12, height = 10)
  print(p)
  dev.off()
  
  cat("Created custom volcano plot:", output_file, "\n")
  
  # Return the plot object
  return(p)
}

# Function to create volcano plots for all resolutions
create_custom_volcano_plots <- function(results_list, output_dir) {
  cat("\nCreating custom volcano plots...\n")
  
  # Create directory for volcano plots
  volcano_dir <- file.path(output_dir, "volcano_plots")
  dir.create(volcano_dir, showWarnings = FALSE)
  
  # Combine all results for overall volcano plot
  all_annotated_results <- NULL
  
  # Generate volcano plots for each resolution
  for (res in names(results_list)) {
    if (!is.null(results_list[[res]]) && 
        !is.null(results_list[[res]]$all_interactions) && 
        nrow(results_list[[res]]$all_interactions) > 0) {
      
      # Get all results and add resolution
      annotated_results <- results_list[[res]]$all_interactions
      annotated_results$resolution <- as.numeric(res)
      
      # Add to combined results
      all_annotated_results <- rbind(all_annotated_results, annotated_results)
      
      # Create overall volcano plot for this resolution
      create_custom_volcano_plot(
        data = annotated_results,
        title = paste0("Differential Interactions at ", res, "bp Resolution"),
        output_file = file.path(volcano_dir, paste0("custom_volcano_", res, "bp_all.pdf"))
      )
      
      # Create separate volcano plots for cis and trans
      create_custom_volcano_plot(
        data = annotated_results,
        title = paste0("Differential Interactions at ", res, "bp Resolution"),
        output_file = file.path(volcano_dir, paste0("custom_volcano_", res, "bp_cis.pdf")),
        interaction_type = "cis"
      )
      
      create_custom_volcano_plot(
        data = annotated_results,
        title = paste0("Differential Interactions at ", res, "bp Resolution"),
        output_file = file.path(volcano_dir, paste0("custom_volcano_", res, "bp_trans.pdf")),
        interaction_type = "trans"
      )
    }
  }
  
  # Create combined volcano plots if there are results
  if (!is.null(all_annotated_results) && nrow(all_annotated_results) > 0) {
    # Create overall combined volcano plot
    create_custom_volcano_plot(
      data = all_annotated_results,
      title = "Combined Differential Interactions Across All Resolutions",
      output_file = file.path(volcano_dir, "custom_volcano_combined_all.pdf")
    )
    
    # Create separate combined volcano plots for cis and trans
    create_custom_volcano_plot(
      data = all_annotated_results,
      title = "Combined Differential Interactions Across All Resolutions",
      output_file = file.path(volcano_dir, "custom_volcano_combined_cis.pdf"),
      interaction_type = "cis"
    )
    
    create_custom_volcano_plot(
      data = all_annotated_results,
      title = "Combined Differential Interactions Across All Resolutions",
      output_file = file.path(volcano_dir, "custom_volcano_combined_trans.pdf"),
      interaction_type = "trans"
    )
    
    # Create an additional plot with more labeled points for the top findings
    create_custom_volcano_plot(
      data = all_annotated_results,
      title = "Top Differential Interactions Across All Resolutions",
      output_file = file.path(volcano_dir, "custom_volcano_combined_top_hits.pdf"),
      fdr_threshold = 0.001,  # More stringent threshold to highlight top findings
      logFC_threshold = 1.5   # Higher log fold change threshold
    )
  }
  
  cat("Completed creating custom volcano plots.\n")
  cat("Plots saved to:", volcano_dir, "\n")
}

#=============================================================================
# Generate Enhanced Visualizations and Tables
#=============================================================================

cat("\n=======================================================\n")
cat("Generating enhanced visualizations and data tables\n")
cat("=======================================================\n")

# First make sure create_interaction_id function exists
if (!exists("create_interaction_id")) {
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
}

# Install and load required packages
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggplot2)
library(ggrepel)

# Create custom volcano plots for all resolutions
create_custom_volcano_plots(results_list, output_dir)

# Generate comprehensive interaction tables
generate_interaction_tables(results_list, output_dir, fdr_threshold)


########################################################################################
#=============================================================================
# 7. Session Information
#=============================================================================


if (length(distance_representation) > 0) {
  cat("- The following resolutions captured interactions across a broad range of distances:", 
      paste(distance_representation, "bp", collapse=", "), "\n")
}

cat("\nRecommended multi-resolution strategy:\n")
cat("1. Use", most_cis_resolution, "bp resolution for detailed analysis of enhancer-promoter interactions\n")
cat("2. Use", most_trans_resolution, "bp resolution for studying inter-chromosomal interactions\n")
cat("3. For a comprehensive view of the 3D chromatin landscape, combine results from all resolutions\n")

sink()

cat("\nMulti-resolution analysis complete! Results saved to", output_dir, "directory.\n")
cat("Summary reports and visualizations saved to", summary_dir, "directory.\n")
cat("=============================================================================\n")

cat("\n=======================================================\n")
cat("Session Information\n")        


# Session info for reproducibility
cat("Session Information:\n")
print(sessionInfo())
