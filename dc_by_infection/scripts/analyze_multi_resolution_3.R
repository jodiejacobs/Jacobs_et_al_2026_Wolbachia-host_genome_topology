#!/usr/bin/env Rscript
# Fixed Multi-Resolution diffHic Analysis for Differential Chromatin Interactions
# Fixes for data parsing, reference levels, and statistical issues

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
  make_option(c("--reference"), type="character", default="JW18DOX",
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

# FIXED: Improved data reading function with better error handling
read_all_data <- function(resolution, data_dir) {
  # List all contact files in the data directory
  contact_files <- list.files(data_dir, pattern = "_contacts\\.tsv$", full.names = TRUE)
  
  cat("Found contact files:", paste(basename(contact_files), collapse=", "), "\n")
  
  # Initialize a list to store data for each sample and replicate
  sample_data <- list()
  
  # Process each file
  for (file in contact_files) {
    cat("Processing file:", basename(file), "for resolution", resolution, "bp\n")
    
    # Extract sample name from filename
    filename <- basename(file)
    sample_name <- gsub("_contacts\\.tsv$", "", filename)
    
    # FIXED: More robust data reading with error handling
    tryCatch({
      # Try reading with different approaches
      data <- NULL
      
      # First, try reading all lines and check structure
      lines <- readLines(file, n = 10)  # Read first 10 lines to check structure
      cat("  First few lines of file:\n")
      cat("  ", head(lines, 3), sep = "\n  ")
      
      # Try fread with more robust settings
      data <- fread(file, header = TRUE, fill = TRUE, sep = "\t")
      
      cat("  Columns found:", ncol(data), "\n")
      cat("  Column names:", paste(colnames(data), collapse = ", "), "\n")
      
      # Check if data is empty after reading
      if (nrow(data) == 0) {
        warning("No data rows found in ", filename)
        next
      }
      
      # Ensure we have the minimum required columns
      required_base_cols <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
      
      # Check for alternative column names
      alt_col_mapping <- list(
        "chrom1" = c("chr1", "chromosome1", "seqnames1"),
        "chrom2" = c("chr2", "chromosome2", "seqnames2"),
        "count" = c("counts", "reads", "frequency", "score")
      )
      
      # Rename columns if needed
      for (std_name in names(alt_col_mapping)) {
        if (!std_name %in% colnames(data)) {
          alt_names <- alt_col_mapping[[std_name]]
          found_alt <- intersect(alt_names, colnames(data))
          if (length(found_alt) > 0) {
            setnames(data, found_alt[1], std_name)
            cat("  Renamed column:", found_alt[1], "->", std_name, "\n")
          }
        }
      }
      
      # Check for count column
      if (!"count" %in% colnames(data)) {
        # Look for any numeric column that could be count
        numeric_cols <- sapply(data, is.numeric)
        if (any(numeric_cols)) {
          count_col <- names(numeric_cols)[numeric_cols][1]  # Take first numeric column
          setnames(data, count_col, "count")
          cat("  Using column", count_col, "as count\n")
        } else {
          warning("No count/numeric column found in ", filename)
          next
        }
      }
      
      # Verify we have all required columns now
      required_cols <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "count")
      missing_cols <- setdiff(required_cols, colnames(data))
      if (length(missing_cols) > 0) {
        warning("Missing required columns in ", filename, ": ", paste(missing_cols, collapse=", "))
        next
      }
      
      # Filter by resolution if column exists
      if ("resolution" %in% colnames(data)) {
        data <- data[data$resolution == resolution, ]
      } else {
        cat("  No resolution column found, assuming all data is for", resolution, "bp\n")
      }
      
      # Apply minimum count filter
      data <- data[data$count >= min_count, ]
      
      if (nrow(data) == 0) {
        warning("No data found for resolution ", resolution, " in file ", filename, " after filtering")
        next
      }
      
      # Convert columns to numeric if they're not already
      cols_to_convert <- c("start1", "end1", "start2", "end2", "count")
      for (col in cols_to_convert) {
        if (col %in% colnames(data)) {
          if (is.character(data[[col]])) {
            data[[col]] <- as.numeric(gsub("[^0-9.-]", "", data[[col]]))
          }
        }
      }
      
      # Remove rows with NA values in essential columns
      data <- data[complete.cases(data[, required_cols, with = FALSE]), ]
      
      if (nrow(data) == 0) {
        warning("No valid data remaining after cleaning in ", filename)
        next
      }
      
      # Split data by existing condition and replicate columns
      condition_rep_combos <- unique(data[, c("condition", "replicate")])
      
      for (i in 1:nrow(condition_rep_combos)) {
        curr_condition <- condition_rep_combos$condition[i]
        curr_replicate <- condition_rep_combos$replicate[i]
        
        subset_data <- data[data$condition == curr_condition & data$replicate == curr_replicate, ]
        key <- paste0(curr_condition, "_rep", curr_replicate)
        sample_data[[key]] <- as.data.frame(subset_data)
        cat("  - Found", nrow(subset_data), "interactions for", key, "\n")
      }
      
    }, error = function(e) {
      warning("Error processing file ", filename, ": ", e$message)
    })
  }
  
  return(sample_data)
}

# FIXED: More robust InteractionSet creation
find_all_interactions <- function(iset_list) {
  cat("Finding all unique interactions across samples...\n")
  
  if (length(iset_list) == 0) {
    warning("No interaction sets provided")
    return(NULL)
  }
  
  # Extract all interactions as data frames
  interactions_data <- list()
  
  for (i in seq_along(iset_list)) {
    sample_name <- names(iset_list)[i]
    cat("Processing", sample_name, "\n")
    
    current_iset <- iset_list[[i]]
    
    # Handle case where iset might be a data frame instead of InteractionSet
    if (is.data.frame(current_iset)) {
      int_df <- current_iset
      int_df$sample <- i
    } else {
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
    }
    
    interactions_data[[sample_name]] <- int_df
  }
  
  # Combine all interactions
  all_ints <- do.call(rbind, interactions_data)
  
  if (nrow(all_ints) == 0) {
    warning("No interactions found across all samples")
    return(NULL)
  }
  
  # Create a unique identifier for each interaction
  all_ints$int_id <- create_interaction_id(
    all_ints$chrom1, all_ints$start1, all_ints$end1,
    all_ints$chrom2, all_ints$start2, all_ints$end2
  )
  
  # Find unique interactions
  unique_ints <- unique(all_ints[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "int_id")])
  cat("Found", nrow(unique_ints), "unique interactions\n")
  
  # Create count matrix
  n_unique <- nrow(unique_ints)
  n_samples <- length(iset_list)
  count_matrix <- matrix(0, nrow=n_unique, ncol=n_samples)
  colnames(count_matrix) <- names(iset_list)
  
  # Create a map from int_id to row index
  int_id_to_index <- setNames(1:nrow(unique_ints), unique_ints$int_id)
  
  cat("Filling count matrix...\n")
  
  # Fill the count matrix
  for (j in 1:n_samples) {
    sample_data <- all_ints[all_ints$sample == j, ]
    
    if (nrow(sample_data) > 0) {
      # Look up the row indices
      row_indices <- int_id_to_index[sample_data$int_id]
      
      # Aggregate counts for duplicate interactions
      for (k in seq_along(row_indices)) {
        if (!is.na(row_indices[k])) {
          count_matrix[row_indices[k], j] <- count_matrix[row_indices[k], j] + sample_data$count[k]
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

# FIXED: Add annotations function with better error handling
add_annotations <- function(results_table, iset) {
  tryCatch({
    # Get the interactions
    gi <- interactions(iset)
    
    # Get anchors
    anchor1 <- anchors(gi, type="first")
    anchor2 <- anchors(gi, type="second")
    
    # Determine if each interaction is cis or trans
    chr1 <- as.character(seqnames(anchor1))
    chr2 <- as.character(seqnames(anchor2))
    is_cis <- chr1 == chr2
    
    # Calculate interaction distance for cis interactions
    interaction_distance <- rep(NA, length(gi))
    if (any(is_cis)) {
      interaction_distance[is_cis] <- end(anchor2)[is_cis] - start(anchor1)[is_cis]
    }
    
    # Create a data frame with interaction info
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
    
  }, error = function(e) {
    warning("Error adding annotations: ", e$message)
    return(results_table)
  })
}

# FIXED: Modified process_resolution function with better error handling
process_resolution <- function(resolution) {
  output_res_dir <- file.path(output_dir, paste0("res_", resolution))
  
  cat("\n=======================================================\n")
  cat("Processing resolution:", resolution, "bp\n")
  cat("=======================================================\n")
  
  tryCatch({
    # Load data from the combined files
    all_data_list <- read_all_data(resolution, data_dir)
    
    # Check if we have any data
    if (length(all_data_list) == 0) {
      warning("No data found for resolution ", resolution, "bp")
      return(NULL)
    }
    
    cat("Found", length(all_data_list), "samples with data\n")
    
    # Create InteractionSet objects for each sample
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
        
        # Create InteractionSet object
        iset_list[[sample_name]] <- InteractionSet(list(counts=matrix(data$count, ncol=1)), gint)
      }
    }
    
    if (length(iset_list) == 0) {
      warning("No valid InteractionSet objects created for resolution ", resolution, "bp")
      return(NULL)
    }
    
    # Create unified InteractionSet with all samples
    cat("Creating combined interaction set...\n")
    combined_iset <- find_all_interactions(iset_list)
    
    if (is.null(combined_iset) || nrow(combined_iset) == 0) {
      warning("Failed to create combined interaction set for resolution ", resolution, "bp")
      return(NULL)
    }
    
    # Calculate library sizes
    lib_sizes <- colSums(assay(combined_iset))
    combined_iset$totals <- lib_sizes
    
    cat("Combined interaction set has", nrow(combined_iset), "interactions across", 
        ncol(assay(combined_iset)), "samples\n")
    
    # Filter out low-abundance interactions
    cat("Filtering low-abundance interactions...\n")
    keep <- rowSums(assay(combined_iset) >= min_count) >= min_samples
    filtered_iset <- combined_iset[keep,]
    
    if (nrow(filtered_iset) == 0) {
      warning("No interactions remain after filtering for resolution ", resolution, "bp")
      return(NULL)
    }
    
    cat("Filtered from", nrow(combined_iset), "to", nrow(filtered_iset), "interactions\n")
    
    # Create DGEList object
    cat("Creating DGEList object for analysis...\n")
    y <- asDGEList(filtered_iset)
    
    # Extract sample information from column names
    samples <- colnames(y)
    
    # FIXED: Better condition parsing
    conditions <- sapply(strsplit(samples, "_rep"), function(x) x[1])
    replicates <- sapply(strsplit(samples, "_rep"), function(x) x[2])
    
    # Create condition factor with proper reference level
    all_conditions <- unique(conditions)
    cat("Found conditions:", paste(all_conditions, collapse=", "), "\n")
    
    if (reference_level %in% all_conditions) {
      condition_levels <- c(reference_level, setdiff(all_conditions, reference_level))
    } else {
      warning("Reference level ", reference_level, " not found. Available conditions: ", 
              paste(all_conditions, collapse=", "))
      condition_levels <- sort(all_conditions)
      reference_level <- condition_levels[1]
      cat("Using", reference_level, "as reference level\n")
    }
    
    infection <- factor(conditions, levels=condition_levels)
    replicate <- factor(replicates)
    
    # Create sample table
    sample_table <- data.frame(
      sample = samples,
      infection = infection,
      replicate = replicate
    )
    rownames(sample_table) <- samples
    print(sample_table)
    
    # FIXED: Check for sufficient samples per condition
    condition_counts <- table(infection)
    if (any(condition_counts < 2)) {
      warning("Some conditions have fewer than 2 replicates. This may cause issues.")
      cat("Condition counts:\n")
      print(condition_counts)
    }
    
    # Define model matrix
    cat("Creating model matrix for infection effect...\n")
    design <- model.matrix(~infection, data=sample_table)
    print(design)
    
    # Calculate normalization factors
    cat("Calculating normalization factors...\n")
    y <- calcNormFactors(y)
    
    # FIXED: Better handling of log-transform with error checking
    cat("Computing log-transformed counts for QC...\n")
    
    # Check for zero library sizes
    if (any(y$samples$lib.size == 0)) {
      warning("Some samples have zero library sizes")
      y$samples$lib.size[y$samples$lib.size == 0] <- 1  # Set minimum lib size
    }
    
    logcounts <- cpm(y, log=TRUE)
    
    # Check for issues with logcounts
    if (any(!is.finite(logcounts))) {
      warning("Non-finite values in log counts. Replacing with minimum values.")
      logcounts[!is.finite(logcounts)] <- min(logcounts[is.finite(logcounts)], na.rm = TRUE)
    }
    
    # FIXED: More robust PCA and correlation analysis
    cat("Creating QC plots...\n")
    
    # PCA plot with error handling
    tryCatch({
      pdf(file.path(output_res_dir, "sample_pca_plot.pdf"), width=8, height=6)
      
      # Check if we have enough variation for PCA
      logcounts_var <- apply(logcounts, 1, var, na.rm = TRUE)
      high_var_genes <- logcounts[logcounts_var > 0 & is.finite(logcounts_var), ]
      
      if (nrow(high_var_genes) > 10) {
        pca <- prcomp(t(high_var_genes), center = TRUE, scale. = FALSE)
        pca_data <- data.frame(
          PC1 = pca$x[,1], 
          PC2 = pca$x[,2], 
          Sample = colnames(logcounts),
          Group = conditions
        )
        
        plot(pca_data$PC1, pca_data$PC2, 
             main = paste0("PCA of Samples (", resolution, "bp)"), 
             xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
             ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
             col = as.factor(pca_data$Group), pch = 19)
        text(pca_data$PC1, pca_data$PC2, labels = pca_data$Sample, pos = 3, cex = 0.7)
        legend("topright", legend = levels(as.factor(pca_data$Group)), 
               col = 1:length(levels(as.factor(pca_data$Group))), pch = 19)
      } else {
        plot(1, 1, main = "Insufficient variation for PCA", type = "n")
        text(1, 1, "Not enough variable interactions\nfor meaningful PCA")
      }
      
      dev.off()
    }, error = function(e) {
      cat("Error creating PCA plot:", e$message, "\n")
      if (dev.cur() > 1) dev.off()
    })
    
    # FIXED: More robust correlation heatmap
    tryCatch({
      pdf(file.path(output_res_dir, "sample_correlation_heatmap.pdf"), width=7, height=6)
      
      # Check if correlation matrix can be computed
      cor_matrix <- cor(logcounts, use = "complete.obs")
      
      # Check for issues with correlation matrix
      if (any(!is.finite(cor_matrix))) {
        warning("Non-finite values in correlation matrix")
        cor_matrix[!is.finite(cor_matrix)] <- 0
        diag(cor_matrix) <- 1
      }
      
      # Only create heatmap if we have valid correlations
      if (all(is.finite(cor_matrix)) && nrow(cor_matrix) > 1) {
        pheatmap(cor_matrix, 
                main = paste0("Sample Correlation Heatmap (", resolution, "bp)"),
                color = colorRampPalette(c("blue", "white", "red"))(100),
                display_numbers = TRUE, 
                number_format = "%.3f",
                cluster_rows = FALSE,  # Disable clustering to avoid hclust error
                cluster_cols = FALSE)
      } else {
        plot(1, 1, main = "Cannot compute sample correlations", type = "n")
        text(1, 1, "Insufficient data for\ncorrelation analysis")
      }
      
      dev.off()
    }, error = function(e) {
      cat("Error creating correlation heatmap:", e$message, "\n")
      if (dev.cur() > 1) dev.off()
    })
    
    # Continue with statistical analysis
    cat("Estimating dispersion...\n")
    y <- estimateDisp(y, design)
    
    # Plot BCV
    tryCatch({
      pdf(file.path(output_res_dir, "dispersion_BCV_plot.pdf"), width=7, height=6)
      plotBCV(y)
      title(paste0("Biological Coefficient of Variation (", resolution, "bp)"))
      dev.off()
    }, error = function(e) {
      cat("Error creating BCV plot:", e$message, "\n")
      if (dev.cur() > 1) dev.off()
    })
    
    # Fit the model
    cat("Fitting GLM model...\n")
    fit <- glmQLFit(y, design, BPPARAM=bpparam())
    
    # Test for differential interactions
    cat("Testing for differential interactions...\n")
    coef_names <- colnames(design)[-1]  # Skip intercept
    qlf_list <- list()
    
    for (i in seq_along(coef_names)) {
      coef_name <- coef_names[i]
      cat("Testing coefficient:", coef_name, "\n")
      qlf_list[[coef_name]] <- glmQLFTest(fit, coef=i+1)
    }
    
    # Extract results
    top_hits_list <- list()
    for (coef_name in names(qlf_list)) {
      top_hits_list[[coef_name]] <- topTags(qlf_list[[coef_name]], n=Inf)
      
      # Write results
      write.csv(top_hits_list[[coef_name]]$table, 
                file.path(output_res_dir, paste0("differential_interactions_", coef_name, ".csv")), 
                row.names=TRUE)
      
      # Filter significant hits
      sig_hits <- top_hits_list[[coef_name]]$table[top_hits_list[[coef_name]]$table$FDR < fdr_threshold, ]
      if (nrow(sig_hits) > 0) {
        write.csv(sig_hits, 
                  file.path(output_res_dir, paste0("significant_interactions_", coef_name, ".csv")), 
                  row.names=TRUE)
      }
      
      cat("Found", nrow(sig_hits), "significant interactions for", coef_name, "\n")
    }
    
    # Create visualizations
    for (coef_name in names(qlf_list)) {
      qlf <- qlf_list[[coef_name]]
      top_hits <- top_hits_list[[coef_name]]
      
      # Create plots
      tryCatch({
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
      }, error = function(e) {
        cat("Error creating plots for", coef_name, ":", e$message, "\n")
        if (dev.cur() > 1) dev.off()
      })
      
      # Add annotations
      tryCatch({
        annotated_results <- add_annotations(top_hits$table, filtered_iset)
        annotated_file <- file.path(output_res_dir, paste0("annotated_differential_interactions_", coef_name, ".csv"))
        write.csv(annotated_results, annotated_file, row.names=TRUE)
      }, error = function(e) {
        cat("Error adding annotations for", coef_name, ":", e$message, "\n")
      })
    }
    
    # Create results object
    result <- list(
      resolution = resolution,
      significant = if (length(top_hits_list) > 0) {
        sig_results <- top_hits_list[[1]]$table[top_hits_list[[1]]$table$FDR < fdr_threshold, ]
        if (nrow(sig_results) > 0) {
          tryCatch({
            add_annotations(sig_results, filtered_iset)
          }, error = function(e) {
            sig_results
          })
        } else {
          NULL
        }
      } else NULL,
      all_interactions = if (length(top_hits_list) > 0) {
        tryCatch({
          add_annotations(top_hits_list[[1]]$table, filtered_iset)
        }, error = function(e) {
          top_hits_list[[1]]$table
        })
      } else NULL,
      filtered_iset = filtered_iset,
      logcounts = logcounts,
      y = y,
      all_results = top_hits_list
    )
    
    cat("\nDifferential analysis at", resolution, "bp resolution complete!\n")
    return(result)
    
  }, error = function(e) {
    cat("Error processing resolution", resolution, "bp:", e$message, "\n")
    return(NULL)
  })
}

#=============================================================================
# 3. Process Each Resolution
#=============================================================================

# Run the analysis for each resolution
results_list <- list()

for (res in resolutions) {
  result <- process_resolution(res)
  if (!is.null(result)) {
    results_list[[as.character(res)]] <- result
  }
}

# Check if we have any results
if (length(results_list) == 0) {
  stop("No successful analyses completed. Check your data files and parameters.")
}

cat("\nCompleted analysis for", length(results_list), "resolutions.\n")

#=============================================================================
# 4. Generate Summary
#=============================================================================

cat("\n=======================================================\n")
cat("Generating summary report\n")
cat("=======================================================\n")

# Create summary table
summary_data <- data.frame(
  Resolution = character(0),
  Total_Interactions = numeric(0),
  Significant_Interactions = numeric(0),
  stringsAsFactors = FALSE
)

for (res_name in names(results_list)) {
  result <- results_list[[res_name]]
  
  total_int <- if (!is.null(result$filtered_iset)) nrow(result$filtered_iset) else 0
  sig_int <- if (!is.null(result$significant)) nrow(result$significant) else 0
  
  summary_data <- rbind(summary_data, data.frame(
    Resolution = res_name,
    Total_Interactions = total_int,
    Significant_Interactions = sig_int,
    stringsAsFactors = FALSE
  ))
}

# Save summary
write.csv(summary_data, file.path(summary_dir, "analysis_summary.csv"), row.names = FALSE)

# Print summary
cat("\nAnalysis Summary:\n")
print(summary_data)

cat("\nAnalysis complete! Check the output directory for results:\n")
cat("Main results:", output_dir, "\n")
cat("Summary:", summary_dir, "\n")

# Session info
cat("\nSession Information:\n")
print(sessionInfo())
