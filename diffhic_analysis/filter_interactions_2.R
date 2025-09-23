#!/usr/bin/env Rscript

# Script to filter differential chromatin interactions and prepare for downstream analysis
# with cdBEST, REDfly, and DroID
# Modified to include comparison between dataset and null_dataset

# Load required libraries
library(data.table)
library(dplyr)

# Set parameters 
FDR_THRESHOLD <- 0.05  # FDR threshold for significant interactions
LOG_FC_THRESHOLD <- 1  # Log2 fold change threshold (2x fold change)

# Function to read the interaction data
read_interaction_data <- function(file_path) {
  # Read the data
  cat("Reading interaction data from:", file_path, "\n")
  interactions <- fread(file_path)
  cat("Total interactions loaded:", nrow(interactions), "\n")
  return(interactions)
}

# Function to compare dataset with null_dataset
compare_with_null <- function(dataset, null_dataset) {
  cat("Comparing dataset with null dataset...\n")
  
  # Ensure both datasets have the same column structure
  common_cols <- intersect(colnames(dataset), colnames(null_dataset))
  if(length(common_cols) < ncol(dataset)) {
    cat("Warning: Some columns in the dataset are not present in the null dataset.\n")
    cat("Proceeding with common columns only.\n")
  }
  
  # Create a unique identifier for interactions in both datasets
  # Using chromosome and position information
  dataset$interaction_id <- paste(
    dataset$chr1, dataset$start1, dataset$end1,
    dataset$chr2, dataset$start2, dataset$end2,
    sep = "_"
  )
  
  null_dataset$interaction_id <- paste(
    null_dataset$chr1, null_dataset$start1, null_dataset$end1,
    null_dataset$chr2, null_dataset$start2, null_dataset$end2,
    sep = "_"
  )
  
  # Find interactions present in both datasets
  common_interactions <- intersect(dataset$interaction_id, null_dataset$interaction_id)
  cat("Common interactions found in both datasets:", length(common_interactions), "\n")
  
  # Get the subset of interactions in both datasets
  dataset_common <- dataset[dataset$interaction_id %in% common_interactions, ]
  null_dataset_common <- null_dataset[null_dataset$interaction_id %in% common_interactions, ]
  
  # For matching interaction IDs, calculate the empirical p-value
  # by comparing the real value to the null distribution
  dataset_common <- dataset_common[order(dataset_common$interaction_id), ]
  null_dataset_common <- null_dataset_common[order(null_dataset_common$interaction_id), ]
  
  # Calculate empirical p-value for each common interaction
  empirical_pvals <- numeric(length(common_interactions))
  for(i in 1:length(common_interactions)) {
    # Assuming the interactionStrength column exists
    # If not, replace with the appropriate column
    if("interactionStrength" %in% colnames(dataset_common) && 
       "interactionStrength" %in% colnames(null_dataset_common)) {
      real_val <- dataset_common$interactionStrength[i]
      null_val <- null_dataset_common$interactionStrength[i]
      empirical_pvals[i] <- ifelse(real_val > null_val, 0.5, 1)
    } else if("logFC" %in% colnames(dataset_common)) {
      # Alternatively use logFC if available
      real_val <- abs(dataset_common$logFC[i])
      null_val <- ifelse("logFC" %in% colnames(null_dataset_common), 
                         abs(null_dataset_common$logFC[i]), 0)
      empirical_pvals[i] <- ifelse(real_val > null_val, 0.5, 1)
    } else {
      # If neither column exists, set p-value to NA
      empirical_pvals[i] <- NA
    }
  }
  
  # Add empirical p-values to the dataset
  dataset_common$empirical_pval <- empirical_pvals
  
  # Find interactions unique to the dataset
  unique_interactions <- setdiff(dataset$interaction_id, null_dataset$interaction_id)
  cat("Unique interactions only in main dataset:", length(unique_interactions), "\n")
  
  # Get the subset of interactions only in the dataset
  dataset_unique <- dataset[dataset$interaction_id %in% unique_interactions, ]
  dataset_unique$empirical_pval <- 0.001  # Assign a low p-value to unique interactions
  
  # Combine common and unique interactions
  dataset_with_empirical <- rbind(dataset_common, dataset_unique)
  
  # Adjust p-values for multiple testing
  dataset_with_empirical$empirical_FDR <- p.adjust(dataset_with_empirical$empirical_pval, method = "BH")
  
  # Clean up the temporary identifier
  dataset_with_empirical$interaction_id <- NULL
  
  cat("Added empirical p-values based on null dataset comparison.\n")
  return(dataset_with_empirical)
}

# Function to filter significant interactions
filter_significant_interactions <- function(interactions, fdr_threshold = FDR_THRESHOLD, logfc_threshold = LOG_FC_THRESHOLD) {
  # Filter by FDR and logFC
  sig_interactions <- interactions[interactions$FDR < fdr_threshold & abs(interactions$logFC) > logfc_threshold, ]
  cat("Significant interactions (FDR <", fdr_threshold, "and |logFC| >", logfc_threshold, "):", nrow(sig_interactions), "\n")
  
  # Also filter by empirical FDR if available
  if("empirical_FDR" %in% colnames(interactions)) {
    emp_sig_interactions <- interactions[interactions$empirical_FDR < fdr_threshold & abs(interactions$logFC) > logfc_threshold, ]
    cat("Significant interactions by empirical FDR <", fdr_threshold, "and |logFC| >", logfc_threshold, "):", nrow(emp_sig_interactions), "\n")
    
    # Intersect the two sets of significant interactions
    combined_sig <- intersect(
      interactions$interaction_id[interactions$FDR < fdr_threshold & abs(interactions$logFC) > logfc_threshold],
      interactions$interaction_id[interactions$empirical_FDR < fdr_threshold & abs(interactions$logFC) > logfc_threshold]
    )
    if(length(combined_sig) > 0) {
      combined_sig_interactions <- interactions[interactions$interaction_id %in% combined_sig, ]
      cat("Significant interactions by both FDR methods:", nrow(combined_sig_interactions), "\n")
      return(combined_sig_interactions)
    }
    
    # If intersection is empty, return original significant interactions
    return(sig_interactions)
  }
  
  return(sig_interactions)
}

# Function to separate cis and trans interactions
separate_cis_trans <- function(sig_interactions) {
  # Separate by interaction_type
  cis_interactions <- sig_interactions[sig_interactions$interaction_type == "cis", ]
  trans_interactions <- sig_interactions[sig_interactions$interaction_type == "trans", ]
  
  cat("Cis interactions:", nrow(cis_interactions), "\n")
  cat("Trans interactions:", nrow(trans_interactions), "\n")
  
  return(list(cis = cis_interactions, trans = trans_interactions))
}

# Function to sort by significance and fold change
sort_by_significance <- function(interactions) {
  # First sort by FDR, then by absolute logFC
  sorted <- interactions[order(interactions$FDR, -abs(interactions$logFC)), ]
  return(sorted)
}

# Function to prepare files for cdBEST analysis
prepare_cdbest_input <- function(interactions, output_prefix) {
  # For each interaction, create a fasta file with the sequence
  cat("Preparing cdBEST input files...\n")
  
  # Create a BED file with the interactions for easier extraction of sequences later
  bed_file <- paste0(output_prefix, "_regions.bed")
  
  # Write regions for anchor1 and anchor2
  regions_bed <- rbind(
    data.frame(
      chrom = interactions$chr1,
      start = interactions$start1,
      end = interactions$end1,
      name = paste0("anchor1_", 1:nrow(interactions)),
      score = -log10(interactions$FDR),
      strand = "."
    ),
    data.frame(
      chrom = interactions$chr2,
      start = interactions$start2,
      end = interactions$end2,
      name = paste0("anchor2_", 1:nrow(interactions)),
      score = -log10(interactions$FDR),
      strand = "."
    )
  )
  
  # Sort regions by chromosome and start position
  regions_bed <- regions_bed[order(regions_bed$chrom, regions_bed$start), ]
  
  # Write the BED file
  write.table(regions_bed, file = bed_file, quote = FALSE, sep = "\t", 
              row.names = FALSE, col.names = FALSE)
  
  cat("Wrote BED file for sequence extraction:", bed_file, "\n")
  cat("Note: You'll need to extract sequences using bedtools getfasta:\n")
  cat("bedtools getfasta -fi reference.fa -bed", bed_file, "-fo", paste0(output_prefix, "_sequences.fa"), "\n")
  
  return(bed_file)
}

# Function to prepare files for REDfly analysis (for cis interactions)
prepare_redfly_input <- function(cis_interactions, output_prefix) {
  cat("Preparing REDfly input files...\n")
  
  # Create a BED file with the cis interactions
  redfly_bed <- paste0(output_prefix, "_redfly.bed")
  
  # Prepare BED format
  redfly_data <- data.frame(
    chrom = cis_interactions$chr1,  # Using chr1 since this is cis
    start = pmin(cis_interactions$start1, cis_interactions$start2),  # Min start
    end = pmax(cis_interactions$end1, cis_interactions$end2),  # Max end
    name = paste0("interaction_", 1:nrow(cis_interactions)),
    score = -log10(cis_interactions$FDR),
    strand = "."
  )
  
  # Write the BED file
  write.table(redfly_data, file = redfly_bed, quote = FALSE, sep = "\t", 
              row.names = FALSE, col.names = FALSE)
  
  cat("Wrote BED file for REDfly analysis:", redfly_bed, "\n")
  cat("Note: Use this file with bedtools intersect to find overlaps with REDfly elements:\n")
  cat("bedtools intersect -a", redfly_bed, "-b redfly_annotations.gff -wa -wb > redfly_overlaps.txt\n")
  
  return(redfly_bed)
}

# Function to prepare files for DroID analysis (for trans interactions)
prepare_droid_input <- function(trans_interactions, output_prefix) {
  cat("Preparing DroID input files...\n")
  
  # Create two BED files - one for each anchor
  droid_bed1 <- paste0(output_prefix, "_droid_anchor1.bed")
  droid_bed2 <- paste0(output_prefix, "_droid_anchor2.bed")
  
  # Prepare BED format for first anchor
  droid_data1 <- data.frame(
    chrom = trans_interactions$chr1,
    start = trans_interactions$start1,
    end = trans_interactions$end1,
    name = paste0("interaction_", 1:nrow(trans_interactions), "_anchor1"),
    score = -log10(trans_interactions$FDR),
    strand = "."
  )
  
  # Prepare BED format for second anchor
  droid_data2 <- data.frame(
    chrom = trans_interactions$chr2,
    start = trans_interactions$start2,
    end = trans_interactions$end2,
    name = paste0("interaction_", 1:nrow(trans_interactions), "_anchor2"),
    score = -log10(trans_interactions$FDR),
    strand = "."
  )
  
  # Write the BED files
  write.table(droid_data1, file = droid_bed1, quote = FALSE, sep = "\t", 
              row.names = FALSE, col.names = FALSE)
  write.table(droid_data2, file = droid_bed2, quote = FALSE, sep = "\t", 
              row.names = FALSE, col.names = FALSE)
  
  cat("Wrote BED files for DroID analysis:\n")
  cat(" - Anchor 1:", droid_bed1, "\n")
  cat(" - Anchor 2:", droid_bed2, "\n")
  cat("Note: Use these files to look up genes in the regions and then query DroID for interactions\n")
  
  return(list(anchor1 = droid_bed1, anchor2 = droid_bed2))
}

# Function to create summary files with formatted output
create_summary_files <- function(cis_interactions, trans_interactions, output_prefix, include_empirical = FALSE) {
  # Create summary files for top interactions
  cat("Creating summary files...\n")
  
  # Top cis interactions
  top_cis_file <- paste0(output_prefix, "_top_cis_interactions.txt")
  if (nrow(cis_interactions) > 0) {
    top_cis <- head(cis_interactions, 100)  # Get top 100 or fewer
    write.table(top_cis, file = top_cis_file, quote = FALSE, sep = "\t", 
                row.names = FALSE)
    cat("Wrote top cis interactions to:", top_cis_file, "\n")
  } else {
    cat("No cis interactions to write\n")
  }
  
  # Top trans interactions
  top_trans_file <- paste0(output_prefix, "_top_trans_interactions.txt")
  if (nrow(trans_interactions) > 0) {
    top_trans <- head(trans_interactions, 100)  # Get top 100 or fewer
    write.table(top_trans, file = top_trans_file, quote = FALSE, sep = "\t", 
                row.names = FALSE)
    cat("Wrote top trans interactions to:", top_trans_file, "\n")
  } else {
    cat("No trans interactions to write\n")
  }
  
  # Create a summary statistics file
  summary_file <- paste0(output_prefix, "_summary_statistics.txt")
  sink(summary_file)
  cat("DIFFERENTIAL INTERACTION ANALYSIS SUMMARY\n")
  cat("=========================================\n\n")
  
  if(include_empirical) {
    cat("Analysis includes empirical p-values from null dataset comparison\n\n")
  }
  
  cat("Cis Interactions Summary:\n")
  cat("-------------------------\n")
  cat("Total significant cis interactions:", nrow(cis_interactions), "\n")
  if (nrow(cis_interactions) > 0) {
    cat("Mean logFC:", mean(cis_interactions$logFC), "\n")
    cat("Mean FDR:", mean(cis_interactions$FDR), "\n")
    if(include_empirical && "empirical_FDR" %in% colnames(cis_interactions)) {
      cat("Mean empirical FDR:", mean(cis_interactions$empirical_FDR, na.rm = TRUE), "\n")
    }
    cat("Mean interaction distance:", mean(cis_interactions$interaction_distance, na.rm = TRUE), "bp\n")
    
    # Count by chromosome
    chr_counts <- table(cis_interactions$chr1)
    cat("\nDistribution by chromosome:\n")
    print(chr_counts)
    
    # Count by distance category
    if (any(!is.na(cis_interactions$interaction_distance))) {
      dist_breaks <- c(0, 10000, 50000, 100000, 500000, 1000000, Inf)
      dist_labels <- c("<10kb", "10-50kb", "50-100kb", "100-500kb", "500kb-1Mb", ">1Mb")
      distance_summary <- table(cut(cis_interactions$interaction_distance, 
                               breaks = dist_breaks, labels = dist_labels))
      cat("\nInteraction distance distribution:\n")
      print(distance_summary)
    }
  }
  
  cat("\n\nTrans Interactions Summary:\n")
  cat("---------------------------\n")
  cat("Total significant trans interactions:", nrow(trans_interactions), "\n")
  if (nrow(trans_interactions) > 0) {
    cat("Mean logFC:", mean(trans_interactions$logFC), "\n")
    cat("Mean FDR:", mean(trans_interactions$FDR), "\n")
    if(include_empirical && "empirical_FDR" %in% colnames(trans_interactions)) {
      cat("Mean empirical FDR:", mean(trans_interactions$empirical_FDR, na.rm = TRUE), "\n")
    }
    
    # Count by chromosome pair
    chr_pairs <- paste(trans_interactions$chr1, trans_interactions$chr2, sep = "-")
    chr_pair_counts <- sort(table(chr_pairs), decreasing = TRUE)
    cat("\nDistribution by chromosome pair:\n")
    print(head(chr_pair_counts, 10))  # Show top 10 chromosome pairs
  }
  
  sink()
  cat("Wrote summary statistics to:", summary_file, "\n")
  
  return(list(top_cis = top_cis_file, top_trans = top_trans_file, summary = summary_file))
}

# Function to plot comparison between dataset and null_dataset
plot_dataset_comparison <- function(dataset, null_dataset, output_prefix) {
  cat("Creating comparison plots...\n")
  
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("Warning: ggplot2 package not available. Skipping plots.\n")
    return(NULL)
  }
  
  library(ggplot2)
  
  # Create output directory for plots if it doesn't exist
  plots_dir <- paste0(output_prefix, "_plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  
  # Plot 1: logFC distribution comparison
  if("logFC" %in% colnames(dataset) && "logFC" %in% colnames(null_dataset)) {
    logfc_plot <- ggplot() +
      geom_density(data = dataset, aes(x = logFC, fill = "Dataset"), alpha = 0.5) +
      geom_density(data = null_dataset, aes(x = logFC, fill = "Null dataset"), alpha = 0.5) +
      scale_fill_manual(values = c("Dataset" = "blue", "Null dataset" = "red")) +
      labs(title = "logFC Distribution Comparison",
           x = "logFC",
           y = "Density",
           fill = "") +
      theme_minimal()
    
    # Save the plot
    logfc_plot_file <- file.path(plots_dir, "logFC_comparison.pdf")
    ggsave(logfc_plot_file, plot = logfc_plot, width = 8, height = 6)
    cat("Saved logFC comparison plot to:", logfc_plot_file, "\n")
  }
  
  # Plot 2: FDR distribution comparison
  if("FDR" %in% colnames(dataset) && "FDR" %in% colnames(null_dataset)) {
    fdr_plot <- ggplot() +
      geom_density(data = dataset, aes(x = -log10(FDR), fill = "Dataset"), alpha = 0.5) +
      geom_density(data = null_dataset, aes(x = -log10(FDR), fill = "Null dataset"), alpha = 0.5) +
      scale_fill_manual(values = c("Dataset" = "blue", "Null dataset" = "red")) +
      labs(title = "FDR Distribution Comparison",
           x = "-log10(FDR)",
           y = "Density",
           fill = "") +
      theme_minimal()
    
    # Save the plot
    fdr_plot_file <- file.path(plots_dir, "FDR_comparison.pdf")
    ggsave(fdr_plot_file, plot = fdr_plot, width = 8, height = 6)
    cat("Saved FDR comparison plot to:", fdr_plot_file, "\n")
  }
  
  # Plot 3: Scatter plot of logFC vs FDR for dataset
  if("logFC" %in% colnames(dataset) && "FDR" %in% colnames(dataset)) {
    volcano_plot <- ggplot(dataset, aes(x = logFC, y = -log10(FDR))) +
      geom_point(aes(color = abs(logFC) > LOG_FC_THRESHOLD & FDR < FDR_THRESHOLD)) +
      scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
      labs(title = "Volcano Plot of Dataset",
           x = "logFC",
           y = "-log10(FDR)",
           color = "Significant") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Save the plot
    volcano_plot_file <- file.path(plots_dir, "volcano_plot.pdf")
    ggsave(volcano_plot_file, plot = volcano_plot, width = 8, height = 8)
    cat("Saved volcano plot to:", volcano_plot_file, "\n")
  }
  
  # Plot 4: Interaction distance distribution comparison for cis interactions
  if("interaction_distance" %in% colnames(dataset) && "interaction_distance" %in% colnames(null_dataset)) {
    # Filter for cis interactions
    dataset_cis <- dataset[dataset$interaction_type == "cis", ]
    null_dataset_cis <- null_dataset[null_dataset$interaction_type == "cis", ]
    
    if(nrow(dataset_cis) > 0 && nrow(null_dataset_cis) > 0) {
      distance_plot <- ggplot() +
        geom_density(data = dataset_cis, aes(x = log10(interaction_distance), fill = "Dataset"), alpha = 0.5) +
        geom_density(data = null_dataset_cis, aes(x = log10(interaction_distance), fill = "Null dataset"), alpha = 0.5) +
        scale_fill_manual(values = c("Dataset" = "blue", "Null dataset" = "red")) +
        labs(title = "Interaction Distance Distribution (cis interactions)",
             x = "log10(Interaction Distance)",
             y = "Density",
             fill = "") +
        theme_minimal()
      
      # Save the plot
      distance_plot_file <- file.path(plots_dir, "distance_comparison.pdf")
      ggsave(distance_plot_file, plot = distance_plot, width = 8, height = 6)
      cat("Saved interaction distance comparison plot to:", distance_plot_file, "\n")
    }
  }
  
  cat("All comparison plots created in directory:", plots_dir, "\n")
  return(plots_dir)
}

# Main function
main <- function() {
  # Get arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript filter_interactions.R <interaction_file> <null_interaction_file> [output_prefix]\n")
    cat("Example: Rscript filter_interactions.R annotated_differential_interactions.csv null_model_interactions.csv my_analysis\n")
    return(1)
  }
  
  input_file <- args[1]
  null_file <- args[2]
  output_prefix <- if (length(args) >= 3) args[3] else "filtered_interactions"
  
  # Read input data
  interactions <- read_interaction_data(input_file)
  null_interactions <- read_interaction_data(null_file)
  
  # Compare dataset with null dataset
  interactions_with_empirical <- compare_with_null(interactions, null_interactions)
  
  # Create comparison plots
  plot_dataset_comparison(interactions_with_empirical, null_interactions, output_prefix)
  
  # Filter significant interactions
  sig_interactions <- filter_significant_interactions(interactions_with_empirical)
  
  # If no significant interactions, exit
  if (nrow(sig_interactions) == 0) {
    cat("No significant interactions found. Exiting.\n")
    return(0)
  }
  
  # Add interaction ID for later processing
  sig_interactions$interaction_id <- paste(
    sig_interactions$chr1, sig_interactions$start1, sig_interactions$end1,
    sig_interactions$chr2, sig_interactions$start2, sig_interactions$end2,
    sep = "_"
  )
  
  # Separate cis and trans interactions
  interaction_types <- separate_cis_trans(sig_interactions)
  cis_interactions <- interaction_types$cis
  trans_interactions <- interaction_types$trans
  
  # Sort by significance
  sorted_cis <- sort_by_significance(cis_interactions)
  sorted_trans <- sort_by_significance(trans_interactions)
  
  # Create BED files for downstream analysis
  if (nrow(sorted_cis) > 0) {
    redfly_bed <- prepare_redfly_input(sorted_cis, output_prefix)
    cdbest_bed_cis <- prepare_cdbest_input(sorted_cis, paste0(output_prefix, "_cis"))
  }
  
  if (nrow(sorted_trans) > 0) {
    droid_beds <- prepare_droid_input(sorted_trans, output_prefix)
    cdbest_bed_trans <- prepare_cdbest_input(sorted_trans, paste0(output_prefix, "_trans"))
  }
  
  # Create summary files
  include_empirical <- "empirical_FDR" %in% colnames(sig_interactions)
  summary_files <- create_summary_files(sorted_cis, sorted_trans, output_prefix, include_empirical)
  
  # Clean up temporary columns
  if("interaction_id" %in% colnames(sorted_cis)) {
    sorted_cis$interaction_id <- NULL
  }
  if("interaction_id" %in% colnames(sorted_trans)) {
    sorted_trans$interaction_id <- NULL
  }
  
  # Save the processed data
  write.table(sorted_cis, file = paste0(output_prefix, "_processed_cis.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(sorted_trans, file = paste0(output_prefix, "_processed_trans.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  cat("\nAnalysis complete!\n")
  return(0)
}

# Run the main function
main()