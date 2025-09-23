#!/usr/bin/env Rscript
# Combine diffHic Results from Multiple Resolutions and Infections
# 
# This script combines the results from individual resolution analyses
# and creates summary reports and visualizations.
#
# Usage: Rscript scripts/combine_diffhic_results.R --results_dir results/diffhic_results --resolutions 1000,8000,32000,128000

#=============================================================================
# 1. Setup Environment and Load Libraries
#=============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(VennDiagram)
  library(gridExtra)
  library(data.table)
  library(optparse)
  library(reshape2)
  library(dplyr)
})

# Define command-line options
option_list <- list(
  make_option(c("--results_dir"), type="character", default="results/diffhic_results",
              help="Directory containing the individual resolution results [default=%default]"),
  make_option(c("--resolutions"), type="character", default="1000,8000,32000,128000",
              help="Comma-separated list of resolutions that were analyzed [default=%default]"),
  make_option(c("--fdr"), type="double", default=0.01,
              help="FDR threshold for significance [default=%default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list, 
                          description="Combine multi-resolution diffHic results")
opt <- parse_args(opt_parser)

# Set parameters
results_dir <- opt$results_dir
resolutions <- as.numeric(unlist(strsplit(opt$resolutions, ",")))
fdr_threshold <- opt$fdr

cat("=================================================================\n")
cat("COMBINING DIFFHIC RESULTS ACROSS RESOLUTIONS AND INFECTIONS\n")
cat("=================================================================\n")
cat("Results directory:", results_dir, "\n")
cat("Resolutions:", paste(resolutions, "bp", collapse=", "), "\n")
cat("FDR threshold:", fdr_threshold, "\n\n")

# Create summary directory
summary_dir <- file.path(results_dir, "summary")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

#=============================================================================
# 2. Define Functions
#=============================================================================

# Function to extract infection metadata from filename
extract_infection_metadata <- function(filename) {
  pattern <- "annotated_differential_interactions_(.+)\\.csv$"
  match <- regmatches(filename, regexpr(pattern, filename, perl=TRUE))
  
  if (length(match) > 0) {
    infection <- gsub("annotated_differential_interactions_", "", match)
    infection <- gsub("\\.csv$", "", infection)
    return(infection)
  } else {
    return("unknown")
  }
}

# Function to create a unique identifier for interactions
create_interaction_id <- function(chr1, start1, end1, chr2, start2, end2) {
  # Sort chromosome pairs to ensure consistency
  is_swap <- as.character(chr2) < as.character(chr1) | 
            (as.character(chr2) == as.character(chr1) & start2 < start1)
  
  result <- character(length(chr1))
  
  for (i in 1:length(chr1)) {
    if (is_swap[i]) {
      id <- paste(chr2[i], start2[i], end2[i], chr1[i], start1[i], end1[i], sep = "_")
    } else {
      id <- paste(chr1[i], start1[i], end1[i], chr2[i], start2[i], end2[i], sep = "_")
    }
    result[i] <- id
  }
  
  return(result)
}

#=============================================================================
# 3. Load and Process All Results
#=============================================================================

all_results_combined <- NULL
summary_table <- data.frame()

cat("Loading results from each resolution...\n")

for (res in resolutions) {
  res_str <- as.character(res)
  res_dir <- file.path(results_dir, paste0("res_", res))
  
  cat("\nProcessing resolution", res, "bp...\n")
  
  # Check if the results directory exists
  if (!dir.exists(res_dir)) {
    cat("  WARNING: Directory not found:", res_dir, "\n")
    next
  }
  
  # Find ALL annotated results files for this resolution
  results_files <- list.files(res_dir, pattern="^annotated_differential_interactions_.*\\.csv$", full.names=TRUE)
  
  if (length(results_files) == 0) {
    cat("  WARNING: No annotated results files found in:", res_dir, "\n")
    next
  }
  
  cat("  Found", length(results_files), "infection condition(s)\n")
  
  # Process each infection condition file
  for (results_file in results_files) {
    filename <- basename(results_file)
    infection <- extract_infection_metadata(filename)
    
    cat("    Loading:", infection, "\n")
    
    # Read the results
    data <- read.csv(results_file, row.names=1)
    
    # Add metadata columns
    data$resolution <- as.integer(res)
    data$infection <- infection
    data$filename <- filename
    
    # Get significant results
    sig_data <- data[data$FDR < fdr_threshold, ]
    
    # Count interaction types
    all_cis <- sum(data$interaction_type == "cis", na.rm=TRUE)
    all_trans <- sum(data$interaction_type == "trans", na.rm=TRUE)
    sig_cis <- sum(sig_data$interaction_type == "cis", na.rm=TRUE)
    sig_trans <- sum(sig_data$interaction_type == "trans", na.rm=TRUE)
    
    # Add to summary table
    summary_row <- data.frame(
      Resolution = as.integer(res),
      Infection = infection,
      Total_Interactions = nrow(data),
      Total_Cis = all_cis,
      Total_Trans = all_trans,
      Significant_Interactions = nrow(sig_data),
      Significant_Cis = sig_cis,
      Significant_Trans = sig_trans,
      Percent_Sig_Cis = ifelse(nrow(sig_data) > 0, round(100 * sig_cis / nrow(sig_data), 1), 0),
      Percent_Sig_Trans = ifelse(nrow(sig_data) > 0, round(100 * sig_trans / nrow(sig_data), 1), 0),
      Min_Trans_FDR = ifelse(all_trans > 0, min(data[data$interaction_type == "trans", "FDR"], na.rm=TRUE), NA),
      Filename = filename,
      stringsAsFactors = FALSE
    )
    
    summary_table <- rbind(summary_table, summary_row)
    
    # Add to combined results
    all_results_combined <- rbind(all_results_combined, data)
    
    cat("      Total:", nrow(data), "interactions")
    cat(" | Significant:", nrow(sig_data))
    cat(" | Cis:", sig_cis, "| Trans:", sig_trans, "\n")
  }
}

# Sort summary table
summary_table <- summary_table[order(summary_table$Resolution, summary_table$Infection), ]

#=============================================================================
# 4. Save Results
#=============================================================================

cat("\n=================================================================\n")
cat("SAVING RESULTS\n")
cat("=================================================================\n")

# Save detailed summary
write.csv(summary_table, file.path(summary_dir, "detailed_summary.csv"), row.names=FALSE)
cat("Saved: detailed_summary.csv\n")

# Save all combined results
if (!is.null(all_results_combined) && nrow(all_results_combined) > 0) {
  write.csv(all_results_combined, 
            file.path(summary_dir, "all_results_combined.csv"), 
            row.names=FALSE)
  cat("Saved: all_results_combined.csv (", nrow(all_results_combined), "total interactions)\n")
  
  # Save only significant results
  all_sig_combined <- all_results_combined[all_results_combined$FDR < fdr_threshold, ]
  if (nrow(all_sig_combined) > 0) {
    write.csv(all_sig_combined, 
              file.path(summary_dir, "significant_results_combined.csv"), 
              row.names=FALSE)
    cat("Saved: significant_results_combined.csv (", nrow(all_sig_combined), "significant interactions)\n")
  }
}

#=============================================================================
# 5. Create Summary Report
#=============================================================================

sink(file.path(summary_dir, "analysis_summary.txt"))

cat("=====================================================\n")
cat("DIFFHIC MULTI-RESOLUTION ANALYSIS SUMMARY\n")
cat("=====================================================\n\n")

cat("Analysis Parameters:\n")
cat("- FDR threshold:", fdr_threshold, "\n")
cat("- Resolutions analyzed:", paste(resolutions, "bp", collapse=", "), "\n")
cat("- Infection conditions:", paste(unique(summary_table$Infection), collapse=", "), "\n\n")

cat("DETAILED RESULTS BY RESOLUTION AND INFECTION:\n")
cat("==============================================\n")
for (res in unique(summary_table$Resolution)) {
  cat("\nResolution", res, "bp:\n")
  res_data <- summary_table[summary_table$Resolution == res, ]
  
  for (i in 1:nrow(res_data)) {
    row <- res_data[i, ]
    cat(sprintf("  %s:\n", row$Infection))
    cat(sprintf("    Total interactions: %s (cis: %s, trans: %s)\n", 
               format(row$Total_Interactions, big.mark=","),
               format(row$Total_Cis, big.mark=","),
               format(row$Total_Trans, big.mark=",")))
    cat(sprintf("    Significant: %s (cis: %s [%.1f%%], trans: %s [%.1f%%])\n",
               row$Significant_Interactions, row$Significant_Cis, row$Percent_Sig_Cis,
               row$Significant_Trans, row$Percent_Sig_Trans))
    if (!is.na(row$Min_Trans_FDR)) {
      cat(sprintf("    Minimum trans FDR: %.4f\n", row$Min_Trans_FDR))
    } else {
      cat("    No trans interactions found\n")
    }
  }
}

cat("\nSUMMARY STATISTICS:\n")
cat("==================\n")
total_interactions <- sum(summary_table$Total_Interactions)
total_significant <- sum(summary_table$Significant_Interactions)
total_sig_cis <- sum(summary_table$Significant_Cis)
total_sig_trans <- sum(summary_table$Significant_Trans)

cat("Total interactions across all conditions:", format(total_interactions, big.mark=","), "\n")
cat("Total significant interactions:", format(total_significant, big.mark=","), "\n")
cat("Significant cis interactions:", format(total_sig_cis, big.mark=","), 
    sprintf(" (%.1f%% of significant)\n", 100 * total_sig_cis / total_significant))
cat("Significant trans interactions:", format(total_sig_trans, big.mark=","), 
    sprintf(" (%.1f%% of significant)\n", 100 * total_sig_trans / total_significant))

# Recommendations
cat("\nRECOMMENDations:\n")
cat("================\n")
if (total_sig_trans == 0) {
  min_trans_fdr <- min(summary_table$Min_Trans_FDR, na.rm=TRUE)
  cat("- No trans interactions are significant at FDR <", fdr_threshold, "\n")
  cat("- Minimum trans FDR observed:", sprintf("%.4f", min_trans_fdr), "\n")
  cat("- Consider using FDR threshold of", sprintf("%.3f", min_trans_fdr), "to include trans interactions\n")
  cat("- This suggests Wolbachia primarily affects local (cis) chromatin organization\n")
} else {
  best_condition <- summary_table[which.max(summary_table$Significant_Interactions), ]
  cat("- Best condition:", best_condition$Infection, "at", best_condition$Resolution, "bp resolution\n")
  cat("- Most significant interactions:", best_condition$Significant_Interactions, "\n")
}

# Best resolution analysis
res_totals <- aggregate(Significant_Interactions ~ Resolution, data=summary_table, FUN=sum)
best_res <- res_totals$Resolution[which.max(res_totals$Significant_Interactions)]
cat("- Resolution with most significant interactions:", best_res, "bp\n")

# Best infection analysis  
inf_totals <- aggregate(Significant_Interactions ~ Infection, data=summary_table, FUN=sum)
best_inf <- inf_totals$Infection[which.max(inf_totals$Significant_Interactions)]
cat("- Infection condition with most significant interactions:", best_inf, "\n")

sink()

#=============================================================================
# 6. Create Visualizations
#=============================================================================

if (nrow(summary_table) > 0) {
  cat("Creating visualizations...\n")
  
  # Create comprehensive plot
  pdf(file.path(summary_dir, "summary_plots.pdf"), width=14, height=10)
  
  # Plot 1: Significant interactions by resolution and infection
  plot_data <- summary_table %>%
    select(Resolution, Infection, Significant_Cis, Significant_Trans) %>%
    reshape2::melt(id.vars=c("Resolution", "Infection"), 
                   variable.name="Interaction_Type", value.name="Count")
  
  plot_data$Interaction_Type <- gsub("Significant_", "", plot_data$Interaction_Type)
  
  p1 <- ggplot(plot_data, aes(x=factor(Resolution), y=Count, fill=Interaction_Type)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=c("Cis"="darkblue", "Trans"="darkred")) +
    labs(x="Resolution (bp)", y="Number of Significant Interactions", 
         title="Significant Interactions by Resolution and Infection",
         fill="Type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~Infection, scales="free_y")
  
  print(p1)
  
  # Plot 2: Total vs Significant interactions
  p2 <- ggplot(summary_table, aes(x=factor(Resolution), y=Significant_Interactions, fill=Infection)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x="Resolution (bp)", y="Significant Interactions", 
         title="Significant Interactions by Resolution and Infection") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p2)
  
  # Plot 3: Heatmap if multiple conditions
  if (nrow(summary_table) > 1) {
    heatmap_data <- summary_table %>%
      select(Resolution, Infection, Significant_Interactions) %>%
      reshape2::dcast(Infection ~ Resolution, value.var="Significant_Interactions", fill=0)
    
    rownames(heatmap_data) <- heatmap_data$Infection
    heatmap_matrix <- as.matrix(heatmap_data[, -1])
    
    pheatmap(heatmap_matrix,
             main="Significant Interactions Heatmap",
             cluster_rows=FALSE, cluster_cols=FALSE,
             display_numbers=TRUE, number_color="white",
             color=colorRampPalette(c("white", "blue", "darkblue"))(50))
  }
  
  dev.off()
  cat("Saved: summary_plots.pdf\n")
}

cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=================================================================\n")
cat("Results saved in:", summary_dir, "\n")
cat("Key files:\n")
cat("- detailed_summary.csv: Complete summary table\n")
cat("- all_results_combined.csv: All interactions combined\n")
cat("- significant_results_combined.csv: Significant interactions only\n")
cat("- analysis_summary.txt: Detailed text report\n")
cat("- summary_plots.pdf: Visualization plots\n")

# Display summary table
cat("\nQUICK SUMMARY:\n")
print(summary_table[, c("Resolution", "Infection", "Significant_Interactions", "Significant_Cis", "Significant_Trans")])

cat("\nSession Information:\n")
sessionInfo()