#!/usr/bin/env Rscript
# Simplified Multi-Resolution diffHic Analysis with Volcano Plots

#=============================================================================
# 1. Setup and Load Libraries
#=============================================================================

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
  library(data.table)
  library(BiocParallel)
  library(optparse)
})

# Install ggrepel if needed
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

# Command-line options
option_list <- list(
  make_option(c("--data_dir"), type="character", default="interactions_data",
              help="Directory containing contact files"),
  make_option(c("--resolutions"), type="character", default="1000,8000,32000,128000",
              help="Comma-separated list of resolutions"),
  make_option(c("--min_count"), type="integer", default=2,
              help="Minimum read count"),
  make_option(c("--min_samples"), type="integer", default=2,
              help="Minimum samples with min_count"),
  make_option(c("--fdr"), type="double", default=0.05,
              help="FDR threshold"),
  make_option(c("--output_dir"), type="character", default="diffhic_results",
              help="Output directory"),
  make_option(c("--reference"), type="character", default="JW18DOX",
              help="Reference condition"),
  make_option(c("--threads"), type="integer", default=4,
              help="Number of threads")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Set parameters
data_dir <- opt$data_dir
resolutions <- as.numeric(unlist(strsplit(opt$resolutions, ",")))
min_count <- opt$min_count
min_samples <- opt$min_samples
fdr_threshold <- opt$fdr
output_dir <- opt$output_dir
reference_level <- opt$reference
threads <- opt$threads

register(MulticoreParam(threads))

cat("Parameters:\n")
cat("  Data dir:", data_dir, "\n")
cat("  Resolutions:", paste(resolutions, collapse=", "), "bp\n")
cat("  Min count:", min_count, "\n")
cat("  Min samples:", min_samples, "\n")
cat("  FDR:", fdr_threshold, "\n")
cat("  Reference:", reference_level, "\n\n")

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
for (res in resolutions) {
  dir.create(file.path(output_dir, paste0("res_", res)), showWarnings = FALSE)
}
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "comprehensive_tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "volcano_plots"), showWarnings = FALSE)

#=============================================================================
# 2. Helper Functions
#=============================================================================

# Read and filter data for a specific resolution
# Read and filter data for a specific resolution
read_resolution_data <- function(resolution, data_dir, min_count) {
  # Look for filtered CSV files instead of original TSV files
  contact_files <- list.files(data_dir, pattern = "_filtered_contacts\\.csv$", full.names = TRUE)
  cat("\nResolution", resolution, "bp - found", length(contact_files), "contact files\n")
  
  all_data <- list()
  
  for (file in contact_files) {
    cat("  Reading:", basename(file), "\n")
    
    # Read with data.table
    dt <- fread(file, header = TRUE)
    
    # Rename columns if they use 'chr' instead of 'chrom' (from filtered CSV files)
    if ("chr1" %in% colnames(dt)) {
      setnames(dt, old = c("chr1", "chr2"), new = c("chrom1", "chrom2"))
    }
    
    # Filter by resolution and count
    dt_filtered <- dt[dt$resolution == resolution & dt$count >= min_count]
    
    cat("    Filtered to", nrow(dt_filtered), "interactions (resolution=", resolution, ", count>=", min_count, ")\n")
    
    if (nrow(dt_filtered) == 0) {
      warning("    No data after filtering!")
      next
    }
    
    # Split by condition and replicate
    for (cond in unique(dt_filtered$condition)) {
      for (rep in unique(dt_filtered[dt_filtered$condition == cond, ]$replicate)) {
        key <- paste0(cond, "_rep", rep)
        subset_dt <- dt_filtered[dt_filtered$condition == cond & dt_filtered$replicate == rep, ]
        
        if (nrow(subset_dt) > 0) {
          all_data[[key]] <- as.data.frame(subset_dt)
          cat("    ", key, ":", nrow(subset_dt), "interactions\n")
        }
      }
    }
  }
  
  return(all_data)
}

# Create InteractionSet from data frame
create_iset <- function(data) {
  anchor1 <- GRanges(
    seqnames = data$chrom1,
    ranges = IRanges(start = data$start1, end = data$end1)
  )
  
  anchor2 <- GRanges(
    seqnames = data$chrom2,
    ranges = IRanges(start = data$start2, end = data$end2)
  )
  
  gi <- GInteractions(anchor1, anchor2)
  iset <- InteractionSet(list(counts = matrix(data$count, ncol = 1)), gi)
  
  return(iset)
}

# Combine all samples into unified InteractionSet
combine_isets <- function(iset_list) {
  cat("Combining", length(iset_list), "samples...\n")
  
  # Extract all unique interactions
  all_interactions <- data.frame()
  
  for (i in seq_along(iset_list)) {
    gi <- interactions(iset_list[[i]])
    a1 <- anchors(gi, type="first")
    a2 <- anchors(gi, type="second")
    
    df <- data.frame(
      chr1 = as.character(seqnames(a1)),
      start1 = start(a1),
      end1 = end(a1),
      chr2 = as.character(seqnames(a2)),
      start2 = start(a2),
      end2 = end(a2),
      count = assay(iset_list[[i]])[, 1],
      sample = i
    )
    
    all_interactions <- rbind(all_interactions, df)
  }
  
  # Create unique interaction IDs
  all_interactions$id <- paste(
    all_interactions$chr1, all_interactions$start1, all_interactions$end1,
    all_interactions$chr2, all_interactions$start2, all_interactions$end2,
    sep = "_"
  )
  
  # Get unique interactions
  unique_ints <- unique(all_interactions[, c("chr1", "start1", "end1", "chr2", "start2", "end2", "id")])
  cat("Found", nrow(unique_ints), "unique interactions\n")
  
  # Build count matrix
  count_matrix <- matrix(0, nrow = nrow(unique_ints), ncol = length(iset_list))
  colnames(count_matrix) <- names(iset_list)
  rownames(unique_ints) <- unique_ints$id
  
  for (i in 1:length(iset_list)) {
    sample_data <- all_interactions[all_interactions$sample == i, ]
    for (j in 1:nrow(sample_data)) {
      row_idx <- which(rownames(unique_ints) == sample_data$id[j])
      count_matrix[row_idx, i] <- count_matrix[row_idx, i] + sample_data$count[j]
    }
  }
  
  # Create GInteractions
  anchor1 <- GRanges(
    seqnames = unique_ints$chr1,
    ranges = IRanges(start = unique_ints$start1, end = unique_ints$end1)
  )
  
  anchor2 <- GRanges(
    seqnames = unique_ints$chr2,
    ranges = IRanges(start = unique_ints$start2, end = unique_ints$end2)
  )
  
  gi <- GInteractions(anchor1, anchor2)
  combined_iset <- InteractionSet(list(counts = count_matrix), interactions = gi)
  
  return(combined_iset)
}

# Add annotations to results
add_annotations <- function(results_table, iset) {
  gi <- interactions(iset)
  a1 <- anchors(gi, type="first")
  a2 <- anchors(gi, type="second")
  
  chr1 <- as.character(seqnames(a1))
  chr2 <- as.character(seqnames(a2))
  is_cis <- chr1 == chr2
  
  distance <- rep(NA, length(gi))
  distance[is_cis] <- end(a2)[is_cis] - start(a1)[is_cis]
  
  annotations <- data.frame(
    chr1 = chr1,
    start1 = start(a1),
    end1 = end(a1),
    chr2 = chr2,
    start2 = start(a2),
    end2 = end(a2),
    interaction_type = ifelse(is_cis, "cis", "trans"),
    interaction_distance = distance,
    row.names = rownames(results_table)
  )
  
  return(cbind(results_table, annotations))
}

# Generate comprehensive interaction tables
generate_interaction_tables <- function(results_list, output_dir, fdr_threshold) {
  cat("\nGenerating comprehensive interaction tables...\n")
  
  tables_dir <- file.path(output_dir, "comprehensive_tables")
  dir.create(tables_dir, showWarnings = FALSE)
  
  all_sig_interactions <- NULL
  
  for (res in names(results_list)) {
    if (!is.null(results_list[[res]]) && 
        !is.null(results_list[[res]]$significant) && 
        nrow(results_list[[res]]$significant) > 0) {
      
      sig_results <- results_list[[res]]$significant
      sig_results$resolution <- as.numeric(res)
      all_sig_interactions <- rbind(all_sig_interactions, sig_results)
      
      write.csv(sig_results, 
                file.path(tables_dir, paste0("sig_interactions_", res, "bp.csv")), 
                row.names = FALSE)
      
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
  
  if (!is.null(all_sig_interactions) && nrow(all_sig_interactions) > 0) {
    write.csv(all_sig_interactions, 
              file.path(tables_dir, "all_sig_interactions_combined.csv"), 
              row.names = FALSE)
    
    all_cis <- all_sig_interactions[all_sig_interactions$interaction_type == "cis", ]
    all_trans <- all_sig_interactions[all_sig_interactions$interaction_type == "trans", ]
    
    if (nrow(all_cis) > 0) {
      write.csv(all_cis, file.path(tables_dir, "all_sig_cis_interactions.csv"), row.names = FALSE)
    }
    
    if (nrow(all_trans) > 0) {
      write.csv(all_trans, file.path(tables_dir, "all_sig_trans_interactions.csv"), row.names = FALSE)
    }
  }
  
  cat("Tables saved to:", tables_dir, "\n")
}

# Create custom volcano plot
create_custom_volcano_plot <- function(data, title, output_file, interaction_type = NULL, 
                                      logFC_col = "logFC", pval_col = "PValue", fdr_col = "FDR", 
                                      fdr_threshold = 0.05, logFC_threshold = 1) {
  
  # Filter by interaction type if specified
  if (!is.null(interaction_type)) {
    plot_data <- data[data$interaction_type == interaction_type, ]
    subtitle <- paste0(" (", interaction_type, " interactions)")
  } else {
    plot_data <- data
    subtitle <- " (all interactions)"
  }
  
  if (nrow(plot_data) == 0) {
    warning(paste0("No data to plot for ", title, " ", interaction_type, " interactions"))
    return(NULL)
  }
  
  # Create labels
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
  
  # Color coding
  plot_data$color <- "grey80"
  plot_data$color[abs(plot_data[[logFC_col]]) > logFC_threshold] <- "grey60"
  plot_data$color[plot_data[[fdr_col]] < fdr_threshold] <- "grey40"
  
  sig_idx <- plot_data[[fdr_col]] < fdr_threshold & abs(plot_data[[logFC_col]]) > logFC_threshold
  plot_data$color[sig_idx & plot_data[[logFC_col]] > 0] <- "#8ecc85"  # UP in reference
  plot_data$color[sig_idx & plot_data[[logFC_col]] < 0] <- "#1bab4b"  # UP in treatment
  
  # Label top 15 points
  plot_data$to_label <- FALSE
  sig_points <- plot_data[sig_idx, ]
  
  if (nrow(sig_points) > 0) {
    sig_points <- sig_points[order(sig_points[[fdr_col]], -abs(sig_points[[logFC_col]])), ]
    top_n <- min(15, nrow(sig_points))
    top_sig_idx <- rownames(sig_points)[1:top_n]
    plot_data$to_label[rownames(plot_data) %in% top_sig_idx] <- TRUE
  }
  
  # Count significant points
  sig_up <- sum(sig_idx & plot_data[[logFC_col]] > 0)
  sig_down <- sum(sig_idx & plot_data[[logFC_col]] < 0)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = .data[[logFC_col]], y = neg_log10_pval)) +
    geom_point(aes(color = color), size = 2, alpha = 0.75) +
    scale_color_identity() +
    geom_text_repel(
      data = subset(plot_data, to_label),
      aes(label = label),
      size = 2.5,
      box.padding = 0.5,
      point.padding = 0.2,
      segment.color = "grey50",
      max.overlaps = 15,
      min.segment.length = 0
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "grey") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(title, subtitle),
      subtitle = paste0("FDR < ", fdr_threshold, ", |log2FC| > ", logFC_threshold),
      x = "log2 Fold Change",
      y = "-log10(P-value)",
      caption = paste0('Total = ', nrow(plot_data), ' | Labeled top ', sum(plot_data$to_label))
    ) +
    annotate(
      "text",
      x = max(plot_data[[logFC_col]], na.rm = TRUE) * 0.85,
      y = max(plot_data$neg_log10_pval, na.rm = TRUE) * 0.9,
      label = paste0("UP (positive): ", sig_up, "\nUP (negative): ", sig_down),
      hjust = 1,
      size = 3.5
    )
  
  pdf(output_file, width = 12, height = 10)
  print(p)
  dev.off()
  
  cat("Created volcano plot:", output_file, "\n")
  
  return(p)
}

# Create volcano plots for all resolutions
create_custom_volcano_plots <- function(results_list, output_dir, fdr_threshold) {
  cat("\nCreating custom volcano plots...\n")
  
  volcano_dir <- file.path(output_dir, "volcano_plots")
  dir.create(volcano_dir, showWarnings = FALSE)
  
  all_annotated_results <- NULL
  
  for (res in names(results_list)) {
    if (!is.null(results_list[[res]]) && 
        !is.null(results_list[[res]]$all_interactions) && 
        nrow(results_list[[res]]$all_interactions) > 0) {
      
      annotated_results <- results_list[[res]]$all_interactions
      annotated_results$resolution <- as.numeric(res)
      all_annotated_results <- rbind(all_annotated_results, annotated_results)
      
      create_custom_volcano_plot(
        data = annotated_results,
        title = paste0("Differential Interactions at ", res, "bp Resolution"),
        output_file = file.path(volcano_dir, paste0("volcano_", res, "bp_all.pdf")),
        fdr_threshold = fdr_threshold
      )
      
      create_custom_volcano_plot(
        data = annotated_results,
        title = paste0("Differential Interactions at ", res, "bp Resolution"),
        output_file = file.path(volcano_dir, paste0("volcano_", res, "bp_cis.pdf")),
        interaction_type = "cis",
        fdr_threshold = fdr_threshold
      )
      
      create_custom_volcano_plot(
        data = annotated_results,
        title = paste0("Differential Interactions at ", res, "bp Resolution"),
        output_file = file.path(volcano_dir, paste0("volcano_", res, "bp_trans.pdf")),
        interaction_type = "trans",
        fdr_threshold = fdr_threshold
      )
    }
  }
  
  if (!is.null(all_annotated_results) && nrow(all_annotated_results) > 0) {
    create_custom_volcano_plot(
      data = all_annotated_results,
      title = "Combined Differential Interactions Across All Resolutions",
      output_file = file.path(volcano_dir, "volcano_combined_all.pdf"),
      fdr_threshold = fdr_threshold
    )
    
    create_custom_volcano_plot(
      data = all_annotated_results,
      title = "Combined Differential Interactions",
      output_file = file.path(volcano_dir, "volcano_combined_cis.pdf"),
      interaction_type = "cis",
      fdr_threshold = fdr_threshold
    )
    
    create_custom_volcano_plot(
      data = all_annotated_results,
      title = "Combined Differential Interactions",
      output_file = file.path(volcano_dir, "volcano_combined_trans.pdf"),
      interaction_type = "trans",
      fdr_threshold = fdr_threshold
    )
    
    create_custom_volcano_plot(
      data = all_annotated_results,
      title = "Top Differential Interactions",
      output_file = file.path(volcano_dir, "volcano_combined_top_hits.pdf"),
      fdr_threshold = fdr_threshold / 10,
      logFC_threshold = 1.5
    )
  }
  
  cat("Volcano plots saved to:", volcano_dir, "\n")
}

#=============================================================================
# 3. Process Each Resolution
#=============================================================================

process_resolution <- function(resolution) {
  output_res_dir <- file.path(output_dir, paste0("res_", resolution))
  
  cat("\n========================================\n")
  cat("Processing resolution:", resolution, "bp\n")
  cat("========================================\n")
  
  # Read data
  sample_data <- read_resolution_data(resolution, data_dir, min_count)
  
  if (length(sample_data) == 0) {
    warning("No data for resolution ", resolution)
    return(NULL)
  }
  
  # Create InteractionSets
  iset_list <- lapply(sample_data, create_iset)
  
  # Combine all samples
  combined_iset <- combine_isets(iset_list)
  combined_iset$totals <- colSums(assay(combined_iset))
  
  cat("Combined set:", nrow(combined_iset), "interactions x", ncol(assay(combined_iset)), "samples\n")
  
  # Filter low-abundance interactions
  keep <- rowSums(assay(combined_iset) >= min_count) >= min_samples
  filtered_iset <- combined_iset[keep, ]
  cat("After filtering:", nrow(filtered_iset), "interactions\n")
  
  # Create DGEList
  y <- asDGEList(filtered_iset)
  
  # Setup design
  samples <- colnames(y)
  conditions <- sapply(strsplit(samples, "_rep"), `[`, 1)
  replicates <- sapply(strsplit(samples, "_rep"), `[`, 2)
  
  # Set factor levels with reference first
  all_conditions <- unique(conditions)
  if (reference_level %in% all_conditions) {
    condition_levels <- c(reference_level, setdiff(all_conditions, reference_level))
  } else {
    condition_levels <- sort(all_conditions)
  }
  
  infection <- factor(conditions, levels = condition_levels)
  
  sample_table <- data.frame(
    sample = samples,
    infection = infection,
    replicate = factor(replicates),
    row.names = samples
  )
  
  print(sample_table)
  
  design <- model.matrix(~infection, data = sample_table)
  
  # Normalize and estimate dispersion
  y <- calcNormFactors(y)
  logcounts <- cpm(y, log = TRUE)
  y <- estimateDisp(y, design)
  
  # QC plots
  pdf(file.path(output_res_dir, "qc_plots.pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))
  
  # PCA
  pca <- prcomp(t(logcounts))
  plot(pca$x[, 1:2], col = as.numeric(infection), pch = 19,
       xlab = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       ylab = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"),
       main = paste0("PCA (", resolution, "bp)"))
  text(pca$x[, 1:2], labels = samples, pos = 3, cex = 0.7)
  legend("topright", legend = levels(infection), col = 1:length(levels(infection)), pch = 19)
  
  # BCV
  plotBCV(y, main = paste0("BCV (", resolution, "bp)"))
  
  dev.off()
  
  # Fit models
  fit <- glmQLFit(y, design)
  
  # Test each coefficient
  coef_names <- colnames(design)[-1]
  results_list <- list()
  
  for (i in seq_along(coef_names)) {
    coef_name <- coef_names[i]
    cat("\nTesting:", coef_name, "\n")
    
    qlf <- glmQLFTest(fit, coef = i + 1)
    top <- topTags(qlf, n = Inf)
    
    # Add annotations
    annotated <- add_annotations(top$table, filtered_iset)
    
    # Count significant
    sig <- annotated[annotated$FDR < fdr_threshold, ]
    cat("  Significant:", nrow(sig), "interactions (FDR <", fdr_threshold, ")\n")
    
    # Save results
    write.csv(top$table, 
              file.path(output_res_dir, paste0("differential_interactions_", coef_name, ".csv")))
    write.csv(annotated, 
              file.path(output_res_dir, paste0("annotated_differential_interactions_", coef_name, ".csv")))
    
    if (nrow(sig) > 0) {
      write.csv(sig, 
                file.path(output_res_dir, paste0("significant_interactions_", coef_name, ".csv")))
    }
    
    # Plots
    pdf(file.path(output_res_dir, paste0("plots_", coef_name, ".pdf")), width = 10, height = 5)
    par(mfrow = c(1, 2))
    
    plotMD(qlf, main = paste0("MA plot: ", coef_name))
    abline(h = c(-1, 0, 1), col = c("blue", "red", "blue"), lty = c(2, 1, 2))
    
    plot(top$table$logFC, -log10(top$table$PValue),
         pch = 20, col = ifelse(top$table$FDR < fdr_threshold, "red", "black"),
         xlab = "log2 FC", ylab = "-log10(P)",
         main = paste0("Volcano: ", coef_name))
    abline(v = c(-1, 1), h = -log10(0.05), lty = 2, col = "blue")
    
    dev.off()
    
    results_list[[coef_name]] <- list(
      top = top,
      annotated = annotated,
      significant = sig
    )
  }
  
  # Null model
  null_fit <- glmFit(y, model.matrix(~1, data = sample_table))
  null_results <- data.frame(
    logFC = rep(0, nrow(y)),
    logCPM = aveLogCPM(y),
    PValue = rep(1, nrow(y)),
    FDR = rep(1, nrow(y)),
    row.names = rownames(y)
  )
  write.csv(null_results, file.path(output_res_dir, "null_model_results.csv"))
  
  # Model comparison
  lr_stat <- sum(null_fit$deviance) - sum(fit$deviance)
  lr_df <- ncol(design) - 1
  lr_pval <- 1 - pchisq(lr_stat, lr_df)
  
  model_comparison <- data.frame(
    Model = c("Full (~infection)", "Null (intercept)"),
    LogLikelihood = c(sum(fit$loglik), sum(null_fit$loglik)),
    Deviance = c(sum(fit$deviance), sum(null_fit$deviance)),
    LRT_stat = c(lr_stat, NA),
    LRT_df = c(lr_df, NA),
    LRT_pval = c(lr_pval, NA)
  )
  write.csv(model_comparison, file.path(output_res_dir, "model_comparison.csv"), row.names = FALSE)
  
  # Return first coefficient results for summary
  first_coef <- results_list[[1]]
  
  return(list(
    resolution = resolution,
    significant = first_coef$significant,
    all_interactions = first_coef$annotated,
    filtered_iset = filtered_iset,
    model_comparison = model_comparison,
    null_results = null_results
  ))
}

#=============================================================================
# 4. Run Analysis
#=============================================================================

results_list <- list()
for (res in resolutions) {
  results_list[[as.character(res)]] <- process_resolution(res)
}

#=============================================================================
# 5. Generate Summary and Plots
#=============================================================================

cat("\n========================================\n")
cat("Generating summary and plots\n")
cat("========================================\n")

# Generate comprehensive tables
generate_interaction_tables(results_list, output_dir, fdr_threshold)

# Generate volcano plots
create_custom_volcano_plots(results_list, output_dir, fdr_threshold)

# Summary table
summary_table <- data.frame(
  Resolution = resolutions,
  Total = sapply(results_list, function(x) if (!is.null(x$filtered_iset)) nrow(x$filtered_iset) else 0),
  Significant = sapply(results_list, function(x) if (!is.null(x$significant)) nrow(x$significant) else 0),
  Cis = sapply(results_list, function(x) if (!is.null(x$significant)) sum(x$significant$interaction_type == "cis") else 0),
  Trans = sapply(results_list, function(x) if (!is.null(x$significant)) sum(x$significant$interaction_type == "trans") else 0)
)

write.csv(summary_table, file.path(summary_dir, "resolution_summary.csv"), row.names = FALSE)

# Text summary
sink(file.path(summary_dir, "comprehensive_summary.txt"))
cat("MULTI-RESOLUTION ANALYSIS SUMMARY\n")
cat("=================================\n\n")
cat("Parameters:\n")
cat("  Min count:", min_count, "\n")
cat("  Min samples:", min_samples, "\n")
cat("  FDR threshold:", fdr_threshold, "\n\n")
cat("Results by Resolution:\n")
print(summary_table)
cat("\n")

for (res in names(results_list)) {
  if (!is.null(results_list[[res]]$model_comparison)) {
    mc <- results_list[[res]]$model_comparison
    cat(res, "bp resolution:\n")
    cat("  LRT p-value:", mc$LRT_pval[1], "\n")
    cat("  Significant:", ifelse(mc$LRT_pval[1] < 0.05, "YES", "NO"), "\n\n")
  }
}
sink()

cat("\nAnalysis complete! Results in:", output_dir, "\n")
cat("Session info:\n")
print(sessionInfo())