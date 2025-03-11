#!/usr/bin/env Rscript

# Script for analyzing PICRUSt2 functional predictions
# This script processes the output from PICRUSt2 and generates visualizations and statistical analyses

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(RColorBrewer)
  library(DESeq2)
  library(pheatmap)
  library(knitr)
  library(vegan)
  library(randomForest)
  library(gridExtra)
  library(KEGGREST)
  library(patchwork)
  library(ggdendro)
  library(ape)
  library(ggridges)
})

# Parse command line arguments
option_list <- list(
  make_option("--ko", type="character", help="Path to KEGG Orthology (KO) predictions"),
  make_option("--ec", type="character", help="Path to Enzyme Commission (EC) predictions"),
  make_option("--pathways", type="character", help="Path to MetaCyc pathway predictions"),
  make_option("--metadata", type="character", help="Path to metadata file"),
  make_option("--group-col", type="character", help="Column name for group comparisons"),
  make_option("--output-dir", type="character", help="Output directory for results")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# Function to read and process PICRUSt2 output files
read_picrust2 <- function(file_path) {
  if (grepl("\\.qza$", file_path)) {
    # For QIIME2 artifacts, we need to extract the file first
    temp_dir <- tempdir()
    system(paste("qiime tools export --input-path", file_path, "--output-path", temp_dir))
    
    # Find the exported biom file
    biom_file <- list.files(temp_dir, pattern="\\.biom$", full.names=TRUE)[1]
    
    if (is.na(biom_file)) {
      stop("Could not find exported biom file from QIIME2 artifact")
    }
    
    # Convert biom to TSV
    system(paste("biom convert -i", biom_file, "-o", paste0(temp_dir, "/table.tsv"), "--to-tsv"))
    
    # Read the TSV file
    data <- read.table(paste0(temp_dir, "/table.tsv"), 
                       sep="\t", header=TRUE, row.names=1, comment.char="", skip=1)
    
  } else if (grepl("\\.tsv$", file_path)) {
    # Direct TSV file from PICRUSt2
    data <- read.table(file_path, sep="\t", header=TRUE, row.names=1, comment.char="")
  } else {
    stop("Unsupported file format. Please provide a .qza or .tsv file.")
  }
  
  # Transpose to get samples as rows
  data_t <- t(data)
  
  return(data_t)
}

# Read metadata
cat("Reading metadata...\n")
metadata <- read.table(opt$metadata, sep="\t", header=TRUE, row.names=1, comment.char="")

# Read PICRUSt2 outputs
cat("Reading PICRUSt2 outputs...\n")
ko_data <- NULL
ec_data <- NULL
pathway_data <- NULL

if (!is.null(opt$ko)) {
  ko_data <- read_picrust2(opt$ko)
  cat("KO predictions loaded with", nrow(ko_data), "samples and", ncol(ko_data), "KOs\n")
}

if (!is.null(opt$ec)) {
  ec_data <- read_picrust2(opt$ec)
  cat("EC predictions loaded with", nrow(ec_data), "samples and", ncol(ec_data), "ECs\n")
}

if (!is.null(opt$pathways)) {
  pathway_data <- read_picrust2(opt$pathways)
  cat("Pathway predictions loaded with", nrow(pathway_data), "samples and", ncol(pathway_data), "pathways\n")
}

# Check that we have at least one dataset
if (is.null(ko_data) && is.null(ec_data) && is.null(pathway_data)) {
  stop("No input data provided. Please provide at least one of KO, EC, or pathway predictions.")
}

# Check that group column exists in metadata
if (!(opt$group_col %in% colnames(metadata))) {
  stop(paste0("Group column '", opt$group_col, "' not found in metadata"))
}

# Ensure metadata row names match sample IDs
fix_metadata_for_data <- function(meta, data) {
  # Get common samples
  common_samples <- intersect(rownames(meta), rownames(data))
  
  if (length(common_samples) == 0) {
    # Try to match ignoring case
    rownames_meta_lower <- tolower(rownames(meta))
    rownames_data_lower <- tolower(rownames(data))
    
    # Create mapping
    meta_to_data <- match(rownames_data_lower, rownames_meta_lower)
    valid_matches <- !is.na(meta_to_data)
    
    if (sum(valid_matches) == 0) {
      stop("No matching samples between metadata and data")
    }
    
    # Create subset of metadata with matching sample IDs
    meta_subset <- meta[meta_to_data[valid_matches], , drop=FALSE]
    rownames(meta_subset) <- rownames(data)[valid_matches]
    
    cat("Found", sum(valid_matches), "samples after case-insensitive matching\n")
    return(list(meta=meta_subset, data=data[valid_matches, , drop=FALSE]))
  } else {
    cat("Found", length(common_samples), "samples in common between metadata and data\n")
    return(list(meta=meta[common_samples, , drop=FALSE], data=data[common_samples, , drop=FALSE]))
  }
}

# Create PDF for all plots
pdf(file.path(opt$output_dir, "picrust2_analysis.pdf"), width=12, height=10)

# Process and analyze KEGG Orthology (KO) data
if (!is.null(ko_data)) {
  cat("Analyzing KEGG Orthology (KO) predictions...\n")
  
  # Match metadata to KO data
  matched_ko <- fix_metadata_for_data(metadata, ko_data)
  ko_meta <- matched_ko$meta
  ko_data <- matched_ko$data
  
  # Remove features with zero variance
  ko_var <- apply(ko_data, 2, var)
  ko_data <- ko_data[, ko_var > 0, drop=FALSE]
  
  # Filter low-abundance KOs
  ko_means <- colMeans(ko_data)
  ko_data <- ko_data[, ko_means > quantile(ko_means, 0.25), drop=FALSE]
  
  # Normalize data (relative abundance)
  ko_rel <- t(apply(ko_data, 1, function(x) x / sum(x)))
  
  # Get group information
  groups <- as.factor(ko_meta[[opt$group_col]])
  
  # Principal component analysis
  ko_pca <- prcomp(ko_rel, scale=TRUE)
  
  # Extract PC scores
  pc_data <- as.data.frame(ko_pca$x)
  pc_data$Group <- groups
  
  # Calculate percent variance explained
  pc_var <- (ko_pca$sdev^2) / sum(ko_pca$sdev^2) * 100
  
  # PCA plot
  p_pca <- ggplot(pc_data, aes(x=PC1, y=PC2, color=Group)) +
    geom_point(size=3) +
    stat_ellipse(aes(group=Group), type="t", linetype=2) +
    theme_bw() +
    labs(title="PCA of KEGG Orthology (KO) Functional Profiles",
         x=paste0("PC1 (", round(pc_var[1], 1), "%)"),
         y=paste0("PC2 (", round(pc_var[2], 1), "%)")) +
    theme(legend.position="bottom")
  
  print(p_pca)
  
  # PERMANOVA to test for group differences
  ko_dist <- vegdist(ko_rel, method="bray")
  
  adonis_res <- adonis2(ko_dist ~ groups)
  
  # Display PERMANOVA results
  grid.newpage()
  grid.text(paste("PERMANOVA Results for KO Profiles\n\n",
                 "R-squared: ", round(adonis_res$R2[1], 3), "\n",
                 "p-value: ", round(adonis_res$`Pr(>F)`[1], 3)),
            x=0.5, y=0.5, gp=gpar(cex=1.2))
  
  # Differential abundance analysis using DESeq2
  # Round counts for DESeq2 (it requires integers)
  ko_counts <- round(ko_data)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = ko_counts,
    colData = data.frame(group=groups, row.names=rownames(ko_counts)),
    design = ~group
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get number of groups
  group_levels <- levels(groups)
  
  if(length(group_levels) == 2) {
    # Simple pairwise comparison for 2 groups
    res <- results(dds)
    
    # Identify significant features
    sig_features <- which(res$padj < 0.05 & !is.na(res$padj))
    
    if(length(sig_features) > 0) {
      # Sort by significance
      sig_res <- res[sig_features,]
      sig_res <- sig_res[order(sig_res$padj),]
      
      # Create volcano plot
      res_df <- as.data.frame(res)
      res_df$feature <- rownames(res_df)
      res_df$significant <- ifelse(res_df$padj < 0.05 & !is.na(res_df$padj), "Yes", "No")
      
      p_volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
        geom_point(alpha=0.7) +
        scale_color_manual(values=c("No"="gray", "Yes"="red")) +
        geom_hline(yintercept=-log10(0.05), linetype="dashed") +
        geom_vline(xintercept=c(-1, 1), linetype="dashed") +
        theme_bw() +
        labs(title="Differential Abundance Analysis of KO Features",
             x="Log2 Fold Change",
             y="-Log10 Adjusted P-value")
      
      print(p_volcano)
      
      # Get KEGG annotations for significant features
      sig_ko_ids <- rownames(sig_res)
      
      # Try to get KO annotations from KEGGREST
      tryCatch({
        ko_list <- keggList("ko")
        ko_info <- ko_list[sig_ko_ids]
        
        # Create data frame with KO information
        ko_info_df <- data.frame(
          ko_id = sig_ko_ids,
          description = ko_info[sig_ko_ids],
          log2FoldChange = sig_res$log2FoldChange,
          padj = sig_res$padj
        )
        
        # Save to file
        write.csv(ko_info_df, file=file.path(opt$output_dir, "significant_ko_features.csv"), row.names=FALSE)
        
      }, error=function(e) {
        warning("Could not retrieve KO annotations from KEGGREST. Continuing without annotations.")
        
        # Save simple table without annotations
        sig_res_df <- as.data.frame(sig_res)
        sig_res_df$ko_id <- rownames(sig_res_df)
        write.csv(sig_res_df, file=file.path(opt$output_dir, "significant_ko_features.csv"), row.names=FALSE)
      })
      
      # Create heatmap of significant features
      sig_ko_data <- ko_rel[, sig_ko_ids]
      
      # Scale data for better visualization
      sig_ko_data_scaled <- t(scale(t(sig_ko_data)))
      
      # Create annotation data frame
      annotation_df <- data.frame(Group=groups, row.names=rownames(sig_ko_data))
      
      # Create color palette for groups
      group_colors <- brewer.pal(max(3, length(unique(groups))), "Set1")[1:length(unique(groups))]
      names(group_colors) <- unique(groups)
      annotation_colors <- list(Group=group_colors)
      
      # Create heatmap
      pheatmap(sig_ko_data_scaled, annotation_row=annotation_df, annotation_colors=annotation_colors,
               fontsize_row=8, fontsize_col=8, main="Significant KO Features",
               show_colnames=FALSE, scale="none")
    } else {
      grid.newpage()
      grid.text("No significant differentially abundant KO features found", x=0.5, y=0.5, gp=gpar(cex=1.2))
    }
  } else {
    # Multiple group comparisons
    # Perform pairwise comparisons
    result_list <- list()
    
    for(i in 1:(length(group_levels)-1)) {
      for(j in (i+1):length(group_levels)) {
        group1 <- group_levels[i]
        group2 <- group_levels[j]
        
        # Subset data for these two groups
        idx <- groups %in% c(group1, group2)
        if(sum(idx) < 3) next  # Skip if fewer than 3 samples
        
        sub_dds <- DESeqDataSetFromMatrix(
          countData = ko_counts[idx, ],
          colData = data.frame(group=groups[idx], row.names=rownames(ko_counts)[idx]),
          design = ~group
        )
        
        # Run DESeq2
        sub_dds <- DESeq(sub_dds)
        
        # Get results
        res <- results(sub_dds, contrast=c("group", group2, group1))
        
        # Find significant features
        sig_idx <- which(res$padj < 0.05 & !is.na(res$padj))
        
        if(length(sig_idx) > 0) {
          sig_res <- res[sig_idx,]
          sig_res <- as.data.frame(sig_res[order(sig_res$padj),])
          sig_res$comparison <- paste(group2, "vs", group1)
          sig_res$feature <- rownames(sig_res)
          
          result_list[[length(result_list) + 1]] <- sig_res
        }
      }
    }
    
    if(length(result_list) > 0) {
      all_results <- do.call(rbind, result_list)
      
      # Save results
      write.csv(all_results, file=file.path(opt$output_dir, "ko_pairwise_comparisons.csv"), row.names=FALSE)
      
      # Create summary plot
      p_summary <- ggplot(all_results, aes(x=comparison, y=feature, fill=log2FoldChange)) +
        geom_tile() +
        scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
        theme_bw() +
        theme(axis.text.y=element_text(size=8),
              axis.text.x=element_text(angle=45, hjust=1)) +
        labs(title="Significant KO Features by Comparison",
             x="Comparison", y="Feature", fill="Log2 Fold Change")
      
      print(p_summary)
    } else {
      grid.newpage()
      grid.text("No significant differentially abundant KO features found in any pairwise comparison", 
                x=0.5, y=0.5, gp=gpar(cex=1.2))
    }
  }
  
  # Random Forest classification
  if(length(unique(groups)) > 1 && length(unique(groups)) < nrow(ko_rel)) {
    cat("Performing Random Forest classification...\n")
    
    set.seed(123)
    rf_model <- randomForest(x=ko_rel, y=groups, importance=TRUE, ntree=500)
    
    # Print model accuracy
    grid.newpage()
    grid.text(paste("Random Forest Classification of Groups Using KO Features\n\n",
                    "Out-of-bag error rate: ", round(rf_model$err.rate[nrow(rf_model$err.rate), "OOB"] * 100, 1), "%\n\n",
                    "Confusion Matrix:"),
              x=0.5, y=0.8, gp=gpar(cex=1.2))
    
    grid.table(rf_model$confusion, vp=viewport(x=0.5, y=0.4, width=0.8, height=0.4))
    
    # Extract feature importance
    importance_df <- as.data.frame(importance(rf_model))
    importance_df$feature <- rownames(importance_df)
    
    # Sort by importance
    importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing=TRUE), ]
    
    # Plot top 20 features
    top_features <- head(importance_df, 20)
    top_features$feature <- factor(top_features$feature, levels=top_features$feature[order(top_features$MeanDecreaseGini)])
    
    p_importance <- ggplot(top_features, aes(x=MeanDecreaseGini, y=feature)) +
      geom_bar(stat="identity", fill="steelblue") +
      theme_bw() +
      labs(title="Top 20 KO Features by Importance in Random Forest Classification",
           x="Mean Decrease in Gini Index", y="Feature")
    
    print(p_importance)
    
    # Save feature importance
    write.csv(importance_df, file=file.path(opt$output_dir, "ko_feature_importance.csv"), row.names=FALSE)
  }
}

# Process and analyze Enzyme Commission (EC) data
if (!is.null(ec_data)) {
  cat("Analyzing Enzyme Commission (EC) predictions...\n")
  
  # Match metadata to EC data
  matched_ec <- fix_metadata_for_data(metadata, ec_data)
  ec_meta <- matched_ec$meta
  ec_data <- matched_ec$data
  
  # Remove features with zero variance
  ec_var <- apply(ec_data, 2, var)
  ec_data <- ec_data[, ec_var > 0, drop=FALSE]
  
  # Filter low-abundance ECs
  ec_means <- colMeans(ec_data)
  ec_data <- ec_data[, ec_means > quantile(ec_means, 0.25), drop=FALSE]
  
  # Normalize data (relative abundance)
  ec_rel <- t(apply(ec_data, 1, function(x) x / sum(x)))
  
  # Get group information
  groups <- as.factor(ec_meta[[opt$group_col]])
  
  # Perform PCA
  ec_pca <- prcomp(ec_rel, scale=TRUE)
  
  # Extract PC scores
  pc_data <- as.data.frame(ec_pca$x)
  pc_data$Group <- groups
  
  # Calculate percent variance explained
  pc_var <- (ec_pca$sdev^2) / sum(ec_pca$sdev^2) * 100
  
  # PCA plot
  p_pca <- ggplot(pc_data, aes(x=PC1, y=PC2, color=Group)) +
    geom_point(size=3) +
    stat_ellipse(aes(group=Group), type="t", linetype=2) +
    theme_bw() +
    labs(title="PCA of Enzyme Commission (EC) Functional Profiles",
         x=paste0("PC1 (", round(pc_var[1], 1), "%)"),
         y=paste0("PC2 (", round(pc_var[2], 1), "%)")) +
    theme(legend.position="bottom")
  
  print(p_pca)
  
  # Top EC numbers by abundance
  ec_means <- rowMeans(t(ec_rel))
  top_ecs <- sort(ec_means, decreasing=TRUE)[1:20]
  
  ec_df <- data.frame(
    ec = names(top_ecs),
    mean_abundance = top_ecs
  )
  
  p_top_ec <- ggplot(ec_df, aes(x=reorder(ec, mean_abundance), y=mean_abundance)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(title="Top 20 EC Numbers by Mean Abundance",
         x="EC Number", y="Mean Relative Abundance")
  
  print(p_top_ec)
  
  # Top EC numbers by group
  top_ec_ids <- names(top_ecs)
  top_ec_data <- ec_rel[, top_ec_ids, drop=FALSE]
  
  top_ec_melt <- melt(top_ec_data)
  colnames(top_ec_melt) <- c("Sample", "EC", "Abundance")
  top_ec_melt$Group <- groups[match(top_ec_melt$Sample, rownames(ec_rel))]
  
  # Calculate mean abundance by EC and group
  ec_group_means <- top_ec_melt %>%
    group_by(EC, Group) %>%
    summarize(Mean = mean(Abundance), .groups='drop')
  
  # Create heatmap of EC abundance by group
  ec_matrix <- ec_group_means %>%
    pivot_wider(names_from=Group, values_from=Mean) %>%
    column_to_rownames("EC") %>%
    as.matrix()
  
  # Scale for better visualization
  ec_matrix_scaled <- scale(ec_matrix)
  
  pheatmap(ec_matrix_scaled, main="EC Abundance by Group (z-score)",
           fontsize_row=8, fontsize_col=10, scale="none")
  
  # Differential abundance analysis similar to KO data
  # (Code follows same pattern as KO analysis)
  
  # NMDS analysis for EC data
  tryCatch({
    ec_nmds <- metaMDS(vegdist(ec_rel, method="bray"), k=2)
    
    # Extract NMDS scores
    nmds_scores <- as.data.frame(ec_nmds$points)
    nmds_scores$Group <- groups
    
    p_nmds <- ggplot(nmds_scores, aes(x=MDS1, y=MDS2, color=Group)) +
      geom_point(size=3) +
      stat_ellipse(aes(group=Group), type="t", linetype=2) +
      theme_bw() +
      labs(title="NMDS of EC Functional Profiles",
           x="NMDS1", y="NMDS2",
           subtitle=paste("Stress =", round(ec_nmds$stress, 3))) +
      theme(legend.position="bottom")
    
    print(p_nmds)
  }, error=function(e) {
    warning("NMDS analysis failed for EC data: ", e$message)
  })
}

# Process and analyze pathway data
if (!is.null(pathway_data)) {
  cat("Analyzing MetaCyc pathway predictions...\n")
  
  # Match metadata to pathway data
  matched_pathway <- fix_metadata_for_data(metadata, pathway_data)
  pathway_meta <- matched_pathway$meta
  pathway_data <- matched_pathway$data
  
  # Remove features with zero variance
  pathway_var <- apply(pathway_data, 2, var)
  pathway_data <- pathway_data[, pathway_var > 0, drop=FALSE]
  
  # Filter low-abundance pathways
  pathway_means <- colMeans(pathway_data)
  pathway_data <- pathway_data[, pathway_means > quantile(pathway_means, 0.25), drop=FALSE]
  
  # Normalize data (relative abundance)
  pathway_rel <- t(apply(pathway_data, 1, function(x) x / sum(x)))
  
  # Get group information
  groups <- as.factor(pathway_meta[[opt$group_col]])
  
  # Perform NMDS
  tryCatch({
    pathway_nmds <- metaMDS(vegdist(pathway_rel, method="bray"), k=2)
    
    # Extract NMDS scores
    nmds_scores <- as.data.frame(pathway_nmds$points)
    nmds_scores$Group <- groups
    
    p_nmds <- ggplot(nmds_scores, aes(x=MDS1, y=MDS2, color=Group)) +
      geom_point(size=3) +
      stat_ellipse(aes(group=Group), type="t", linetype=2) +
      theme_bw() +
      labs(title="NMDS of MetaCyc Pathway Profiles",
           x="NMDS1", y="NMDS2",
           subtitle=paste("Stress =", round(pathway_nmds$stress, 3))) +
      theme(legend.position="bottom")
    
    print(p_nmds)
  }, error=function(e) {
    warning("NMDS analysis failed for pathway data: ", e$message)
  })
  
  # Top pathways by abundance
  pathway_means <- rowMeans(t(pathway_rel))
  top_pathways <- sort(pathway_means, decreasing=TRUE)[1:20]
  
  pathway_df <- data.frame(
    pathway = names(top_pathways),
    mean_abundance = top_pathways
  )
  
  # Clean pathway names
  pathway_df$pathway_short <- gsub("^.+\\|", "", pathway_df$pathway)
  pathway_df$pathway_short <- gsub("\\-PWY$", "", pathway_df$pathway_short)
  pathway_df$pathway_short <- gsub("\\_", " ", pathway_df$pathway_short)
  
  # If pathway names are too long, truncate them
  if(max(nchar(pathway_df$pathway_short)) > 30) {
    pathway_df$pathway_short <- substr(pathway_df$pathway_short, 1, 30)
    pathway_df$pathway_short <- paste0(pathway_df$pathway_short, "...")
  }
  
  p_top_pathway <- ggplot(pathway_df, aes(x=reorder(pathway_short, mean_abundance), y=mean_abundance)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(title="Top 20 MetaCyc Pathways by Mean Abundance",
         x="Pathway", y="Mean Relative Abundance")
  
  print(p_top_pathway)
  
  # Differential abundance analysis
  # Round counts for DESeq2
  pathway_counts <- round(pathway_data)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = pathway_counts,
    colData = data.frame(group=groups, row.names=rownames(pathway_counts)),
    design = ~group
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get number of groups
  group_levels <- levels(groups)
  
  if(length(group_levels) == 2) {
    # Simple pairwise comparison for 2 groups
    res <- results(dds)
    
    # Identify significant features
    sig_features <- which(res$padj < 0.05 & !is.na(res$padj))
    
    if(length(sig_features) > 0) {
      # Sort by significance
      sig_res <- res[sig_features,]
      sig_res <- sig_res[order(sig_res$padj),]
      
      # Create volcano plot
      res_df <- as.data.frame(res)
      res_df$feature <- rownames(res_df)
      res_df$significant <- ifelse(res_df$padj < 0.05 & !is.na(res_df$padj), "Yes", "No")
      
      # Clean pathway names for plotting
      res_df$pathway_short <- gsub("^.+\\|", "", res_df$feature)
      res_df$pathway_short <- gsub("\\-PWY$", "", res_df$pathway_short)
      res_df$pathway_short <- gsub("\\_", " ", res_df$pathway_short)
      
      p_volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
        geom_point(alpha=0.7) +
        scale_color_manual(values=c("No"="gray", "Yes"="red")) +
        geom_hline(yintercept=-log10(0.05), linetype="dashed") +
        geom_vline(xintercept=c(-1, 1), linetype="dashed") +
        theme_bw() +
        labs(title="Differential Abundance Analysis of Pathways",
             x="Log2 Fold Change",
             y="-Log10 Adjusted P-value")
      
      print(p_volcano)
      
      # Create a labeled volcano plot with pathway names
      top_sig <- head(res_df[res_df$significant == "Yes", ], 10)
      
      p_volcano_labeled <- p_volcano +
        geom_text_repel(data=top_sig, aes(label=pathway_short), size=3, box.padding=0.5)
      
      print(p_volcano_labeled)
      
      # Create heatmap of significant pathways
      sig_pathway_ids <- rownames(sig_res)
      sig_pathway_data <- pathway_rel[, sig_pathway_ids, drop=FALSE]
      
      # Clean pathway names for heatmap
      sig_pathway_names <- gsub("^.+\\|", "", sig_pathway_ids)
      sig_pathway_names <- gsub("\\-PWY$", "", sig_pathway_names)
      sig_pathway_names <- gsub("\\_", " ", sig_pathway_names)
      
      # Truncate long names
      if(max(nchar(sig_pathway_names)) > 40) {
        sig_pathway_names <- substr(sig_pathway_names, 1, 40)
        sig_pathway_names <- paste0(sig_pathway_names, "...")
      }
      
      colnames(sig_pathway_data) <- sig_pathway_names
      
      # Scale data for better visualization
      sig_pathway_data_scaled <- t(scale(t(sig_pathway_data)))
      
      # Create annotation data frame
      annotation_df <- data.frame(Group=groups, row.names=rownames(sig_pathway_data))
      
      # Create color palette for groups
      group_colors <- brewer.pal(max(3, length(unique(groups))), "Set1")[1:length(unique(groups))]
      names(group_colors) <- unique(groups)
      annotation_colors <- list(Group=group_colors)
      
      # Create heatmap
      pheatmap(sig_pathway_data_scaled, annotation_row=annotation_df, annotation_colors=annotation_colors,
               fontsize_row=8, fontsize_col=8, main="Significant Pathways",
               cluster_cols=TRUE, cluster_rows=TRUE)
      
      # Save significant pathways
      sig_res_df <- as.data.frame(sig_res)
      sig_res_df$pathway_id <- rownames(sig_res_df)
      sig_res_df$pathway_name <- sig_pathway_names
      write.csv(sig_res_df, file=file.path(opt$output_dir, "significant_pathways.csv"), row.names=FALSE)
    } else {
      grid.newpage()
      grid.text("No significant differentially abundant pathways found", x=0.5, y=0.5, gp=gpar(cex=1.2))
    }
  } else {
    # Multiple group comparisons - same pattern as for KO analysis
    # Implement pairwise comparisons between groups
    result_list <- list()
    
    for(i in 1:(length(group_levels)-1)) {
      for(j in (i+1):length(group_levels)) {
        group1 <- group_levels[i]
        group2 <- group_levels[j]
        
        # Subset data for these two groups
        idx <- groups %in% c(group1, group2)
        if(sum(idx) < 3) next  # Skip if fewer than 3 samples
        
        sub_dds <- DESeqDataSetFromMatrix(
          countData = pathway_counts[idx, ],
          colData = data.frame(group=groups[idx], row.names=rownames(pathway_counts)[idx]),
          design = ~group
        )
        
        # Run DESeq2
        sub_dds <- DESeq(sub_dds)
        
        # Get results
        res <- results(sub_dds, contrast=c("group", group2, group1))
        
        # Find significant features
        sig_idx <- which(res$padj < 0.05 & !is.na(res$padj))
        
        if(length(sig_idx) > 0) {
          sig_res <- res[sig_idx,]
          sig_res <- as.data.frame(sig_res[order(sig_res$padj),])
          sig_res$comparison <- paste(group2, "vs", group1)
          sig_res$feature <- rownames(sig_res)
          
          # Clean pathway names
          sig_res$pathway_name <- gsub("^.+\\|", "", sig_res$feature)
          sig_res$pathway_name <- gsub("\\-PWY$", "", sig_res$pathway_name)
          sig_res$pathway_name <- gsub("\\_", " ", sig_res$pathway_name)
          
          result_list[[length(result_list) + 1]] <- sig_res
        }
      }
    }
    
    if(length(result_list) > 0) {
      all_results <- do.call(rbind, result_list)
      
      # Save results
      write.csv(all_results, file=file.path(opt$output_dir, "pathway_pairwise_comparisons.csv"), row.names=FALSE)
      
      # Create summary plot
      p_summary <- ggplot(all_results, aes(x=comparison, y=pathway_name, fill=log2FoldChange)) +
        geom_tile() +
        scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
        theme_bw() +
        theme(axis.text.y=element_text(size=8),
              axis.text.x=element_text(angle=45, hjust=1)) +
        labs(title="Significant Pathways by Comparison",
             x="Comparison", y="Pathway", fill="Log2 Fold Change")
      
      print(p_summary)
    } else {
      grid.newpage()
      grid.text("No significant differentially abundant pathways found in any pairwise comparison", 
                x=0.5, y=0.5, gp=gpar(cex=1.2))
    }
  }
  
  # Pathway enrichment visualization
  # Create a ridge plot of pathway distributions across groups
  if(length(sig_features) > 0 && length(sig_features) < 30) {
    top_pathways <- rownames(sig_res)[1:min(length(sig_features), 15)]
    
    # Prepare data for ridge plot
    ridge_data <- pathway_rel[, top_pathways, drop=FALSE]
    ridge_df <- melt(ridge_data)
    colnames(ridge_df) <- c("Sample", "Pathway", "Abundance")
    ridge_df$Group <- groups[match(ridge_df$Sample, rownames(pathway_rel))]
    
    # Clean pathway names
    ridge_df$Pathway_name <- gsub("^.+\\|", "", ridge_df$Pathway)
    ridge_df$Pathway_name <- gsub("\\-PWY$", "", ridge_df$Pathway_name)
    ridge_df$Pathway_name <- gsub("\\_", " ", ridge_df$Pathway_name)
    
    # Create ridge plot
    p_ridge <- ggplot(ridge_df, aes(x=Abundance, y=Pathway_name, fill=Group)) +
      geom_density_ridges(alpha=0.7, scale=0.9, quantile_lines=TRUE, quantiles=2) +
      theme_bw() +
      labs(title="Distribution of Top Significant Pathways",
           x="Relative Abundance", y="Pathway")
    
    print(p_ridge)
  }
}

# Integrated analysis (if multiple data types are available)
if(!is.null(ko_data) && !is.null(pathway_data)) {
  cat("Performing integrated analysis of KO and pathway data...\n")
  
  # Find common samples
  common_samples <- intersect(rownames(ko_rel), rownames(pathway_rel))
  
  if(length(common_samples) > 0) {
    # Get common group information
    group_info <- groups[match(common_samples, rownames(pathway_rel))]
    
    # Procrustes analysis to compare KO and pathway ordinations
    ko_pca <- prcomp(ko_rel[common_samples, ], scale=TRUE)
    pathway_pca <- prcomp(pathway_rel[common_samples, ], scale=TRUE)
    
    # Extract first 2 PCs
    ko_coords <- ko_pca$x[, 1:2]
    pathway_coords <- pathway_pca$x[, 1:2]
    
    # Perform Procrustes analysis
    proc <- protest(ko_coords, pathway_coords)
    
    # Create plot data
    procrustes_data <- data.frame(
      ko1 = ko_coords[, 1],
      ko2 = ko_coords[, 2],
      pathway1 = pathway_coords[, 1] * proc$scale,
      pathway2 = pathway_coords[, 2] * proc$scale,
      sample = common_samples,
      group = group_info
    )
    
    # Create Procrustes plot
    p_proc <- ggplot(procrustes_data) +
      geom_point(aes(x=ko1, y=ko2, color=group), size=3) +
      geom_point(aes(x=pathway1, y=pathway2, color=group), size=3, shape=17) +
      geom_segment(aes(x=ko1, y=ko2, xend=pathway1, yend=pathway2, color=group), alpha=0.5) +
      theme_bw() +
      labs(title="Procrustes Analysis of KO and Pathway Ordinations",
           subtitle=paste("Correlation =", round(sqrt(1 - proc$ss), 3), 
                          ", p-value =", round(proc$signif, 3)))
    
    print(p_proc)
    
    # Create a legend explaining the symbols
    grid.newpage()
    grid.text("Procrustes Analysis Legend", x=0.5, y=0.9, gp=gpar(cex=1.2))
    grid.text("Circles: KO data ordination", x=0.5, y=0.7, gp=gpar(cex=1))
    grid.text("Triangles: Pathway data ordination", x=0.5, y=0.6, gp=gpar(cex=1))
    grid.text("Lines: Connect corresponding samples", x=0.5, y=0.5, gp=gpar(cex=1))
    grid.text(paste("Procrustes correlation:", round(sqrt(1 - proc$ss), 3)), x=0.5, y=0.3, gp=gpar(cex=1))
    grid.text(paste("p-value:", round(proc$signif, 3)), x=0.5, y=0.2, gp=gpar(cex=1))
  }
}

# Create summary plots
# Function distribution by group
if(!is.null(pathway_data)) {
  # Calculate mean abundance by pathway category and group
  # Extract pathway categories from names (first part before |)
  pathway_categories <- gsub("\\|.+$", "", colnames(pathway_data))
  
  # Create data frame with pathway categories
  category_data <- data.frame(
    pathway = colnames(pathway_data),
    category = pathway_categories,
    stringsAsFactors = FALSE
  )
  
  # Calculate mean abundance by category and group
  category_abundance <- data.frame(row.names=rownames(pathway_rel))
  
  for(cat in unique(pathway_categories)) {
    cat_pathways <- category_data$pathway[category_data$category == cat]
    if(length(cat_pathways) > 0) {
      category_abundance[[cat]] <- rowSums(pathway_rel[, cat_pathways, drop=FALSE])
    }
  }
  
  # Add group information
  category_abundance$Group <- groups
  
  # Melt for plotting
  category_melt <- melt(category_abundance, id.vars="Group", variable.name="Category", value.name="Abundance")
  
  # Calculate mean by group and category
  category_means <- category_melt %>%
    group_by(Group, Category) %>%
    summarize(Mean = mean(Abundance), .groups='drop')
  
  # Stacked bar plot
  p_category <- ggplot(category_means, aes(x=Group, y=Mean, fill=Category)) +
    geom_bar(stat="identity", position="stack") +
    theme_bw() +
    labs(title="Pathway Categories by Group",
         x="Group", y="Mean Relative Abundance")
  
  print(p_category)
}

# Create final summary page
grid.newpage()
grid.text("PICRUSt2 Functional Analysis Summary", x=0.5, y=0.95, gp=gpar(cex=1.5))

summary_text <- paste(
  "Data summary:",
  paste("- KO predictions:", ifelse(!is.null(ko_data), "Analyzed", "Not provided")),
  paste("- EC predictions:", ifelse(!is.null(ec_data), "Analyzed", "Not provided")),
  paste("- Pathway predictions:", ifelse(!is.null(pathway_data), "Analyzed", "Not provided")),
  paste("- Number of samples:", ifelse(!is.null(pathway_data), nrow(pathway_data), 
                               ifelse(!is.null(ko_data), nrow(ko_data), 
                                      ifelse(!is.null(ec_data), nrow(ec_data), "Unknown")))),
  paste("- Group variable:", opt$group_col),
  paste("- Number of groups:", length(unique(groups))),
  "",
  "Key findings:",
  ifelse(!is.null(ko_data) && exists("sig_features") && length(sig_features) > 0,
         paste("- Found", length(sig_features), "differentially abundant KO features"),
         "- No significant KO features found"),
  ifelse(!is.null(pathway_data) && exists("sig_features") && length(sig_features) > 0,
         paste("- Found", length(sig_features), "differentially abundant pathways"),
         "- No significant pathways found"),
  sep="\n"
)

grid.text(summary_text, x=0.5, y=0.5, gp=gpar(cex=1))

# Close PDF
dev.off()

# Create output files list
output_files <- list.files(opt$output_dir, full.names=TRUE)
cat("Analysis complete! Output files:\n")
cat(paste(" -", output_files), sep="\n")

