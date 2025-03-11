#!/usr/bin/env Rscript

# Comprehensive Phyloseq analysis script for microbiome data
# This script generates various visualizations and analyses using phyloseq

# Load required libraries
suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(biomformat)
  library(optparse)
  library(ape)
  library(gridExtra)
  library(knitr)
  library(rmarkdown)
  library(ggrepel)
  library(RColorBrewer)
  library(microbiome)
  library(patchwork)
})

# Parse command line arguments
option_list <- list(
  make_option("--biom", type="character", help="Path to BIOM file"),
  make_option("--taxonomy", type="character", help="Path to taxonomy file"),
  make_option("--tree", type="character", help="Path to phylogenetic tree file (optional)", default=NULL),
  make_option("--metadata", type="character", help="Path to metadata file"),
  make_option("--group-col", type="character", help="Column name for group comparisons"),
  make_option("--output-dir", type="character", help="Output directory for results")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# Load data into phyloseq object
cat("Loading data into phyloseq...\n")

# Load biom file
biom_data <- read_biom(opt$biom)
otu_table <- as(biom_data, "matrix")
OTU <- otu_table(otu_table, taxa_are_rows = TRUE)

# Load taxonomy
tax_data <- read.table(opt$taxonomy, sep="\t", header=TRUE, row.names=1, comment.char="")
if("Taxon" %in% colnames(tax_data)) {
  # Parse taxonomy string
  tax_split <- strsplit(as.character(tax_data$Taxon), "; ")
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_matrix <- matrix(NA, nrow=nrow(tax_data), ncol=length(tax_levels))
  
  for(i in 1:nrow(tax_data)) {
    tax_vector <- tax_split[[i]]
    for(j in 1:length(tax_vector)) {
      if(j <= ncol(tax_matrix)) {
        tax_matrix[i,j] <- tax_vector[j]
      }
    }
  }
  
  rownames(tax_matrix) <- rownames(tax_data)
  colnames(tax_matrix) <- tax_levels
  TAX <- tax_table(tax_matrix)
} else {
  # Assume taxonomy is already in the right format
  TAX <- tax_table(as.matrix(tax_data))
}

# Load phylogenetic tree if provided
if(!is.null(opt$tree) && file.exists(opt$tree)) {
  tree <- read.tree(opt$tree)
  
  # Check if tree tip labels match OTU IDs
  tree_tips <- tree$tip.label
  otu_ids <- rownames(otu_table)
  shared_ids <- intersect(tree_tips, otu_ids)
  
  if(length(shared_ids) == 0) {
    warning("No matching IDs between tree and OTU table. Proceeding without tree.")
    tree <- NULL
  } else if(length(shared_ids) < length(otu_ids)) {
    warning(paste("Only", length(shared_ids), "of", length(otu_ids), "OTUs match tree tips. Pruning data."))
    OTU <- OTU[shared_ids, ]
    TAX <- TAX[shared_ids, ]
    tree <- prune_taxa(shared_ids, tree)
  }
}

# Load metadata
metadata <- read.table(opt$metadata, sep="\t", header=TRUE, row.names=1, comment.char="")
SAMP <- sample_data(metadata)

# Combine data into phyloseq object
if(!is.null(opt$tree) && exists("tree")) {
  ps <- phyloseq(OTU, TAX, SAMP, tree)
} else {
  ps <- phyloseq(OTU, TAX, SAMP)
}

# Check if OTU IDs match taxonomy IDs
missing_taxa <- setdiff(taxa_names(ps), rownames(TAX))
if(length(missing_taxa) > 0) {
  warning(paste("Some OTU IDs are missing from taxonomy. Removing", length(missing_taxa), "OTUs."))
  ps <- prune_taxa(setdiff(taxa_names(ps), missing_taxa), ps)
}

# Basic data summary
cat("Creating basic data summary...\n")
summary_file <- file.path(opt$output_dir, "data_summary.txt")
sink(summary_file)
cat("Number of samples:", nsamples(ps), "\n")
cat("Number of taxa:", ntaxa(ps), "\n")
cat("Taxonomic ranks:", paste(colnames(tax_table(ps)), collapse=", "), "\n")
cat("\nSample variables:\n")
print(sample_variables(ps))

cat("\nRead counts per sample:\n")
samp_sums <- sort(sample_sums(ps))
print(summary(samp_sums))
sink()

# Create alpha diversity plots
cat("Generating alpha diversity plots...\n")
alpha_div <- estimate_richness(ps, measures=c("Observed", "Shannon", "Simpson", "InvSimpson", "Chao1"))
alpha_div$Sample <- rownames(alpha_div)
alpha_div_long <- melt(alpha_div, id.vars="Sample")
alpha_div_long$Group <- sample_data(ps)[[opt$group_col]][match(alpha_div_long$Sample, rownames(sample_data(ps)))]

# Add Faith's PD if tree is present
if(!is.null(opt$tree) && exists("tree")) {
  pd <- pd(ps, include.root=FALSE)
  pd_df <- data.frame(Sample=names(pd), PD=pd)
  pd_long <- melt(pd_df, id.vars="Sample", variable.name="measure", value.name="value")
  pd_long$Group <- sample_data(ps)[[opt$group_col]][match(pd_long$Sample, rownames(sample_data(ps)))]
  alpha_div_long <- rbind(alpha_div_long, pd_long)
}

# Alpha diversity plots
p_alpha <- ggplot(alpha_div_long, aes(x=Group, y=value, color=Group)) +
  geom_boxplot(alpha=0.6) +
  geom_jitter(width=0.2) +
  facet_wrap(~variable, scales="free_y") +
  theme_bw() +
  labs(title="Alpha Diversity Measures", x="Group", y="Value") +
  theme(legend.position="bottom", 
        axis.text.x=element_text(angle=45, hjust=1))

ggsave(file.path(opt$output_dir, "alpha_diversity_boxplots.pdf"), p_alpha, width=10, height=8)

# Statistical testing for alpha diversity
alpha_stats_file <- file.path(opt$output_dir, "alpha_diversity_stats.txt")
sink(alpha_stats_file)
cat("Alpha Diversity Statistical Tests\n")
cat("================================\n\n")

for(measure in unique(alpha_div_long$variable)) {
  cat("Measure:", measure, "\n")
  subset_data <- alpha_div_long[alpha_div_long$variable == measure, ]
  
  # Check normality
  shapiro_test <- shapiro.test(subset_data$value)
  cat("Shapiro-Wilk normality test: W =", shapiro_test$statistic, ", p-value =", shapiro_test$p.value, "\n")
  
  # If more than 2 groups, use ANOVA or Kruskal-Wallis
  if(length(unique(subset_data$Group)) > 2) {
    if(shapiro_test$p.value > 0.05) {
      # Use ANOVA for normal data
      aov_result <- aov(value ~ Group, data = subset_data)
      cat("ANOVA results:\n")
      print(summary(aov_result))
      
      # Post-hoc test
      tukey_result <- TukeyHSD(aov_result)
      cat("Tukey's HSD post-hoc test:\n")
      print(tukey_result)
    } else {
      # Use Kruskal-Wallis for non-normal data
      kw_result <- kruskal.test(value ~ Group, data = subset_data)
      cat("Kruskal-Wallis test:\n")
      print(kw_result)
      
      # Post-hoc Dunn test
      if(kw_result$p.value < 0.05) {
        dunn_result <- dunnTest(value ~ Group, data = subset_data, method="bonferroni")
        cat("Dunn's post-hoc test:\n")
        print(dunn_result)
      }
    }
  } else {
    # For 2 groups, use t-test or Wilcoxon
    if(shapiro_test$p.value > 0.05) {
      # Use t-test for normal data
      t_result <- t.test(value ~ Group, data = subset_data)
      cat("t-test results:\n")
      print(t_result)
    } else {
      # Use Wilcoxon for non-normal data
      w_result <- wilcox.test(value ~ Group, data = subset_data)
      cat("Wilcoxon test results:\n")
      print(w_result)
    }
  }
  cat("\n\n")
}
sink()

# Beta diversity analysis
cat("Performing beta diversity analysis...\n")

# Ordination plots
set.seed(123)
ord_methods <- c("NMDS", "PCoA")
dist_methods <- c("bray", "jaccard")

# Create plots for each combination
for(ord_method in ord_methods) {
  for(dist_method in dist_methods) {
    # Skip if phylogeny is required but not available
    if(ord_method == "UniFrac" && (!exists("tree") || is.null(tree))) {
      next
    }
    
    tryCatch({
      if(dist_method %in% c("bray", "jaccard")) {
        ord <- ordinate(ps, method=ord_method, distance=dist_method)
      } else {
        ord <- ordinate(ps, method=ord_method, distance=dist_method)
      }
      
      # Create plot
      p_beta <- plot_ordination(ps, ord, color=opt$group_col) + 
        geom_point(size=3, alpha=0.7) + 
        theme_bw() +
        ggtitle(paste(ord_method, "using", dist_method, "distances"))
      
      # Add ellipses
      p_beta <- p_beta + 
        stat_ellipse(aes(group=get(opt$group_col)), type="t", linetype=2) +
        theme(legend.position="bottom")
      
      # Save plot
      filename <- paste0("beta_diversity_", tolower(ord_method), "_", dist_method, ".pdf")
      ggsave(file.path(opt$output_dir, filename), p_beta, width=8, height=6)
      
      # Statistical testing (PERMANOVA)
      perm_result <- adonis2(distance(ps, method=dist_method) ~ sample_data(ps)[[opt$group_col]], permutations=999)
      
      # Save stats
      stats_filename <- paste0("beta_diversity_", tolower(ord_method), "_", dist_method, "_stats.txt")
      sink(file.path(opt$output_dir, stats_filename))
      cat("PERMANOVA Results for", ord_method, "using", dist_method, "distances\n")
      cat("========================================================\n\n")
      print(perm_result)
      sink()
      
    }, error=function(e) {
      warning(paste("Error in ordination with", ord_method, "and", dist_method, ":", e$message))
    })
  }
}

# Taxonomic composition barplots
cat("Creating taxonomic composition barplots...\n")

# Function to create barplots at different taxonomic levels
create_taxa_barplot <- function(ps, level) {
  ps_glom <- tax_glom(ps, taxrank=level)
  ps_rel <- transform_sample_counts(ps_glom, function(x) x/sum(x))
  
  # Extract the top taxa
  top_taxa_count <- 20
  top_taxa <- names(sort(taxa_sums(ps_rel), decreasing=TRUE)[1:min(top_taxa_count, ntaxa(ps_rel))])
  ps_rel_top <- prune_taxa(top_taxa, ps_rel)
  
  # Melt data for plotting
  ps_melt <- psmelt(ps_rel_top)
  
  # Create a more pleasing color palette
  n_colors <- length(unique(ps_melt[[level]]))
  if(n_colors <= 12) {
    colors <- brewer.pal(max(n_colors, 3), "Set3")
    if(n_colors < 3) colors <- colors[1:n_colors]
  } else {
    colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
  }
  
  # Create barplot
  p <- ggplot(ps_melt, aes_string(x="Sample", y="Abundance", fill=level)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=colors) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
          legend.text=element_text(size=8),
          legend.key.size=unit(0.8, "lines")) +
    labs(x="Sample", y="Relative Abundance", title=paste("Taxonomic Composition at", level, "Level")) +
    facet_grid(~get(opt$group_col), scales="free_x", space="free")
  
  # Save plot
  filename <- paste0("taxa_barplot_", level, ".pdf")
  ggsave(file.path(opt$output_dir, filename), p, width=12, height=8)
  
  return(p)
}

# Create barplots for each taxonomic level
tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
for(level in tax_levels) {
  if(level %in% colnames(tax_table(ps))) {
    create_taxa_barplot(ps, level)
  }
}

# Create heatmap of top taxa
cat("Creating taxa heatmap...\n")

# Filter to keep only the most abundant taxa
top_taxa_count <- 30
rel_abund <- transform_sample_counts(ps, function(x) x/sum(x))
top_taxa <- names(sort(taxa_sums(rel_abund), decreasing=TRUE)[1:min(top_taxa_count, ntaxa(rel_abund))])
ps_top <- prune_taxa(top_taxa, rel_abund)

# Get taxonomic information
tax_info <- tax_table(ps_top)[, "Genus"]
if("Species" %in% colnames(tax_table(ps_top))) {
  tax_info <- paste(tax_info, tax_table(ps_top)[, "Species"])
}

# Create heatmap
otu_matrix <- as.matrix(otu_table(ps_top))
rownames(otu_matrix) <- tax_info

# Order samples by group
sample_groups <- as.factor(sample_data(ps_top)[[opt$group_col]])
sample_order <- order(sample_groups)

# Create heatmap using pheatmap
pdf(file.path(opt$output_dir, "heatmap_plots.pdf"), width=12, height=10)
pheatmap(otu_matrix[, sample_order], 
         scale="row", 
         clustering_distance_rows="correlation",
         clustering_method="average",
         annotation_col=data.frame(Group=sample_groups, row.names=names(sample_groups)),
         fontsize_row=8, 
         fontsize_col=8)
dev.off()

# Differential abundance analysis
cat("Performing differential abundance analysis with DESeq2...\n")

# Function to run DESeq2 differential abundance
run_deseq2 <- function(ps, group_var) {
  # Convert to DESeq2 format
  dds <- phyloseq_to_deseq2(ps, design = formula(paste0("~", group_var)))
  
  # Run DESeq2
  dds <- DESeq(dds, fitType="local")
  
  # Get results
  res <- results(dds)
  
  # Add taxonomy information
  res_df <- as.data.frame(res)
  res_df$OTU <- rownames(res_df)
  tax_table_df <- as.data.frame(tax_table(ps))
  res_df <- merge(res_df, tax_table_df, by.x="OTU", by.y="row.names")
  
  # Sort by p-value
  res_df <- res_df[order(res_df$padj, na.last=TRUE), ]
  
  return(res_df)
}

# Run DESeq2 if there are at least 2 groups
group_levels <- unique(sample_data(ps)[[opt$group_col]])
if(length(group_levels) >= 2) {
  deseq_res <- run_deseq2(ps, opt$group_col)
  
  # Save results
  write.csv(deseq_res, file.path(opt$output_dir, "deseq2_results.csv"), row.names=FALSE)
  
  # Create volcano plot
  if(nrow(deseq_res) > 0 && !all(is.na(deseq_res$padj))) {
    # Add labels for significant features
    deseq_res$Significant <- ifelse(deseq_res$padj < 0.05 & abs(deseq_res$log2FoldChange) > 1, "Yes", "No")
    deseq_res$Label <- ifelse(deseq_res$Significant == "Yes", 
                             paste(deseq_res$Genus, deseq_res$Species, sep=" "), 
                             "")
    
    # Create plot
    p_volcano <- ggplot(deseq_res, aes(x=log2FoldChange, y=-log10(padj), color=Significant, label=Label)) +
      geom_point(alpha=0.7) +
      scale_color_manual(values=c("No"="gray", "Yes"="red")) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed") +
      geom_vline(xintercept=c(-1, 1), linetype="dashed") +
      theme_bw() +
      labs(title="DESeq2 Differential Abundance",
           x="Log2 Fold Change",
           y="-Log10 Adjusted P-value") +
      geom_text_repel(size=3, max.overlaps=25)
    
    ggsave(file.path(opt$output_dir, "deseq2_volcano_plot.pdf"), p_volcano, width=10, height=8)
  }
}

# Network analysis
cat("Performing network analysis...\n")

# Filter data to keep more abundant taxa
min_prevalence <- 0.1  # Present in at least 10% of samples
min_abundance <- 0.001  # At least 0.1% relative abundance
ps_filt <- filter_taxa(rel_abund, function(x) sum(x > min_abundance) > (min_prevalence * length(x)), TRUE)

# Create network at genus level if available
if("Genus" %in% colnames(tax_table(ps))) {
  ps_genus <- tax_glom(ps_filt, taxrank="Genus")
  
  # Export network data
  otu_table_genus <- as.data.frame(t(otu_table(ps_genus)))
  write.csv(otu_table_genus, file.path(opt$output_dir, "genus_abundance_for_network.csv"))
  
  # Basic correlation network
  genus_cor <- cor(otu_table_genus, method="spearman")
  genus_cor[abs(genus_cor) < 0.7] <- 0  # Keep only strong correlations
  
  # Save correlation matrix
  write.csv(genus_cor, file.path(opt$output_dir, "genus_correlation_matrix.csv"))
  
  # Create simple network visualization
  pdf(file.path(opt$output_dir, "network_analysis.pdf"), width=12, height=10)
  net <- igraph::graph_from_adjacency_matrix(
    genus_cor, 
    mode="undirected", 
    weighted=TRUE, 
    diag=FALSE
  )
  
  # Remove isolated nodes
  net <- delete.vertices(net, which(degree(net) == 0))
  
  # Plot network
  set.seed(42)
  plot(net, 
       vertex.label=igraph::V(net)$name,
       vertex.label.cex=0.7,
       vertex.size=5,
       edge.width=abs(E(net)$weight) * 2,
       layout=layout_with_fr(net))
  dev.off()
}

# Create RMarkdown report
cat("Generating comprehensive report...\n")

# Create R markdown template
rmd_template <- '
---
title: "Phyloseq Microbiome Analysis Report"
date: "`r format(Sys.time(), "%B %d, %Y")`"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
library(ggplot2)
library(knitr)
library(DT)
```

# Overview

This report provides a comprehensive analysis of 16S rRNA microbiome data using the phyloseq R package.

## Dataset Summary

```{r}
summary_data <- read.table("data_summary.txt", sep="\n", header=FALSE, quote="")
cat(paste(summary_data$V1, collapse="\n"))
```

# Alpha Diversity

Alpha diversity represents the species diversity within each sample.

```{r, fig.width=10, fig.height=8}
knitr::include_graphics("alpha_diversity_boxplots.pdf")
```

## Statistical Analysis

```{r}
alpha_stats <- read.table("alpha_diversity_stats.txt", sep="\n", header=FALSE, quote="")
cat(paste(alpha_stats$V1, collapse="\n"))
```

# Beta Diversity

Beta diversity represents the comparison of microbial communities between samples.

```{r, fig.width=8, fig.height=6}
beta_files <- list.files(pattern="beta_diversity_.*\\.pdf$")
for(file in beta_files) {
  if(!grepl("stats", file)) {
    knitr::include_graphics(file)
  }
}
```

## PERMANOVA Results

```{r}
for(file in list.files(pattern="beta_diversity_.*_stats.txt$")) {
  cat(paste0("### ", gsub("_stats.txt", "", gsub("beta_diversity_", "", file)), "\n\n"))
  stats_data <- read.table(file, sep="\n", header=FALSE, quote="")
  cat(paste(stats_data$V1, collapse="\n"))
  cat("\n\n")
}
```

# Taxonomic Composition

```{r, fig.width=12, fig.height=8}
for(file in list.files(pattern="taxa_barplot_.*\\.pdf$")) {
  level <- gsub("taxa_barplot_", "", gsub(".pdf", "", file))
  cat(paste0("## ", level, " Level\n\n"))
  knitr::include_graphics(file)
  cat("\n\n")
}
```

# Heatmap of Top Taxa

```{r, fig.width=12, fig.height=10}
knitr::include_graphics("heatmap_plots.pdf")
```

# Differential Abundance

Differential abundance analysis identifies taxa that significantly differ between groups.

```{r}
if(file.exists("deseq2_results.csv")) {
  deseq_res <- read.csv("deseq2_results.csv")
  
  # Display significant results (if any)
  sig_results <- deseq_res[!is.na(deseq_res$padj) & deseq_res$padj < 0.05, ]
  
  if(nrow(sig_results) > 0) {
    cat("## Significant Differentially Abundant Taxa\n\n")
    datatable(sig_results[, c("OTU", "Phylum", "Class", "Order", "Family", "Genus", "log2FoldChange", "padj")], 
              options = list(pageLength = 10, scrollX = TRUE))
    
    cat("\n\n## Volcano Plot\n\n")
    if(file.exists("deseq2_volcano_plot.pdf")) {
      knitr::include_graphics("deseq2_volcano_plot.pdf")
    }
  } else {
    cat("No significantly differentially abundant taxa were found (adjusted p-value < 0.05).\n")
  }
} else {
  cat("Differential abundance analysis results not available.\n")
}
```

# Network Analysis

Correlation network showing relationships between taxa.

```{r, fig.width=12, fig.height=10}
if(file.exists("network_analysis.pdf")) {
  knitr::include_graphics("network_analysis.pdf")
  
  cat("\n\nThis network visualization shows taxa (nodes) with strong correlations (edges). Positive correlations indicate taxa that tend to occur together, while negative correlations indicate taxa that tend to exclude each other.\n")
} else {
  cat("Network analysis results not available.\n")
}
```

# Conclusion

This comprehensive analysis provides insights into the microbial community structure, diversity patterns, and differentially abundant taxa in the dataset. The results can be further interpreted in the context of the specific biological questions and experimental design.
'

# Write RMarkdown template to file
rmd_file <- file.path(opt$output_dir, "phyloseq_analysis.Rmd")
writeLines(rmd_template, rmd_file)

# Render RMarkdown to HTML
rmarkdown::render(rmd_file, 
                  output_file = "phyloseq_analysis_report.html",
                  output_dir = opt$output_dir,
                  quiet = TRUE)

cat("Analysis complete! Results are available in:", opt$output_dir, "\n")
