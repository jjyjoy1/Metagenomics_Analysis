#!/usr/bin/env Rscript

# Script for generating microbiome sample correlation network analysis from Phyloseq data
# This focuses specifically on network analysis to complement the main phyloseq script

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(biomformat)
  library(igraph)
  library(ggplot2)
  library(network)
  library(sna)
  library(GGally)
  library(RColorBrewer)
  library(intergraph)
})

# Parse command line arguments
option_list <- list(
  make_option("--biom", type="character", help="Path to BIOM file"),
  make_option("--taxonomy", type="character", help="Path to taxonomy file"),
  make_option("--metadata", type="character", help="Path to metadata file"),
  make_option("--method", type="character", default="spearman", help="Correlation method (spearman or pearson)"),
  make_option("--threshold", type="double", default=0.7, help="Correlation threshold"),
  make_option("--output", type="character", help="Output file path for network plot")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory if needed
output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
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

# Load metadata
metadata <- read.table(opt$metadata, sep="\t", header=TRUE, row.names=1, comment.char="")
SAMP <- sample_data(metadata)

# Create phyloseq object
ps <- phyloseq(OTU, TAX, SAMP)

# Normalize to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Filter to keep more abundant taxa
min_prevalence <- 0.1  # Present in at least 10% of samples
min_abundance <- 0.001  # At least 0.1% relative abundance
ps_filt <- filter_taxa(ps_rel, function(x) sum(x > min_abundance) > (min_prevalence * length(x)), TRUE)

cat("Creating sample correlation network analysis...\n")

# Get metadata for coloring nodes
metadata_df <- data.frame(sample_data(ps_filt))

# Extract the first available grouping column from metadata if not specified
if(ncol(metadata_df) > 0) {
  group_column <- colnames(metadata_df)[1]
  group_data <- metadata_df[, group_column]
  sample_groups <- as.factor(group_data)
} else {
  sample_groups <- rep("unknown", nsamples(ps_filt))
}

# Generate sample-sample network
# Get OTU table
otu_df <- as.data.frame(t(otu_table(ps_filt)))

# Calculate correlations between samples
sample_correlations <- cor(otu_df, method = opt$method)

# Apply threshold
sample_correlations_thresholded <- sample_correlations
sample_correlations_thresholded[abs(sample_correlations_thresholded) < opt$threshold] <- 0

# Create network from correlation matrix
sample_network <- graph_from_adjacency_matrix(
  sample_correlations_thresholded,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

# Remove any isolated nodes (samples with no connections)
sample_network <- delete.vertices(sample_network, which(degree(sample_network) == 0))

# Check if network has any edges
if (ecount(sample_network) == 0) {
  cat("Warning: No edges in the network after filtering with threshold", opt$threshold, "\n")
  cat("Try lowering the threshold value.\n")
  
  # Create a simple warning plot
  pdf(opt$output, width=10, height=8)
  plot(1, 1, type="n", xlab="", ylab="", main="Empty Network (No connections above threshold)")
  text(1, 1, paste("No edges found with correlation threshold of", opt$threshold))
  text(1, 0.8, "Try lowering the threshold value")
  dev.off()
  
  quit(status=0)
}

# Add vertex attributes for visualization
V(sample_network)$name <- names(V(sample_network))
V(sample_network)$group <- sample_groups[match(V(sample_network)$name, names(sample_groups))]

# Get layout
set.seed(42)  # For reproducibility
layout <- layout_with_fr(sample_network)

# Create output PDF
pdf(opt$output, width=12, height=10)

# Plot network with ggnet2
net <- network::network(sample_correlations_thresholded, directed=FALSE)

# Check if we have group data 
group_attr_name <- "group"
if(length(sample_groups) > 0) {
  net %v% group_attr_name <- sample_groups[match(network.vertex.names(net), names(sample_groups))]
}

# Create network plot
ggnet2_plot <- ggnet2(
  net,
  mode = "fruchtermanreingold",
  size = 6,
  color = group_attr_name,
  palette = "Set1",
  edge.alpha = 0.5,
  edge.size = "weight",
  edge.color = "grey70",
  label = TRUE,
  label.size = 3
)

print(ggnet2_plot)

# Also create a plot using igraph for better control
plot(
  sample_network,
  vertex.color = rainbow(length(unique(V(sample_network)$group)))[as.numeric(as.factor(V(sample_network)$group))],
  vertex.label = V(sample_network)$name,
  vertex.label.cex = 0.8,
  vertex.size = 12,
  edge.width = abs(E(sample_network)$weight) * 3,
  edge.color = ifelse(E(sample_network)$weight > 0, "blue", "red"),
  layout = layout,
  main = paste("Sample Correlation Network (", opt$method, ", threshold =", opt$threshold, ")")
)

legend(
  "topright",
  legend = levels(as.factor(V(sample_network)$group)),
  col = rainbow(length(unique(V(sample_network)$group))),
  pch = 19,
  title = "Sample Groups",
  cex = 0.8
)

# Additional analysis: Create taxa co-occurrence network
cat("Creating taxa co-occurrence network...\n")

# Filter to genus level if available
if("Genus" %in% colnames(tax_table(ps_filt))) {
  ps_genus <- tax_glom(ps_filt, taxrank="Genus")
  
  # Get abundance table
  genus_abund <- as.data.frame(otu_table(ps_genus))
  
  # Get taxonomic information for labeling
  taxa_info <- tax_table(ps_genus)[, "Genus"]
  rownames(genus_abund) <- taxa_info
  
  # Calculate correlations
  genus_cor <- cor(t(genus_abund), method=opt$method)
  
  # Apply threshold
  genus_cor_threshold <- genus_cor
  genus_cor_threshold[abs(genus_cor_threshold) < opt$threshold] <- 0
  
  # Remove self-correlations
  diag(genus_cor_threshold) <- 0
  
  # Create network
  taxa_network <- graph_from_adjacency_matrix(
    genus_cor_threshold,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  # Remove isolated taxa
  taxa_network <- delete.vertices(taxa_network, which(degree(taxa_network) == 0))
  
  # Add taxonomic information
  V(taxa_network)$phylum <- tax_table(ps_genus)[V(taxa_network)$name, "Phylum"]
  
  # Create layout
  set.seed(42)
  taxa_layout <- layout_with_fr(taxa_network)
  
  # Plot the taxa network
  if(ecount(taxa_network) > 0) {
    # Color nodes by phylum
    phylum_colors <- rainbow(length(unique(V(taxa_network)$phylum)))
    node_colors <- phylum_colors[as.numeric(as.factor(V(taxa_network)$phylum))]
    
    # Plot taxa network
    plot(
      taxa_network,
      vertex.color = node_colors,
      vertex.label = V(taxa_network)$name,
      vertex.label.cex = 0.7,
      vertex.size = 8,
      edge.width = abs(E(taxa_network)$weight) * 2,
      edge.color = ifelse(E(taxa_network)$weight > 0, "blue", "red"),
      layout = taxa_layout,
      main = "Taxa Co-occurrence Network"
    )
    
    legend(
      "topright",
      legend = levels(as.factor(V(taxa_network)$phylum)),
      col = phylum_colors,
      pch = 19,
      title = "Phylum",
      cex = 0.7
    )
  } else {
    plot(1, 1, type="n", xlab="", ylab="", main="Empty Taxa Network (No connections above threshold)")
    text(1, 1, paste("No edges found with correlation threshold of", opt$threshold))
  }
}

# Network statistics
cat("Calculating network statistics...\n")

# Sample network statistics
sample_density <- edge_density(sample_network)
sample_transitivity <- transitivity(sample_network)
sample_diameter <- diameter(sample_network)
sample_average_path <- mean_distance(sample_network)
sample_modularity <- modularity(cluster_louvain(sample_network))

# Create a table of network statistics
stats_df <- data.frame(
  Metric = c("Number of nodes", "Number of edges", "Network density", 
             "Clustering coefficient", "Diameter", "Average path length",
             "Modularity"),
  Value = c(vcount(sample_network), ecount(sample_network), sample_density,
            sample_transitivity, sample_diameter, sample_average_path,
            sample_modularity)
)

# Plot the statistics table
grid.newpage()
grid.table(stats_df, rows=NULL)

# Close PDF device
dev.off()

# Save network data for further analysis
network_data <- list(
  sample_correlations = sample_correlations,
  sample_network = sample_network,
  threshold = opt$threshold,
  method = opt$method
)

saveRDS(network_data, file=paste0(tools::file_path_sans_ext(opt$output), "_data.rds"))

cat("Network analysis complete! Results saved to:", opt$output, "\n")
