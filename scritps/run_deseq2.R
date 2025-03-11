#!/usr/bin/env Rscript

# Script to perform DESeq2 analysis on QIIME2 feature table
# Usage: Rscript run_deseq2.R biom_file metadata_file group_column control_group output_file

# Load required libraries
library(biomformat)
library(DESeq2)
library(phyloseq)
library(tidyverse)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Please provide all required arguments:
       1. BIOM file path
       2. Metadata file path
       3. Group column name
       4. Control group name
       5. Output file path")
}

biom_file <- args[1]
metadata_file <- args[2]
group_col <- args[3]
control_group <- args[4]
output_file <- args[5]

# Read biom file
biom_data <- read_biom(biom_file)
otu_table <- as.matrix(biom_data_matrix(biom_data))

# Read metadata
metadata <- read.table(metadata_file, header = TRUE, sep = "\t", row.names = 1)

# Ensure samples match between OTU table and metadata
shared_samples <- intersect(colnames(otu_table), rownames(metadata))
if (length(shared_samples) == 0) {
  stop("No matching samples between feature table and metadata")
}

otu_table <- otu_table[, shared_samples]
metadata <- metadata[shared_samples, , drop = FALSE]

# Create phyloseq object
OTU <- otu_table(otu_table, taxa_are_rows = TRUE)
META <- sample_data(metadata)
sample_names(META) <- rownames(metadata)
ps <- phyloseq(OTU, META)

# Convert to DESeq2 object
group_var <- metadata[[group_col]]
deseq_object <- phyloseq_to_deseq2(ps, design = as.formula(paste0("~", group_col)))

# Run DESeq2
deseq_object <- DESeq(deseq_object)

# Get results for each group compared to the control
groups <- unique(group_var)
groups <- groups[groups != control_group]

results_list <- list()
for (group in groups) {
  # Get results
  res <- results(deseq_object, contrast = c(group_col, group, control_group))
  res_df <- as.data.frame(res)
  res_df$feature_id <- rownames(res_df)
  res_df$comparison <- paste(group, "vs", control_group)
  
  # Add to list
  results_list[[length(results_list) + 1]] <- res_df
}

# Combine all results
all_results <- bind_rows(results_list)

# Clean up p-values
all_results$padj[is.na(all_results$padj)] <- 1

# Write results to file
write.table(all_results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("DESeq2 analysis completed. Results written to:", output_file, "\n")


