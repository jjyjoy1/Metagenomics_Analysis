# scripts/calculate_viral_richness.R
# Script to calculate viral richness from Kraken2 reports and viral contigs

# Load required libraries
library(tidyverse)
library(vegan)
library(reshape2)
library(ggplot2)
library(viridis)

# Function to parse Kraken2 reports
parse_kraken_reports <- function(filenames) {
  sample_names <- basename(filenames) %>%
    str_remove("_kraken2_report.txt")
  
  all_taxa <- data.frame()
  
  for (i in seq_along(filenames)) {
    file <- filenames[i]
    sample <- sample_names[i]
    
    # Read Kraken2 report
    kraken_data <- read.delim(file, header = FALSE, 
                             col.names = c("percentage", "clade_reads", "direct_reads", 
                                          "rank", "taxid", "name"))
    
    # Filter for viral species
    viral_species <- kraken_data %>%
      filter(rank == "S" & str_detect(name, "Viruses"))
    
    if (nrow(viral_species) > 0) {
      viral_species$sample <- sample
      all_taxa <- bind_rows(all_taxa, viral_species)
    }
  }
  
  return(all_taxa)
}

# Function to calculate viral contig statistics
calculate_contig_stats <- function(fasta_files) {
  sample_names <- basename(dirname(fasta_files))
  
  contig_stats <- data.frame(sample = character(),
                            contig_count = numeric(),
                            total_length = numeric(),
                            min_length = numeric(),
                            max_length = numeric(),
                            mean_length = numeric(),
                            n50 = numeric(),
                            stringsAsFactors = FALSE)
  
  for (i in seq_along(fasta_files)) {
    file <- fasta_files[i]
    sample <- sample_names[i]
    
    # Read FASTA headers to get lengths
    fasta_headers <- system(paste("grep '>' ", file), intern = TRUE)
    
    if (length(fasta_headers) > 0) {
      # Extract lengths if available in headers, otherwise count sequences
      has_length <- str_detect(fasta_headers[1], "length=")
      
      if (has_length) {
        lengths <- str_extract(fasta_headers, "length=\\d+") %>%
          str_remove("length=") %>%
          as.numeric()
      } else {
        # Calculate lengths by reading the file
        system(paste("bioawk -c fastx '{print $name, length($seq)}' ", file, " > temp_lengths.txt"))
        length_data <- read.delim("temp_lengths.txt", header = FALSE, col.names = c("header", "length"))
        lengths <- length_data$length
        file.remove("temp_lengths.txt")
      }
      
      # Calculate statistics
      contig_count <- length(lengths)
      total_length <- sum(lengths)
      min_length <- min(lengths)
      max_length <- max(lengths)
      mean_length <- mean(lengths)
      
      # Calculate N50
      sorted_lengths <- sort(lengths, decreasing = TRUE)
      cumulative_length <- cumsum(sorted_lengths)
      n50 <- sorted_lengths[which(cumulative_length >= total_length/2)[1]]
      
      new_row <- data.frame(sample = sample,
                           contig_count = contig_count,
                           total_length = total_length,
                           min_length = min_length,
                           max_length = max_length,
                           mean_length = mean_length,
                           n50 = n50)
      
      contig_stats <- rbind(contig_stats, new_row)
    } else {
      # No contigs found
      new_row <- data.frame(sample = sample,
                           contig_count = 0,
                           total_length = 0,
                           min_length = NA,
                           max_length = NA,
                           mean_length = NA,
                           n50 = NA)
      
      contig_stats <- rbind(contig_stats, new_row)
    }
  }
  
  return(contig_stats)
}

# Function to calculate richness indices
calculate_richness <- function(taxa_data) {
  # Create abundance matrix for vegan
  abundance_matrix <- taxa_data %>%
    select(name, direct_reads, sample) %>%
    pivot_wider(names_from = name, values_from = direct_reads, values_fill = 0) %>%
    column_to_rownames("sample")
  
  # Calculate richness indices
  richness <- data.frame(
    sample = rownames(abundance_matrix),
    observed_species = specnumber(abundance_matrix),
    shannon = diversity(abundance_matrix, index = "shannon"),
    simpson = diversity(abundance_matrix, index = "simpson"),
    invsimpson = diversity(abundance_matrix, index = "invsimpson"),
    chao1 = estimateR(abundance_matrix)[2,],
    ACE = estimateR(abundance_matrix)[4,]
  )
  
  return(richness)
}

# Function to plot richness results
plot_richness <- function(richness_data, contig_stats) {
  # Combine data
  combined_data <- richness_data %>%
    left_join(contig_stats, by = "sample")
  
  # Create plot
  richness_plot <- ggplot(combined_data, aes(x = reorder(sample, observed_species))) +
    geom_bar(aes(y = observed_species), stat = "identity", fill = "#56B4E9", alpha = 0.7) +
    geom_point(aes(y = shannon * 10), color = "#D55E00", size = 3) +
    geom_text(aes(y = observed_species + 2, label = contig_count), 
              size = 3, vjust = 0) +
    scale_y_continuous(
      name = "Observed Viral Species",
      sec.axis = sec_axis(~./10, name = "Shannon Diversity")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y.left = element_text(color = "#56B4E9"),
      axis.title.y.right = element_text(color = "#D55E00")
    ) +
    labs(x = "Sample",
         title = "Viral Richness by Sample",
         subtitle = "Bars = Observed Species, Points = Shannon Diversity, Numbers = Viral Contigs")
  
  # Create heatmap of viral taxa
  taxa_heatmap <- taxa_data %>%
    select(name, direct_reads, sample) %>%
    pivot_wider(names_from = name, values_from = direct_reads, values_fill = 0) %>%
    column_to_rownames("sample") %>%
    as.matrix() %>%
    log1p() %>%
    t() %>%
    pheatmap::pheatmap(
      main = "Viral Taxa Abundance Heatmap",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,
      show_colnames = TRUE,
      color = viridis(100),
      border_color = NA,
      fontsize_col = 8,
      silent = TRUE
    )
  
  # Save plots
  ggsave("results/viral_richness/richness_barplot.pdf", richness_plot, width = 10, height = 6)
  pdf("results/viral_richness/viral_heatmap.pdf", width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(taxa_heatmap$gtable)
  dev.off()
  
  # Combine plots for final output
  combined_plot <- gridExtra::grid.arrange(
    richness_plot, 
    grid::rasterGrob(taxa_heatmap$gtable),
    ncol = 1,
    heights = c(1, 2)
  )
  
  return(combined_plot)
}

# Main execution
main <- function() {
  # Get input file paths from Snakemake
  kraken_reports <- snakemake@input[["kraken_reports"]]
  viral_contigs <- snakemake@input[["viral_contigs"]]
  
  # Output files
  richness_report <- snakemake@output[["richness_report"]]
  richness_plot <- snakemake@output[["richness_plot"]]
  
  # Create output directory if it doesn't exist
  dir.create(dirname(richness_report), showWarnings = FALSE, recursive = TRUE)
  
  # Parse Kraken2 reports
  cat("Parsing Kraken2 reports...\n")
  taxa_data <- parse_kraken_reports(kraken_reports)
  
  # Calculate contig statistics
  cat("Calculating viral contig statistics...\n")
  contig_stats <- calculate_contig_stats(viral_contigs)
  
  # Calculate richness indices
  cat("Calculating richness indices...\n")
  richness_data <- calculate_richness(taxa_data)
  
  # Combine richness data with contig stats
  combined_data <- richness_data %>%
    left_join(contig_stats, by = "sample")
  
  # Write report to file
  write.table(combined_data, richness_report, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Create and save plots
  cat("Creating richness plots...\n")
  final_plot <- plot_richness(richness_data, contig_stats)
  
  # Save final plot
  pdf(richness_plot, width = 10, height = 14)
  grid::grid.draw(final_plot)
  dev.off()
  
  cat("Completed viral richness analysis!\n")
}

# Run main function
main()


