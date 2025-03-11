#!/usr/bin/env Rscript

# Script for time series analysis of microbiome data using Phyloseq
# This script is designed to analyze temporal changes in microbiome composition

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(biomformat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(vegan)
  library(RColorBrewer)
  library(nlme)      # For linear mixed effects models
  library(lme4)      # For generalized linear mixed effects models
  library(gridExtra) # For arranging multiple plots
})

# Parse command line arguments
option_list <- list(
  make_option("--biom", type="character", help="Path to BIOM file"),
  make_option("--taxonomy", type="character", help="Path to taxonomy file"),
  make_option("--metadata", type="character", help="Path to metadata file"),
  make_option("--time-col", type="character", help="Column name for time points"),
  make_option("--subject-col", type="character", help="Column name for subject IDs"),
  make_option("--taxa-level", type="character", default="Genus", help="Taxonomic level for analysis"),
  make_option("--output", type="character", help="Output file path for plots")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory if needed
output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load data into phyloseq object
cat("Loading data into phyloseq for time series analysis...\n")

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

# Verify that time and subject columns exist
if(!(opt$time_col %in% colnames(sample_data(ps)))) {
  stop(paste0("Time column '", opt$time_col, "' not found in metadata"))
}

if(!(opt$subject_col %in% colnames(sample_data(ps)))) {
  stop(paste0("Subject column '", opt$subject_col, "' not found in metadata"))
}

# Convert time variable to numeric if possible
if(is.character(sample_data(ps)[[opt$time_col]])) {
  if(all(grepl("^\\d+$", na.omit(sample_data(ps)[[opt$time_col]])))) {
    sample_data(ps)[[opt$time_col]] <- as.numeric(sample_data(ps)[[opt$time_col]])
  } else {
    # Try to convert date strings to days from first date
    tryCatch({
      dates <- as.Date(sample_data(ps)[[opt$time_col]])
      if(!all(is.na(dates))) {
        first_date <- min(dates, na.rm=TRUE)
        sample_data(ps)[[paste0(opt$time_col, "_numeric")]] <- as.numeric(dates - first_date)
        opt$time_col <- paste0(opt$time_col, "_numeric")
        cat("Converted date column to numeric days from first date\n")
      }
    }, error=function(e) {
      warning("Could not convert time column to numeric. Using as is.")
    })
  }
}

# Ensure subject column is a factor
sample_data(ps)[[opt$subject_col]] <- as.factor(sample_data(ps)[[opt$subject_col]])

# Create PDF for all plots
pdf(opt$output, width=12, height=10)

# Normalize to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Aggregate at the specified taxonomic level
if(opt$taxa_level %in% colnames(tax_table(ps_rel))) {
  ps_agg <- tax_glom(ps_rel, taxrank=opt$taxa_level)
} else {
  warning(paste0("Taxonomic level '", opt$taxa_level, "' not found. Using raw OTUs."))
  ps_agg <- ps_rel
}

cat("Analyzing temporal patterns in microbiome data...\n")

# Extract time and subject data
time_data <- sample_data(ps_agg)[[opt$time_col]]
subject_data <- sample_data(ps_agg)[[opt$subject_col]]

# Check if we have sufficient time points
unique_times <- unique(time_data)
if(length(unique_times) < 2) {
  stop("Need at least 2 time points for time series analysis")
}

# 1. Alpha diversity over time
alpha_div <- estimate_richness(ps_agg, measures=c("Observed", "Shannon", "Simpson"))
alpha_div$Time <- time_data
alpha_div$Subject <- subject_data

# Melt alpha diversity measures for plotting
alpha_long <- melt(alpha_div, id.vars=c("Time", "Subject"), 
                   variable.name="Measure", value.name="Value")

# Plot alpha diversity over time
p_alpha_time <- ggplot(alpha_long, aes(x=Time, y=Value, color=Subject, group=Subject)) +
  geom_line() +
  geom_point(size=2) +
  facet_wrap(~Measure, scales="free_y") +
  theme_bw() +
  labs(title="Alpha Diversity Changes Over Time", x="Time", y="Alpha Diversity") +
  theme(legend.position="bottom")

print(p_alpha_time)

# 2. Beta diversity changes over time
# Calculate pairwise distances
dist_methods <- c("bray", "jaccard")

for(method in dist_methods) {
  # Calculate distances
  distances <- phyloseq::distance(ps_agg, method=method)
  
  # Create distance matrix
  dist_mat <- as.matrix(distances)
  
  # Create data frame for plotting
  dist_df <- data.frame()
  
  # For each subject, get consecutive time points
  for(subj in unique(subject_data)) {
    subj_idx <- which(subject_data == subj)
    if(length(subj_idx) > 1) {
      subj_times <- time_data[subj_idx]
      subj_samples <- sample_names(ps_agg)[subj_idx]
      
      # Order by time
      time_order <- order(subj_times)
      subj_times <- subj_times[time_order]
      subj_samples <- subj_samples[time_order]
      
      # Get consecutive distances
      for(i in 1:(length(subj_samples)-1)) {
        time_diff <- subj_times[i+1] - subj_times[i]
        sample_dist <- dist_mat[subj_samples[i], subj_samples[i+1]]
        
        dist_df <- rbind(dist_df, data.frame(
          Subject = subj,
          TimePoint1 = subj_times[i],
          TimePoint2 = subj_times[i+1],
          TimeDiff = time_diff,
          Distance = sample_dist
        ))
      }
    }
  }
  
  if(nrow(dist_df) > 0) {
    # Plot distance vs time difference
    p_dist <- ggplot(dist_df, aes(x=TimeDiff, y=Distance, color=Subject)) +
      geom_point(size=3) +
      geom_smooth(method="lm", se=TRUE, color="black") +
      theme_bw() +
      labs(title=paste("Community Dissimilarity Over Time -", method),
           x="Time Difference", y=paste(method, "Distance")) +
      theme(legend.position="bottom")
    
    print(p_dist)
    
    # Test if distance increases with time (across all subjects)
    lm_result <- lm(Distance ~ TimeDiff, data=dist_df)
    lm_summary <- summary(lm_result)
    
    # Create text to show regression results
    plot_text <- data.frame(
      x = min(dist_df$TimeDiff) + 0.1 * (max(dist_df$TimeDiff) - min(dist_df$TimeDiff)),
      y = max(dist_df$Distance) - 0.1 * (max(dist_df$Distance) - min(dist_df$Distance)),
      label = paste("Slope =", round(lm_result$coefficients[2], 4),
                    "\np-value =", round(summary(lm_result)$coefficients[2,4], 4),
                    "\nR² =", round(lm_summary$r.squared, 4))
    )
    
    p_dist_stats <- p_dist + 
      geom_text(data=plot_text, aes(x=x, y=y, label=label), 
                inherit.aes=FALSE, hjust=0, size=3.5)
    
    print(p_dist_stats)
    
    # Mixed-effects model to account for subject variation
    tryCatch({
      me_model <- lme(Distance ~ TimeDiff, random = ~1|Subject, data=dist_df)
      me_summary <- summary(me_model)
      
      # Create text to show mixed effects model results
      me_text <- data.frame(
        x = min(dist_df$TimeDiff) + 0.1 * (max(dist_df$TimeDiff) - min(dist_df$TimeDiff)),
        y = max(dist_df$Distance) - 0.2 * (max(dist_df$Distance) - min(dist_df$Distance)),
        label = paste("Mixed Effects Model:",
                      "\nFixed effect (TimeDiff) =", round(fixef(me_model)[2], 4),
                      "\np-value =", round(summary(me_model)$tTable[2,5], 4))
      )
      
      p_dist_stats <- p_dist_stats + 
        geom_text(data=me_text, aes(x=x, y=y, label=label), 
                  inherit.aes=FALSE, hjust=0, size=3.5, color="blue")
      
      print(p_dist_stats)
    }, error=function(e) {
      warning("Could not fit mixed effects model, continuing with fixed effects only")
    })
  }
}

# 3. Taxa abundance changes over time
# Get top taxa for plotting
top_taxa_count <- 20
top_taxa <- names(sort(taxa_sums(ps_agg), decreasing=TRUE)[1:min(top_taxa_count, ntaxa(ps_agg))])
ps_top <- prune_taxa(top_taxa, ps_agg)

# Melt data for plotting
ps_melt <- psmelt(ps_top)
ps_melt$Time <- ps_melt[[opt$time_col]]
ps_melt$Subject <- ps_melt[[opt$subject_col]]

# Create taxa names for plotting
ps_melt$TaxaLabel <- ps_melt[[opt$taxa_level]]
if("Species" %in% colnames(ps_melt) && opt$taxa_level != "Species") {
  ps_melt$TaxaLabel <- paste(ps_melt$TaxaLabel, ps_melt$Species)
}

# Remove NAs and clean labels
ps_melt$TaxaLabel <- gsub("NA NA", "Unclassified", ps_melt$TaxaLabel)
ps_melt$TaxaLabel <- gsub("NA$", "", ps_melt$TaxaLabel)

# Plot top taxa over time for each subject
for(subj in unique(ps_melt$Subject)) {
  subj_data <- ps_melt[ps_melt$Subject == subj, ]
  
  # Skip subjects with only one time point
  if(length(unique(subj_data$Time)) <= 1) {
    next
  }
  
  # Plot taxa abundances over time
  p_taxa <- ggplot(subj_data, aes(x=Time, y=Abundance, color=TaxaLabel, group=TaxaLabel)) +
    geom_line(size=1) +
    geom_point(size=2) +
    theme_bw() +
    scale_color_brewer(palette="Paired") +
    labs(title=paste("Taxa Abundance Over Time -", subj),
         x="Time", y="Relative Abundance", color=opt$taxa_level) +
    theme(legend.position="right",
          legend.text=element_text(size=8),
          legend.key.size=unit(0.8, "lines"))
  
  print(p_taxa)
  
  # Stacked area plot
  p_area <- ggplot(subj_data, aes(x=Time, y=Abundance, fill=TaxaLabel)) +
    geom_area(position="stack") +
    theme_bw() +
    scale_fill_brewer(palette="Paired") +
    labs(title=paste("Taxa Composition Over Time -", subj),
         x="Time", y="Relative Abundance", fill=opt$taxa_level) +
    theme(legend.position="right",
          legend.text=element_text(size=8),
          legend.key.size=unit(0.8, "lines"))
  
  print(p_area)
}

# 4. Identify taxa with significant temporal trends
# Only do this if there are at least 3 time points
if(length(unique(time_data)) >= 3) {
  cat("Identifying taxa with significant temporal trends...\n")
  
  # Get all taxa abundances
  ps_all_melt <- psmelt(ps_agg)
  ps_all_melt$Time <- ps_all_melt[[opt$time_col]]
  ps_all_melt$Subject <- ps_all_melt[[opt$subject_col]]
  
  # Create taxa names
  ps_all_melt$TaxaLabel <- ps_all_melt[[opt$taxa_level]]
  ps_all_melt$TaxaLabel <- gsub("NA$", "", ps_all_melt$TaxaLabel)
  
  # Iterate through all taxa and test for temporal trends
  taxa_list <- unique(ps_all_melt$TaxaLabel)
  results_df <- data.frame()
  
  for(taxon in taxa_list) {
    taxon_data <- ps_all_melt[ps_all_melt$TaxaLabel == taxon, ]
    
    # Skip taxa with very low abundance
    if(all(taxon_data$Abundance < 0.01)) {
      next
    }
    
    # Linear model with time as predictor
    tryCatch({
      # Simple linear model (ignoring subject)
      lm_simple <- lm(Abundance ~ Time, data=taxon_data)
      
      # Mixed effects model (accounting for subject)
      mixed_model <- lmer(Abundance ~ Time + (1|Subject), data=taxon_data)
      mixed_summary <- summary(mixed_model)
      
      # Extract statistics
      time_effect <- fixef(mixed_model)[2]
      se <- sqrt(diag(vcov(mixed_model)))[2]
      t_value <- time_effect / se
      
      # Calculate p-value
      # Approximation using degrees of freedom
      df <- length(taxon_data$Abundance) - length(fixef(mixed_model)) - 1
      p_value <- 2 * pt(abs(t_value), df, lower.tail=FALSE)
      
      # Store results
      results_df <- rbind(results_df, data.frame(
        Taxon = taxon,
        TimeEffect = time_effect,
        SE = se,
        Tvalue = t_value,
        Pvalue = p_value
      ))
    }, error=function(e) {
      warning(paste("Could not fit model for taxon:", taxon))
    })
  }
  
  # Adjust p-values for multiple testing
  if(nrow(results_df) > 0) {
    results_df$Padj <- p.adjust(results_df$Pvalue, method="BH")
    
    # Sort by significance
    results_df <- results_df[order(results_df$Padj), ]
    
    # Find significant taxa
    sig_taxa <- results_df$Taxon[results_df$Padj < 0.05]
    
    # Plot significant taxa
    if(length(sig_taxa) > 0) {
      sig_data <- ps_all_melt[ps_all_melt$TaxaLabel %in% sig_taxa, ]
      
      # Create a plot for each significant taxon
      for(taxon in sig_taxa) {
        taxon_data <- sig_data[sig_data$TaxaLabel == taxon, ]
        
        time_effect <- results_df$TimeEffect[results_df$Taxon == taxon]
        p_value <- results_df$Padj[results_df$Taxon == taxon]
        
        direction <- ifelse(time_effect > 0, "Increasing", "Decreasing")
        
        p_sig <- ggplot(taxon_data, aes(x=Time, y=Abundance, color=Subject, group=Subject)) +
          geom_line() +
          geom_point(size=2) +
          geom_smooth(method="lm", color="black", se=TRUE, aes(group=1)) +
          theme_bw() +
          labs(title=paste(taxon, "-", direction, "over time"),
               subtitle=paste("Time effect =", round(time_effect, 4), ", Adj. p-value =", 
                              round(p_value, 4)),
               x="Time", y="Relative Abundance") +
          theme(legend.position="bottom")
        
        print(p_sig)
      }
      
      # Heatmap of significant taxa
      if(length(sig_taxa) > 1) {
        # Create taxa x subject matrix
        taxa_matrix <- matrix(NA, nrow=length(sig_taxa), ncol=length(unique(time_data)))
        rownames(taxa_matrix) <- sig_taxa
        colnames(taxa_matrix) <- sort(unique(time_data))
        
        # Fill in mean abundances at each time point
        for(i in 1:nrow(taxa_matrix)) {
          taxon <- rownames(taxa_matrix)[i]
          taxon_data <- sig_data[sig_data$TaxaLabel == taxon, ]
          
          for(j in 1:ncol(taxa_matrix)) {
            time_point <- as.numeric(colnames(taxa_matrix)[j])
            time_data <- taxon_data[taxon_data$Time == time_point, ]
            
            if(nrow(time_data) > 0) {
              taxa_matrix[i, j] <- mean(time_data$Abundance)
            }
          }
        }
        
        # Scale rows for better visualization
        taxa_matrix_scaled <- t(scale(t(taxa_matrix)))
        
        # Plot heatmap
        # Convert to data frame for ggplot
        taxa_df <- as.data.frame(taxa_matrix)
        taxa_df$Taxon <- rownames(taxa_df)
        taxa_long <- melt(taxa_df, id.vars="Taxon", variable.name="TimePoint", value.name="Abundance")
        taxa_long$TimePoint <- as.numeric(as.character(taxa_long$TimePoint))
        
        # Get trend direction for each taxon
        taxa_long$Direction <- "Stable"
        for(taxon in unique(taxa_long$Taxon)) {
          effect <- results_df$TimeEffect[results_df$Taxon == taxon]
          if(!is.na(effect)) {
            direction <- ifelse(effect > 0, "Increasing", "Decreasing")
            taxa_long$Direction[taxa_long$Taxon == taxon] <- direction
          }
        }
        
        # Plot heatmap
        p_heatmap <- ggplot(taxa_long, aes(x=TimePoint, y=Taxon, fill=Abundance)) +
          geom_tile() +
          scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
          theme_bw() +
          theme(axis.text.y=element_text(size=8)) +
          labs(title="Significant Taxa Abundance Changes Over Time",
               x="Time Point", y="Taxon", fill="Z-score") +
          facet_wrap(~Direction, scales="free_y")
        
        print(p_heatmap)
      }
    } else {
      grid.newpage()
      grid.text("No taxa with significant temporal trends found (adjusted p < 0.05)",
                x=0.5, y=0.5, gp=gpar(cex=1.2))
    }
    
    # Save results to file
    write.csv(results_df, file=paste0(tools::file_path_sans_ext(opt$output), "_taxa_time_effects.csv"),
             row.names=FALSE)
  } else {
    grid.newpage()
    grid.text("No taxa with sufficient data for temporal trend analysis",
              x=0.5, y=0.5, gp=gpar(cex=1.2))
  }
}

# 5. Longitudinal sample stability analysis
# Calculate the stability of each subject's microbiome over time

# For each subject, calculate mean distance between consecutive time points
stability_df <- data.frame()

for(subj in unique(subject_data)) {
  subj_idx <- which(subject_data == subj)
  
  if(length(subj_idx) > 1) {
    subj_times <- time_data[subj_idx]
    subj_samples <- sample_names(ps_agg)[subj_idx]
    
    # Order by time
    time_order <- order(subj_times)
    subj_times <- subj_times[time_order]
    subj_samples <- subj_samples[time_order]
    
    # Calculate distances between consecutive samples
    consecutive_dists <- c()
    
    for(dist_method in dist_methods) {
      dist_mat <- as.matrix(phyloseq::distance(ps_agg, method=dist_method))
      
      subj_dists <- c()
      for(i in 1:(length(subj_samples)-1)) {
        subj_dists <- c(subj_dists, dist_mat[subj_samples[i], subj_samples[i+1]])
      }
      
      stability_df <- rbind(stability_df, data.frame(
        Subject = subj,
        Method = dist_method,
        MeanDistance = mean(subj_dists),
        MedianDistance = median(subj_dists),
        SDDistance = sd(subj_dists),
        MinDistance = min(subj_dists),
        MaxDistance = max(subj_dists),
        NumTimePoints = length(subj_samples)
      ))
    }
  }
}

if(nrow(stability_df) > 0) {
  # Plot stability by subject
  for(dist_method in unique(stability_df$Method)) {
    method_data <- stability_df[stability_df$Method == dist_method, ]
    
    # Sort by mean distance
    method_data <- method_data[order(method_data$MeanDistance), ]
    method_data$Subject <- factor(method_data$Subject, levels=method_data$Subject)
    
    p_stability <- ggplot(method_data, aes(x=Subject, y=MeanDistance)) +
      geom_bar(stat="identity", fill="steelblue") +
      geom_errorbar(aes(ymin=MeanDistance-SDDistance, ymax=MeanDistance+SDDistance), width=0.2) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      labs(title=paste("Microbiome Stability by Subject -", dist_method),
           subtitle="Lower values indicate more stable microbiomes over time",
           x="Subject", y=paste("Mean", dist_method, "Distance"))
    
    print(p_stability)
  }
  
  # Save stability results
  write.csv(stability_df, file=paste0(tools::file_path_sans_ext(opt$output), "_subject_stability.csv"),
           row.names=FALSE)
}

# 6. Principal coordinates analysis with time trajectories
for(dist_method in dist_methods) {
  # Run ordination
  pcoa <- ordinate(ps_agg, method="PCoA", distance=dist_method)
  
  # Extract PC scores
  pc_data <- as.data.frame(pcoa$vectors)
  pc_data$Sample <- rownames(pc_data)
  
  # Add metadata
  meta_df <- as.data.frame(sample_data(ps_agg))
  meta_df$Sample <- rownames(meta_df)
  pc_data <- merge(pc_data, meta_df, by="Sample")
  
  # Add time and subject info
  pc_data$Time <- pc_data[[opt$time_col]]
  pc_data$Subject <- pc_data[[opt$subject_col]]
  
  # Calculate percent variation explained
  var_explained <- pcoa$values$Relative_eig * 100
  
  # Plot PCoA with time as color
  p_pcoa <- ggplot(pc_data, aes(x=Axis.1, y=Axis.2, color=Time)) +
    geom_point(size=3) +
    theme_bw() +
    labs(title=paste("PCoA of", dist_method, "Distances with Time Gradient"),
         x=paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y=paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    scale_color_gradient(low="blue", high="red")
  
  print(p_pcoa)
  
  # Plot PCoA with subject trajectories
  p_traj <- ggplot(pc_data, aes(x=Axis.1, y=Axis.2, color=Subject)) +
    geom_point(size=3) +
    geom_path(aes(group=Subject), alpha=0.7) +
    theme_bw() +
    labs(title=paste("PCoA of", dist_method, "Distances with Subject Trajectories"),
         x=paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y=paste0("PC2 (", round(var_explained[2], 1), "%)"))
  
  print(p_traj)
  
  # Add time point labels to show direction
  pc_data$TimeLabel <- as.character(round(pc_data$Time, 1))
  
  p_traj_labels <- p_traj +
    geom_text(aes(label=TimeLabel), hjust=0, vjust=0, size=3, nudge_x=0.01, nudge_y=0.01)
  
  print(p_traj_labels)
}

# Close PDF
dev.off()

# Create summary text file
summary_file <- paste0(tools::file_path_sans_ext(opt$output), "_summary.txt")
sink(summary_file)
cat("Time Series Analysis Summary\n")
cat("==========================\n\n")
cat("Data information:\n")
cat("- Number of samples:", nsamples(ps), "\n")
cat("- Number of taxa:", ntaxa(ps), "\n")
cat("- Number of subjects:", length(unique(subject_data)), "\n")
cat("- Number of time points:", length(unique(time_data)), "\n")
cat("- Time range:", min(time_data), "to", max(time_data), "\n\n")

if(exists("results_df") && nrow(results_df) > 0) {
  cat("Significant temporal trends:\n")
  sig_results <- results_df[results_df$Padj < 0.05, ]
  if(nrow(sig_results) > 0) {
    for(i in 1:nrow(sig_results)) {
      direction <- ifelse(sig_results$TimeEffect[i] > 0, "increasing", "decreasing")
      cat("- ", sig_results$Taxon[i], ": ", direction, 
          " (effect = ", round(sig_results$TimeEffect[i], 4),
          ", p-adj = ", round(sig_results$Padj[i], 4), ")\n", sep="")
    }
  } else {
    cat("- No significant temporal trends detected\n")
  }
}
sink()

cat("Time series analysis complete! Results saved to:", opt$output, "\n")

