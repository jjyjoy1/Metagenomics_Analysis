This Snakemake extension(Snakefile_cross_sample.smk) adds the following components to the single sample snakemake(Snakefile) pipeline:

#1. Cross-Sample Analysis Structure:

Establishes a checkpoint to ensure all individual samples are processed
Creates cross-sample directory structure
Defines a final rule all_cross_sample to run all cross-sample analyses


#2. Pan-genome Matrix Generation:

collect_genes: Combines gene predictions from all samples
cluster_genes: Clusters genes using CD-HIT to identify orthologs
generate_pangenome_matrix: Creates a presence/absence matrix and identifies core/accessory/unique genes


#3. Functional Profile Matrix Generation:

eggnog_annotation: Adds functional annotation using eggNOG-mapper
extract_ko_annotations and combine_ko_counts: Creates KO abundance matrix
extract_pathway_annotations and combine_pathway_counts: Creates pathway abundance matrix


#4. Taxonomic Abundance Matrix Generation:

process_cat_taxonomy: Processes CAT taxonomy results with coverage info
combine_taxonomy: Creates abundance matrices at phylum, genus, and species levels
process_bin_taxonomy and combine_bin_taxonomy: Creates bin-level taxonomy matrix


#5. Integration and Visualization:

integrate_with_metadata: Combines all matrices with sample metadata
visualize_matrices: Creates heatmaps, PCA plots, and summary visualizations
