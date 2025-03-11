# Snakefile for QIIME2 16S rRNA V3-V4 Metagenomics Pipeline
# Author: Claude
# Date: 2025-03-09

import os
import pandas as pd
from snakemake.utils import min_version

# Set minimum Snakemake version
min_version("7.0")

# Configuration
configfile: "config.yaml"

# Load manifest file
manifest = pd.read_csv(config["manifest_file"], sep=",")
samples = manifest["sample-id"].tolist()

# Output directories
RESULTS_DIR = config["output_dir"]
QC_DIR = os.path.join(RESULTS_DIR, "1_qc")
TAXONOMY_DIR = os.path.join(RESULTS_DIR, "2_taxonomy")
DIVERSITY_DIR = os.path.join(RESULTS_DIR, "3_diversity")
DIFF_ABUNDANCE_DIR = os.path.join(RESULTS_DIR, "4_differential_abundance")

# Define final outputs to ensure all rules are run
rule all:
    input:
        # QC outputs
        os.path.join(QC_DIR, "dada2", "rep-seqs.qza"),
        os.path.join(QC_DIR, "dada2", "table.qza"),
        os.path.join(QC_DIR, "dada2", "stats.qza"),
        os.path.join(QC_DIR, "dada2", "rep-seqs.qzv"),
        os.path.join(QC_DIR, "dada2", "table.qzv"),
        os.path.join(QC_DIR, "dada2", "stats.qzv"),
        
        # Alternate QC method (Deblur)
        os.path.join(QC_DIR, "deblur", "rep-seqs.qza"),
        os.path.join(QC_DIR, "deblur", "table.qza"),
        os.path.join(QC_DIR, "deblur", "stats.qza"),
        
        # Taxonomy outputs
        os.path.join(TAXONOMY_DIR, "taxonomy.qza"),
        os.path.join(TAXONOMY_DIR, "taxonomy.qzv"),
        os.path.join(TAXONOMY_DIR, "taxa-bar-plots.qzv"),
        
        # Diversity outputs
        os.path.join(DIVERSITY_DIR, "core-metrics-results"),
        os.path.join(DIVERSITY_DIR, "alpha-group-significance.qzv"),
        os.path.join(DIVERSITY_DIR, "beta-group-significance.qzv"),
        
        # Differential abundance outputs
        os.path.join(DIFF_ABUNDANCE_DIR, "ancom_results.qzv"),
        os.path.join(DIFF_ABUNDANCE_DIR, "deseq2_results.qzv"),
        os.path.join(DIFF_ABUNDANCE_DIR, "lefse_results.tsv")

# Import data from manifest file
rule import_data:
    input:
        manifest = config["manifest_file"]
    output:
        demux = os.path.join(QC_DIR, "demux.qza")
    params:
        type = config["input_type"]
    shell:
        """
        qiime tools import \
          --type {params.type} \
          --input-path {input.manifest} \
          --output-path {output.demux} \
          --input-format PairedEndFastqManifestPhred33V2
        """

# Visualize demultiplexed sequences
rule visualize_demux:
    input:
        demux = os.path.join(QC_DIR, "demux.qza")
    output:
        demux_viz = os.path.join(QC_DIR, "demux.qzv")
    shell:
        """
        qiime demux summarize \
          --i-data {input.demux} \
          --o-visualization {output.demux_viz}
        """

# Download reference databases (SILVA and Greengenes)
rule download_references:
    output:
        silva_seqs = os.path.join(TAXONOMY_DIR, "silva-138-99-seqs.qza"),
        silva_tax = os.path.join(TAXONOMY_DIR, "silva-138-99-tax.qza"),
        gg_seqs = os.path.join(TAXONOMY_DIR, "gg-13-8-99-515-806-nb-classifier.qza")
    shell:
        """
        # Download SILVA reference sequences
        wget -O silva-138-99-seqs.qza https://data.qiime2.org/2023.5/common/silva-138-99-seqs.qza
        wget -O silva-138-99-tax.qza https://data.qiime2.org/2023.5/common/silva-138-99-tax.qza
        
        # Download pre-trained Greengenes classifier for V4 region (515F-806R)
        wget -O gg-13-8-99-515-806-nb-classifier.qza https://data.qiime2.org/2023.5/common/gg-13-8-99-515-806-nb-classifier.qza
        
        # Move to destination
        mv silva-138-99-seqs.qza {output.silva_seqs}
        mv silva-138-99-tax.qza {output.silva_tax}
        mv gg-13-8-99-515-806-nb-classifier.qza {output.gg_seqs}
        """

# Extract V3-V4 regions from SILVA for classifier training
rule extract_v3v4_regions:
    input:
        silva_seqs = os.path.join(TAXONOMY_DIR, "silva-138-99-seqs.qza"),
        silva_tax = os.path.join(TAXONOMY_DIR, "silva-138-99-tax.qza")
    output:
        silva_v3v4_seqs = os.path.join(TAXONOMY_DIR, "silva-138-99-v3v4-seqs.qza")
    params:
        f_primer = config["v3v4_forward_primer"],
        r_primer = config["v3v4_reverse_primer"]
    shell:
        """
        # Extract specific V3-V4 region from reference database
        qiime feature-classifier extract-reads \
          --i-sequences {input.silva_seqs} \
          --p-f-primer {params.f_primer} \
          --p-r-primer {params.r_primer} \
          --p-trunc-len 0 \
          --o-reads {output.silva_v3v4_seqs}
        """

# Train Naive Bayes classifier for V3-V4 using SILVA
rule train_v3v4_classifier:
    input:
        silva_v3v4_seqs = os.path.join(TAXONOMY_DIR, "silva-138-99-v3v4-seqs.qza"),
        silva_tax = os.path.join(TAXONOMY_DIR, "silva-138-99-tax.qza")
    output:
        silva_v3v4_classifier = os.path.join(TAXONOMY_DIR, "silva-138-99-v3v4-classifier.qza")
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes \
          --i-reference-reads {input.silva_v3v4_seqs} \
          --i-reference-taxonomy {input.silva_tax} \
          --o-classifier {output.silva_v3v4_classifier}
        """

# Quality Control using DADA2
rule qc_dada2:
    input:
        demux = os.path.join(QC_DIR, "demux.qza")
    output:
        table = os.path.join(QC_DIR, "dada2", "table.qza"),
        rep_seqs = os.path.join(QC_DIR, "dada2", "rep-seqs.qza"),
        stats = os.path.join(QC_DIR, "dada2", "stats.qza")
    params:
        trim_left_f = config["dada2_trim_left_f"],
        trim_left_r = config["dada2_trim_left_r"],
        trunc_len_f = config["dada2_trunc_len_f"],
        trunc_len_r = config["dada2_trunc_len_r"],
        n_threads = config["threads"]
    shell:
        """
        qiime dada2 denoise-paired \
          --i-demultiplexed-seqs {input.demux} \
          --p-trim-left-f {params.trim_left_f} \
          --p-trim-left-r {params.trim_left_r} \
          --p-trunc-len-f {params.trunc_len_f} \
          --p-trunc-len-r {params.trunc_len_r} \
          --p-n-threads {params.n_threads} \
          --o-table {output.table} \
          --o-representative-sequences {output.rep_seqs} \
          --o-denoising-stats {output.stats}
        """

# Visualize DADA2 outputs
rule visualize_dada2:
    input:
        table = os.path.join(QC_DIR, "dada2", "table.qza"),
        rep_seqs = os.path.join(QC_DIR, "dada2", "rep-seqs.qza"),
        stats = os.path.join(QC_DIR, "dada2", "stats.qza")
    output:
        table_viz = os.path.join(QC_DIR, "dada2", "table.qzv"),
        rep_seqs_viz = os.path.join(QC_DIR, "dada2", "rep-seqs.qzv"),
        stats_viz = os.path.join(QC_DIR, "dada2", "stats.qzv")
    shell:
        """
        qiime feature-table summarize \
          --i-table {input.table} \
          --o-visualization {output.table_viz}
          
        qiime feature-table tabulate-seqs \
          --i-data {input.rep_seqs} \
          --o-visualization {output.rep_seqs_viz}
          
        qiime metadata tabulate \
          --m-input-file {input.stats} \
          --o-visualization {output.stats_viz}
        """

# Alternative QC method: Deblur
rule qc_deblur:
    input:
        demux = os.path.join(QC_DIR, "demux.qza")
    output:
        table = os.path.join(QC_DIR, "deblur", "table.qza"),
        rep_seqs = os.path.join(QC_DIR, "deblur", "rep-seqs.qza"),
        stats = os.path.join(QC_DIR, "deblur", "stats.qza")
    params:
        trim_length = config["deblur_trim_length"],
        sample_stats = os.path.join(QC_DIR, "deblur", "sample-stats.qza"),
        n_threads = config["threads"]
    shell:
        """
        qiime quality-filter q-score \
          --i-demux {input.demux} \
          --o-filtered-sequences {QC_DIR}/deblur/filtered.qza \
          --o-filter-stats {params.sample_stats}
          
        qiime deblur denoise-16S \
          --i-demultiplexed-seqs {QC_DIR}/deblur/filtered.qza \
          --p-trim-length {params.trim_length} \
          --p-sample-stats \
          --p-jobs-to-start {params.n_threads} \
          --o-table {output.table} \
          --o-representative-sequences {output.rep_seqs} \
          --o-stats {output.stats}
        """

# Taxonomic classification using the trained classifier
rule taxonomic_classification:
    input:
        rep_seqs = os.path.join(QC_DIR, "dada2", "rep-seqs.qza"),
        classifier = os.path.join(TAXONOMY_DIR, "silva-138-99-v3v4-classifier.qza")
    output:
        taxonomy = os.path.join(TAXONOMY_DIR, "taxonomy.qza")
    shell:
        """
        qiime feature-classifier classify-sklearn \
          --i-classifier {input.classifier} \
          --i-reads {input.rep_seqs} \
          --o-classification {output.taxonomy}
        """

# Visualize taxonomy results
rule visualize_taxonomy:
    input:
        taxonomy = os.path.join(TAXONOMY_DIR, "taxonomy.qza"),
        table = os.path.join(QC_DIR, "dada2", "table.qza")
    output:
        taxonomy_viz = os.path.join(TAXONOMY_DIR, "taxonomy.qzv"),
        taxa_bar_plots = os.path.join(TAXONOMY_DIR, "taxa-bar-plots.qzv")
    params:
        metadata = config["metadata_file"]
    shell:
        """
        qiime metadata tabulate \
          --m-input-file {input.taxonomy} \
          --o-visualization {output.taxonomy_viz}
          
        qiime taxa barplot \
          --i-table {input.table} \
          --i-taxonomy {input.taxonomy} \
          --m-metadata-file {params.metadata} \
          --o-visualization {output.taxa_bar_plots}
        """

# Core diversity analyses
rule diversity_analyses:
    input:
        table = os.path.join(QC_DIR, "dada2", "table.qza"),
        rep_seqs = os.path.join(QC_DIR, "dada2", "rep-seqs.qza"),
        taxonomy = os.path.join(TAXONOMY_DIR, "taxonomy.qza")
    output:
        core_metrics_dir = directory(os.path.join(DIVERSITY_DIR, "core-metrics-results")),
        alpha_viz = os.path.join(DIVERSITY_DIR, "alpha-group-significance.qzv"),
        beta_viz = os.path.join(DIVERSITY_DIR, "beta-group-significance.qzv")
    params:
        sampling_depth = config["diversity_sampling_depth"],
        metadata = config["metadata_file"],
        metric = "faith_pd", # for alpha diversity example
        group_column = config["group_column"],
        beta_metric = "weighted_unifrac" # for beta diversity example
    shell:
        """
        # First filter to remove mitochondria and chloroplasts
        qiime taxa filter-table \
          --i-table {input.table} \
          --i-taxonomy {input.taxonomy} \
          --p-exclude "mitochondria,chloroplast" \
          --o-filtered-table {DIVERSITY_DIR}/filtered-table.qza
          
        # Core diversity metrics
        qiime diversity core-metrics-phylogenetic \
          --i-phylogeny {DIVERSITY_DIR}/rooted-tree.qza \
          --i-table {DIVERSITY_DIR}/filtered-table.qza \
          --p-sampling-depth {params.sampling_depth} \
          --m-metadata-file {params.metadata} \
          --output-dir {output.core_metrics_dir}
          
        # Alpha diversity statistical tests
        qiime diversity alpha-group-significance \
          --i-alpha-diversity {output.core_metrics_dir}/faith_pd_vector.qza \
          --m-metadata-file {params.metadata} \
          --o-visualization {output.alpha_viz}
          
        # Beta diversity statistical tests
        qiime diversity beta-group-significance \
          --i-distance-matrix {output.core_metrics_dir}/weighted_unifrac_distance_matrix.qza \
          --m-metadata-file {params.metadata} \
          --m-metadata-column {params.group_column} \
          --o-visualization {output.beta_viz} \
          --p-pairwise
        """

# Phylogenetic tree (needed for UniFrac)
rule phylogenetic_tree:
    input:
        rep_seqs = os.path.join(QC_DIR, "dada2", "rep-seqs.qza")
    output:
        aligned_seqs = os.path.join(DIVERSITY_DIR, "aligned-rep-seqs.qza"),
        masked_aligned_seqs = os.path.join(DIVERSITY_DIR, "masked-aligned-rep-seqs.qza"),
        unrooted_tree = os.path.join(DIVERSITY_DIR, "unrooted-tree.qza"),
        rooted_tree = os.path.join(DIVERSITY_DIR, "rooted-tree.qza")
    shell:
        """
        # Alignment
        qiime alignment mafft \
          --i-sequences {input.rep_seqs} \
          --o-alignment {output.aligned_seqs}
          
        # Mask highly variable positions
        qiime alignment mask \
          --i-alignment {output.aligned_seqs} \
          --o-masked-alignment {output.masked_aligned_seqs}
          
        # Build unrooted tree
        qiime phylogeny fasttree \
          --i-alignment {output.masked_aligned_seqs} \
          --o-tree {output.unrooted_tree}
          
        # Root the tree
        qiime phylogeny midpoint-root \
          --i-tree {output.unrooted_tree} \
          --o-rooted-tree {output.rooted_tree}
        """

# Differential abundance analysis using ANCOM
rule ancom_analysis:
    input:
        table = os.path.join(QC_DIR, "dada2", "table.qza"),
        taxonomy = os.path.join(TAXONOMY_DIR, "taxonomy.qza")
    output:
        ancom_results = os.path.join(DIFF_ABUNDANCE_DIR, "ancom_results.qzv")
    params:
        metadata = config["metadata_file"],
        group_column = config["group_column"],
        taxa_level = config["ancom_taxa_level"]
    shell:
        """
        # Collapse features at specified taxonomic level
        qiime taxa collapse \
          --i-table {input.table} \
          --i-taxonomy {input.taxonomy} \
          --p-level {params.taxa_level} \
          --o-collapsed-table {DIFF_ABUNDANCE_DIR}/collapsed-table.qza
        
        # Convert to composition
        qiime composition add-pseudocount \
          --i-table {DIFF_ABUNDANCE_DIR}/collapsed-table.qza \
          --o-composition-table {DIFF_ABUNDANCE_DIR}/comp-table.qza
        
        # ANCOM analysis
        qiime composition ancom \
          --i-table {DIFF_ABUNDANCE_DIR}/comp-table.qza \
          --m-metadata-file {params.metadata} \
          --m-metadata-column {params.group_column} \
          --o-visualization {output.ancom_results}
        """

# DESeq2 analysis
rule deseq2_analysis:
    input:
        table = os.path.join(QC_DIR, "dada2", "table.qza")
    output:
        deseq2_results = os.path.join(DIFF_ABUNDANCE_DIR, "deseq2_results.qzv")
    params:
        metadata = config["metadata_file"],
        group_column = config["group_column"],
        control_group = config["deseq2_control_group"]
    shell:
        """
        # Convert to biom format for QIIME2 DESeq2 plugin
        qiime tools export \
          --input-path {input.table} \
          --output-path {DIFF_ABUNDANCE_DIR}/exported_biom
        
        # Run DESeq2 analysis
        Rscript {workflow.basedir}/scripts/run_deseq2.R \
          {DIFF_ABUNDANCE_DIR}/exported_biom/feature-table.biom \
          {params.metadata} \
          {params.group_column} \
          {params.control_group} \
          {DIFF_ABUNDANCE_DIR}/deseq2_results.tsv
        
        # Import results back to QIIME2 visualization
        qiime metadata tabulate \
          --m-input-file {DIFF_ABUNDANCE_DIR}/deseq2_results.tsv \
          --o-visualization {output.deseq2_results}
        """

# LEfSe analysis
rule lefse_analysis:
    input:
        table = os.path.join(QC_DIR, "dada2", "table.qza"),
        taxonomy = os.path.join(TAXONOMY_DIR, "taxonomy.qza")
    output:
        lefse_results = os.path.join(DIFF_ABUNDANCE_DIR, "lefse_results.tsv")
    params:
        metadata = config["metadata_file"],
        group_column = config["group_column"]
    shell:
        """
        # Prepare data for LEfSe
        python {workflow.basedir}/scripts/prep_for_lefse.py \
          --table {input.table} \
          --taxonomy {input.taxonomy} \
          --metadata {params.metadata} \
          --group-col {params.group_column} \
          --output {DIFF_ABUNDANCE_DIR}/lefse_input.txt
        
        # Run LEfSe
        run_lefse.py {DIFF_ABUNDANCE_DIR}/lefse_input.txt {output.lefse_results}
        """

