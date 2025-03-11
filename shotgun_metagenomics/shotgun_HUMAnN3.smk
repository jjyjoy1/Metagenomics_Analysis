# Functional Profiling Pipeline
# Using HUMAnN3 for metabolic pathway analysis and functional profiling

# Configuration
configfile: "functional_config.yaml"

# Import sample information from manifest file
import pandas as pd
import os

# Read manifest file
manifest = pd.read_csv(config["manifest_file"])
SAMPLES = manifest["sample_id"].tolist()
READ_DIRS = manifest["read_dir"].tolist()

# Get file paths for paired-end reads
def get_read_paths(wildcards):
    idx = SAMPLES.index(wildcards.sample)
    read_dir = READ_DIRS[idx]
    return {
        "r1": os.path.join(read_dir, f"{wildcards.sample}_R1.fastq.gz"),
        "r2": os.path.join(read_dir, f"{wildcards.sample}_R2.fastq.gz")
    }

# Define final output files
rule all:
    input:
        # HUMAnN3 outputs
        expand("results/humann3/{sample}/{sample}_genefamilies.tsv", sample=SAMPLES),
        expand("results/humann3/{sample}/{sample}_pathabundance.tsv", sample=SAMPLES),
        expand("results/humann3/{sample}/{sample}_pathcoverage.tsv", sample=SAMPLES),
        
        # Merged tables
        "results/humann3/merged/genefamilies.tsv",
        "results/humann3/merged/pathabundance.tsv",
        "results/humann3/merged/pathcoverage.tsv",
        
        # Regrouped tables
        "results/humann3/regrouped/genefamilies_ko.tsv",
        "results/humann3/regrouped/genefamilies_ec.tsv",
        "results/humann3/regrouped/genefamilies_go.tsv",
        
        # Stratified tables (per species contribution)
        "results/humann3/stratified/pathabundance_stratified.tsv",
        
        # Visualization files
        "results/humann3/visualization/top_pathways_heatmap.png",
        "results/humann3/visualization/functional_pca.png",
        
        # Final report
        "results/humann3/reports/functional_profiling_report.html"

# Concatenate paired-end reads for HUMAnN3
rule concat_reads:
    input:
        unpack(get_read_paths)
    output:
        concat="results/humann3_prep/{sample}_concat.fastq.gz"
    shell:
        """
        mkdir -p results/humann3_prep
        cat {input.r1} {input.r2} > {output.concat}
        """

# Run HUMAnN3 for functional profiling
rule humann3:
    input:
        concat="results/humann3_prep/{sample}_concat.fastq.gz",
        metaphlan_profile="results/metaphlan/{sample}.metaphlan.profile.txt"  # Use MetaPhlAn results from main pipeline if available
    output:
        genefamilies="results/humann3/{sample}/{sample}_genefamilies.tsv",
        pathabundance="results/humann3/{sample}/{sample}_pathabundance.tsv",
        pathcoverage="results/humann3/{sample}/{sample}_pathcoverage.tsv"
    params:
        outdir="results/humann3/{sample}"
    threads: 16
    shell:
        """
        mkdir -p {params.outdir}
        
        # Check if MetaPhlAn profile exists, use it if available for faster processing
        if [ -s {input.metaphlan_profile} ]; then
            humann --input {input.concat} \
                   --output {params.outdir} \
                   --threads {threads} \
                   --taxonomic-profile {input.metaphlan_profile} \
                   --nucleotide-database {config[humann_nucleotide_db]} \
                   --protein-database {config[humann_protein_db]} \
                   --output-basename {wildcards.sample}
        else
            humann --input {input.concat} \
                   --output {params.outdir} \
                   --threads {threads} \
                   --nucleotide-database {config[humann_nucleotide_db]} \
                   --protein-database {config[humann_protein_db]} \
                   --output-basename {wildcards.sample}
        fi
        """

# Normalize HUMAnN3 output to relative abundance
rule normalize_humann:
    input:
        genefamilies="results/humann3/{sample}/{sample}_genefamilies.tsv",
        pathabundance="results/humann3/{sample}/{sample}_pathabundance.tsv"
    output:
        genefamilies_relab="results/humann3/{sample}/{sample}_genefamilies_relab.tsv",
        pathabundance_relab="results/humann3/{sample}/{sample}_pathabundance_relab.tsv"
    shell:
        """
        # Normalize gene families to relative abundance
        humann_renorm_table --input {input.genefamilies} \
                           --output {output.genefamilies_relab} \
                           --units relab
        
        # Normalize pathway abundance to relative abundance
        humann_renorm_table --input {input.pathabundance} \
                           --output {output.pathabundance_relab} \
                           --units relab
        """

# Merge results from all samples
rule merge_humann:
    input:
        genefamilies=expand("results/humann3/{sample}/{sample}_genefamilies_relab.tsv", sample=SAMPLES),
        pathabundance=expand("results/humann3/{sample}/{sample}_pathabundance_relab.tsv", sample=SAMPLES),
        pathcoverage=expand("results/humann3/{sample}/{sample}_pathcoverage.tsv", sample=SAMPLES)
    output:
        genefamilies="results/humann3/merged/genefamilies.tsv",
        pathabundance="results/humann3/merged/pathabundance.tsv",
        pathcoverage="results/humann3/merged/pathcoverage.tsv"
    shell:
        """
        mkdir -p results/humann3/merged
        
        # Merge gene families
        humann_join_tables --input results/humann3/ \
                          --output {output.genefamilies} \
                          --file_name genefamilies_relab
        
        # Merge pathway abundance
        humann_join_tables --input results/humann3/ \
                          --output {output.pathabundance} \
                          --file_name pathabundance_relab
        
        # Merge pathway coverage
        humann_join_tables --input results/humann3/ \
                          --output {output.pathcoverage} \
                          --file_name pathcoverage
        """

# Regroup gene families to different functional categories
rule regroup_humann:
    input:
        genefamilies="results/humann3/merged/genefamilies.tsv"
    output:
        ko="results/humann3/regrouped/genefamilies_ko.tsv",
        ec="results/humann3/regrouped/genefamilies_ec.tsv",
        go="results/humann3/regrouped/genefamilies_go.tsv"
    shell:
        """
        mkdir -p results/humann3/regrouped
        
        # Regroup to KEGG Orthology (KO)
        humann_regroup_table --input {input.genefamilies} \
                             --output {output.ko} \
                             --groups uniref90_ko
        
        # Regroup to Enzyme Commission (EC)
        humann_regroup_table --input {input.genefamilies} \
                             --output {output.ec} \
                             --groups uniref90_level4ec
        
        # Regroup to Gene Ontology (GO)
        humann_regroup_table --input {input.genefamilies} \
                             --output {output.go} \
                             --groups uniref90_go
        """

# Split stratified tables to see species contribution
rule stratify_humann:
    input:
        pathabundance="results/humann3/merged/pathabundance.tsv"
    output:
        stratified="results/humann3/stratified/pathabundance_stratified.tsv"
    shell:
        """
        mkdir -p results/humann3/stratified
        
        # Split stratified table
        humann_split_stratified_table --input {input.pathabundance} \
                                     --output results/humann3/stratified
        
        # Rename the output file
        mv results/humann3/stratified/pathabundance_stratified.tsv {output.stratified}
        """

# Create visualizations
rule visualize_humann:
    input:
        pathabundance="results/humann3/merged/pathabundance.tsv",
        genefamilies_ko="results/humann3/regrouped/genefamilies_ko.tsv"
    output:
        pathway_heatmap="results/humann3/visualization/top_pathways_heatmap.png",
        pca_plot="results/humann3/visualization/functional_pca.png"
    script:
        "scripts/visualize_functional.py"

# Generate HTML report
rule functional_report:
    input:
        pathabundance="results/humann3/merged/pathabundance.tsv",
        ko="results/humann3/regrouped/genefamilies_ko.tsv",
        pathway_heatmap="results/humann3/visualization/top_pathways_heatmap.png",
        pca_plot="results/humann3/visualization/functional_pca.png"
    output:
        report="results/humann3/reports/functional_profiling_report.html"
    script:
        "scripts/functional_report.py"


