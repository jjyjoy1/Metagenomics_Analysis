# Virulence Factor Detection Pipeline
# Identifies virulence factors in metagenomic samples for pathogen characterization

# Configuration
configfile: "virulence_config.yaml"

# Import sample information from manifest file
import pandas as pd
import os

# Read manifest file
manifest = pd.read_csv(config["manifest_file"])
SAMPLES = manifest["sample_id"].tolist()

# Define the final output files
rule all:
    input:
        # Virulence factor detection using ABRicate
        expand("results/virulence/abricate/{sample}.vfdb.tsv", sample=SAMPLES),
        expand("results/virulence/abricate/{sample}.victors.tsv", sample=SAMPLES),
        
        # Virulence factor detection using SRST2 (direct from reads)
        expand("results/virulence/srst2/{sample}.virulence.scores", sample=SAMPLES),
        
        # Protein-based virulence detection using DIAMOND and VFDB
        expand("results/virulence/diamond/{sample}.virulence.tsv", sample=SAMPLES),
        
        # Combined virulence factor report
        "results/virulence/reports/combined_virulence_factors.tsv",
        
        # Virulome visualization
        "results/virulence/visualization/virulome_heatmap.png",
        
        # Final virulence report
        "results/virulence/reports/virulence_profiling_report.html"

# Use assemblies from the main pipeline
rule get_assemblies:
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/virulence/assemblies/{sample}.contigs.fasta"
    shell:
        """
        mkdir -p results/virulence/assemblies
        cp {input} {output}
        """

# ABRicate for virulence factor detection using VFDB
rule abricate_vfdb:
    input:
        contigs="results/virulence/assemblies/{sample}.contigs.fasta"
    output:
        report="results/virulence/abricate/{sample}.vfdb.tsv"
    threads: 8
    shell:
        """
        mkdir -p results/virulence/abricate
        
        # Run ABRicate with VFDB database
        abricate --db vfdb --minid 80 --mincov 80 --threads {threads} {input.contigs} > {output.report}
        """

# ABRicate for virulence factor detection using VICTORS
rule abricate_victors:
    input:
        contigs="results/virulence/assemblies/{sample}.contigs.fasta"
    output:
        report="results/virulence/abricate/{sample}.victors.tsv"
    threads: 8
    shell:
        """
        mkdir -p results/virulence/abricate
        
        # Run ABRicate with VICTORS database
        abricate --db victors --minid 80 --mincov 80 --threads {threads} {input.contigs} > {output.report}
        """

# Summarize ABRicate results
rule abricate_summary:
    input:
        vfdb=expand("results/virulence/abricate/{sample}.vfdb.tsv", sample=SAMPLES),
        victors=expand("results/virulence/abricate/{sample}.victors.tsv", sample=SAMPLES)
    output:
        vfdb_summary="results/virulence/abricate/summary_vfdb.tsv",
        victors_summary="results/virulence/abricate/summary_victors.tsv"
    shell:
        """
        # Summarize VFDB results
        abricate --summary {input.vfdb} > {output.vfdb_summary}
        
        # Summarize VICTORS results
        abricate --summary {input.victors} > {output.victors_summary}
        """

# SRST2 for virulence factor detection directly from reads
rule srst2_virulence:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        scores="results/virulence/srst2/{sample}.virulence.scores",
        report="results/virulence/srst2/{sample}.virulence.gene_report.tsv"
    params:
        prefix="results/virulence/srst2/{sample}.virulence",
        vf_db=config["srst2_virulence_db"]
    threads: 8
    shell:
        """
        mkdir -p results/virulence/srst2
        
        # Run SRST2 for virulence gene detection
        srst2 --input_pe {input.r1} {input.r2} \
              --output {params.prefix} \
              --log \
              --gene_db {params.vf_db} \
              --min_coverage 80 \
              --threads {threads}
        """

# DIAMOND for protein-based virulence factor detection
rule diamond_virulence:
    input:
        contigs="results/virulence/assemblies/{sample}.contigs.fasta"
    output:
        proteins="results/virulence/diamond/{sample}.proteins.faa",
        diamond="results/virulence/diamond/{sample}.virulence.tsv"
    params:
        vfdb_protein=config["vfdb_protein_db"]
    threads: 8
    shell:
        """
        mkdir -p results/virulence/diamond
        
        # Predict proteins from contigs
        prodigal -i {input.contigs} -a {output.proteins} -p meta
        
        # Run DIAMOND against VFDB protein database
        diamond blastp --query {output.proteins} \
                      --db {params.vfdb_protein} \
                      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
                      --threads {threads} \
                      --evalue 1e-10 \
                      --id 80 \
                      --query-cover 80 \
                      --out {output.diamond}
        """

# Combine virulence factor results
rule combine_virulence:
    input:
        vfdb=expand("results/virulence/abricate/{sample}.vfdb.tsv", sample=SAMPLES),
        victors=expand("results/virulence/abricate/{sample}.victors.tsv", sample=SAMPLES),
        srst2=expand("results/virulence/srst2/{sample}.virulence.gene_report.tsv", sample=SAMPLES),
        diamond=expand("results/virulence/diamond/{sample}.virulence.tsv", sample=SAMPLES)
    output:
        combined="results/virulence/reports/combined_virulence_factors.tsv"
    script:
        "scripts/combine_virulence.py"

# Create virulome visualization
rule visualize_virulome:
    input:
        combined="results/virulence/reports/combined_virulence_factors.tsv"
    output:
        heatmap="results/virulence/visualization/virulome_heatmap.png",
        barplot="results/virulence/visualization/virulence_barplot.png",
        network="results/virulence/visualization/virulence_network.png"
    script:
        "scripts/visualize_virulome.py"

# Generate HTML report for virulence analysis
rule virulence_report:
    input:
        combined="results/virulence/reports/combined_virulence_factors.tsv",
        heatmap="results/virulence/visualization/virulome_heatmap.png",
        barplot="results/virulence/visualization/virulence_barplot.png",
        network="results/virulence/visualization/virulence_network.png"
    output:
        report="results/virulence/reports/virulence_profiling_report.html"
    script:
        "scripts/virulence_report.py"
