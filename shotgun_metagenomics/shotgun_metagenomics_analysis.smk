# Shotgun Metagenomics Analysis Pipeline
# Includes:
# 1. Taxonomic profiling (Kraken2, MetaPhlAn)
# 2. Antibiotic resistance profiling (ResFinder, DeepARG)
# 3. Strain-level analysis and outbreak tracking

# Configuration
configfile: "config.yaml"

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
        # Quality Control
        expand("results/fastqc/{sample}_R{read}_fastqc.html", sample=SAMPLES, read=[1, 2]),
        expand("results/trimmed/{sample}_R{read}.trimmed.fastq.gz", sample=SAMPLES, read=[1, 2]),
        
        # Taxonomic profiling
        expand("results/kraken2/{sample}.kraken2.report", sample=SAMPLES),
        expand("results/metaphlan/{sample}.metaphlan.profile.txt", sample=SAMPLES),
        "results/metaphlan/merged_abundance_table.txt",
        
        # Antimicrobial resistance
        expand("results/resfinder/{sample}.resfinder.json", sample=SAMPLES),
        expand("results/deeparg/{sample}.deeparg.mapping.ARG", sample=SAMPLES),
        "results/deeparg/combined_deeparg_abundance.txt",
        
        # Assembly for strain analysis
        expand("results/assembly/{sample}/contigs.fasta", sample=SAMPLES),
        
        # Strain-level analysis
        expand("results/strain_analysis/{sample}.straingst.json", sample=SAMPLES),
        expand("results/mlst/{sample}.mlst.tsv", sample=SAMPLES),
        "results/strain_tracking/outbreak_clusters.tsv"

# Quality Control
rule fastqc:
    input:
        unpack(get_read_paths)
    output:
        html_r1="results/fastqc/{sample}_R1_fastqc.html",
        html_r2="results/fastqc/{sample}_R2_fastqc.html",
        zip_r1="results/fastqc/{sample}_R1_fastqc.zip",
        zip_r2="results/fastqc/{sample}_R2_fastqc.zip"
    threads: 2
    shell:
        """
        mkdir -p results/fastqc
        fastqc {input.r1} {input.r2} -o results/fastqc -t {threads}
        """

rule trim_reads:
    input:
        unpack(get_read_paths)
    output:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        unpaired_r1="results/trimmed/{sample}_R1.unpaired.fastq.gz",
        unpaired_r2="results/trimmed/{sample}_R2.unpaired.fastq.gz"
    threads: 4
    shell:
        """
        mkdir -p results/trimmed
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {output.unpaired_r1} \
            {output.r2} {output.unpaired_r2} \
            ILLUMINACLIP:{config[adapter_file]}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

# Taxonomic Profiling

## 1. Kraken2 for taxonomic classification
rule kraken2:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        report="results/kraken2/{sample}.kraken2.report",
        output="results/kraken2/{sample}.kraken2.output"
    threads: 8
    shell:
        """
        mkdir -p results/kraken2
        kraken2 --paired \
            --db {config[kraken2_db]} \
            --threads {threads} \
            --output {output.output} \
            --report {output.report} \
            {input.r1} {input.r2}
        """

## 2. MetaPhlAn for species abundance
rule metaphlan:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        profile="results/metaphlan/{sample}.metaphlan.profile.txt",
        bowtie2out="results/metaphlan/{sample}.metaphlan.bowtie2.bz2"
    threads: 8
    shell:
        """
        mkdir -p results/metaphlan
        metaphlan {input.r1},{input.r2} \
            --input_type fastq \
            --bowtie2out {output.bowtie2out} \
            --nproc {threads} \
            --output_file {output.profile} \
            --bowtie2db {config[metaphlan_db]} \
            --unknown_estimation \
            --biom results/metaphlan/{wildcards.sample}.metaphlan.biom
        """

## 3. Merge MetaPhlAn results
rule merge_metaphlan:
    input:
        expand("results/metaphlan/{sample}.metaphlan.profile.txt", sample=SAMPLES)
    output:
        "results/metaphlan/merged_abundance_table.txt"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """

# Antimicrobial Resistance Profiling

## 1. ResFinder for acquired resistance
rule resfinder:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        json="results/resfinder/{sample}.resfinder.json"
    params:
        outdir="results/resfinder/{sample}"
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}
        run_resfinder.py \
            --inputfastq_1 {input.r1} \
            --inputfastq_2 {input.r2} \
            --acquired \
            --outputPath {params.outdir} \
            --db_path_res {config[resfinder_db]} \
            --species {config[species]} \
            --threads {threads}
        mv {params.outdir}/ResFinder_results_tab.txt {params.outdir}/{wildcards.sample}.resfinder.txt
        mv {params.outdir}/ResFinder_results.json {output.json}
        """

## 2. DeepARG for broader resistance gene detection
rule deeparg:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        merged="results/deeparg/{sample}.merged.fastq",
        arg="results/deeparg/{sample}.deeparg.mapping.ARG",
        potential_arg="results/deeparg/{sample}.deeparg.mapping.potential.ARG"
    params:
        prefix="results/deeparg/{sample}.deeparg"
    threads: 8
    shell:
        """
        mkdir -p results/deeparg
        
        # Merge paired reads for DeepARG
        seqtk mergepe {input.r1} {input.r2} > {output.merged}
        
        # Run DeepARG
        deeparg predict \
            --model LS \
            --reads {output.merged} \
            --out {params.prefix} \
            --type nucl \
            --min-prob 0.8 \
            --arg-alignment-identity 90 \
            --arg-alignment-evalue 1e-10 \
            --arg-num-alignments-per-entry 1000 \
            --arg-alignment-overlap 0.8 \
            --threads {threads}
        """

## 3. Combine DeepARG results
rule combine_deeparg:
    input:
        expand("results/deeparg/{sample}.deeparg.mapping.ARG", sample=SAMPLES)
    output:
        "results/deeparg/combined_deeparg_abundance.txt"
    script:
        "scripts/combine_deeparg.py"

# Assembly for strain-level analysis

## 1. Assemble reads using SPAdes
rule assembly:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        contigs="results/assembly/{sample}/contigs.fasta",
        scaffolds="results/assembly/{sample}/scaffolds.fasta"
    threads: 16
    params:
        outdir="results/assembly/{sample}"
    shell:
        """
        mkdir -p {params.outdir}
        spades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            --meta \
            -o {params.outdir} \
            -t {threads}
        """

# Strain-level analysis

## 1. StrainGST for strain identification
rule straingst:
    input:
        profile="results/metaphlan/{sample}.metaphlan.profile.txt"
    output:
        json="results/strain_analysis/{sample}.straingst.json"
    threads: 4
    shell:
        """
        mkdir -p results/strain_analysis
        straingst \
            --metaphlan_profile {input.profile} \
            --output {output.json} \
            --database {config[straingst_db]} \
            --threads {threads}
        """

## 2. MLST for sequence typing
rule mlst:
    input:
        contigs="results/assembly/{sample}/contigs.fasta"
    output:
        tsv="results/mlst/{sample}.mlst.tsv"
    shell:
        """
        mkdir -p results/mlst
        mlst {input.contigs} > {output.tsv}
        """

## 3. Outbreak tracking via SNP analysis (for specific pathogens)
rule outbreak_tracking:
    input:
        straingst=expand("results/strain_analysis/{sample}.straingst.json", sample=SAMPLES),
        mlst=expand("results/mlst/{sample}.mlst.tsv", sample=SAMPLES),
        assemblies=expand("results/assembly/{sample}/contigs.fasta", sample=SAMPLES)
    output:
        clusters="results/strain_tracking/outbreak_clusters.tsv",
        phylogeny="results/strain_tracking/outbreak_phylogeny.nwk"
    params:
        target_species=config["outbreak_species_of_interest"]
    script:
        "scripts/outbreak_tracking.py"

# Produce summary reports
rule generate_reports:
    input:
        kraken=expand("results/kraken2/{sample}.kraken2.report", sample=SAMPLES),
        metaphlan="results/metaphlan/merged_abundance_table.txt",
        resfinder=expand("results/resfinder/{sample}.resfinder.json", sample=SAMPLES),
        deeparg="results/deeparg/combined_deeparg_abundance.txt",
        strain="results/strain_tracking/outbreak_clusters.tsv"
    output:
        taxonomy_report="results/reports/taxonomy_summary.html",
        amr_report="results/reports/amr_summary.html",
        strain_report="results/reports/strain_summary.html",
        outbreak_report="results/reports/outbreak_tracking.html"
    script:
        "scripts/generate_reports.py"
