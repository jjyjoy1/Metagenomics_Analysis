# MetaVirome Pipeline: User Tutorial

Welcome to the MetaVirome Pipeline tutorial. This document provides a step-by-step guide to using the pipeline for your metagenomics studies focused on virus detection.

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Step-by-Step Tutorial](#step-by-step-tutorial)
   - [Preparing Input Data](#preparing-input-data)
   - [Creating Databases](#creating-databases)
   - [Running the Pipeline](#running-the-pipeline)
   - [Understanding the Results](#understanding-the-results)
   - [Analyzing the Results](#analyzing-the-results)
5. [Case Studies](#case-studies)
   - [Environmental Metagenomics](#environmental-metagenomics)
   - [Clinical Virome Analysis](#clinical-virome-analysis)
   - [Phage Discovery](#phage-discovery)
6. [Troubleshooting](#troubleshooting)
7. [FAQs](#faqs)

## Introduction

The MetaVirome Pipeline is a comprehensive Nextflow-based workflow for virus detection and analysis in metagenomic data. It integrates multiple state-of-the-art tools to provide a complete solution for viral metagenomics:

1. **Host DNA Removal**: Removes host contamination to focus on microbial and viral content
2. **Known Virus Detection**: Rapidly screens for known viruses using Kraken2
3. **Sensitive Viral Detection**: Detects divergent viral proteins using DIAMOND
4. **Viral Assembly**: Assembles viral genomes using metaSPAdes
5. **Viral Contig Identification**: Identifies viral sequences using VirSorter2
6. **Viral Quality Control**: Assesses viral genome quality using CheckV
7. **Viral Annotation**: Annotates viral genomes using Pharokka
8. **Viral Taxonomy**: Places novel viruses in taxonomic context using vConTACT2

## Installation

### Prerequisites

- Nextflow (v20.04.0 or later)
- Java 8 or later
- Docker or Singularity (recommended for reproducibility)
- At least 16GB of RAM for standard datasets (64GB+ recommended for complex samples)
- At least 100GB of disk space for databases and outputs

### Setting Up the Environment

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/metavirome.git
   cd metavirome
   ```

2. **Install Nextflow** (if not already installed):
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

3. **Choose a container technology**:
   
   The pipeline supports Docker, Singularity, and Conda. We recommend using containers for maximum reproducibility.
   
   **Docker**:
   ```bash
   # Install Docker
   sudo apt-get update
   sudo apt-get install docker-ce docker-ce-cli containerd.io
   ```
   
   **Singularity**:
   ```bash
   # Install Singularity
   sudo apt-get update
   sudo apt-get install singularity-container
   ```

## Quick Start

The quickest way to try the pipeline is with the included test dataset:

```bash
# Run the test dataset
./run_pipeline.sh --test

# Run test and generate analysis report
./run_pipeline.sh --test --analyze
```

This will download a small test dataset, set up minimal databases, and run the complete pipeline to verify everything is working correctly.

## Step-by-Step Tutorial

### Preparing Input Data

The pipeline accepts Illumina paired-end or single-end sequencing data in FASTQ format (compressed or uncompressed).

1. **Organize your data**:
   
   Place your FASTQ files in a dedicated directory with a consistent naming pattern:
   ```
   data/
   ├── sample1_R1.fastq.gz
   ├── sample1_R2.fastq.gz
   ├── sample2_R1.fastq.gz
   ├── sample2_R2.fastq.gz
   ...
   ```

2. **Prepare host genome**:
   
   You'll need a FASTA file of the host genome for the host DNA removal step. For human samples, download the latest reference:
   ```bash
   mkdir -p reference
   wget -O reference/GRCh38.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
   gunzip reference/GRCh38.fa.gz
   ```

### Creating Databases

The pipeline requires two main databases: Kraken2 and DIAMOND.

1. **Kraken2 viral database**:
   
   ```bash
   mkdir -p databases/kraken2
   cd databases/kraken2
   
   # Download NCBI taxonomy
   kraken2-build --download-taxonomy --db kraken2_viral_db
   
   # Download viral library
   kraken2-build --download-library viral --db kraken2_viral_db
   
   # Build the database
   kraken2-build --build --db kraken2_viral_db --threads 8
   
   cd ../..
   ```

2. **DIAMOND viral protein database**:
   
   ```bash
   mkdir -p databases/diamond
   cd databases/diamond
   
   # Download NCBI viral proteins
   wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
   wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
   gunzip viral.*.protein.faa.gz
   
   # Concatenate files
   cat viral.*.protein.faa > viral_proteins.faa
   
   # Build DIAMOND database
   diamond makedb --in viral_proteins.faa --db viral_proteins
   
   cd ../..
   ```

### Running the Pipeline

Once your data and databases are prepared, you can run the pipeline:

```bash
./run_pipeline.sh \
  --input "data/*_R{1,2}.fastq.gz" \
  --host reference/GRCh38.fa \
  --kraken databases/kraken2/kraken2_viral_db \
  --diamond databases/diamond/viral_proteins \
  --outdir results \
  --profile docker \
  --cpus 16 \
  --memory 64.GB
```

For single-end data, add the `--single-end` flag:

```bash
./run_pipeline.sh \
  --input "data/*.fastq.gz" \
  --single-end \
  --host reference/GRCh38.fa \
  --kraken databases/kraken2/kraken2_viral_db \
  --diamond databases/diamond/viral_proteins \
  --outdir results
```

### Understanding the Results

The pipeline produces a structured output directory:

```
results/
├── 01_host_removed/            # Host-filtered reads
├── 02_kraken2/                 # Kraken2 classification results
├── 03_diamond/                 # DIAMOND protein alignments
├── 04_assembly/                # metaSPAdes assemblies
├── 05_virsorter2/              # VirSorter2 viral predictions
├── 06_checkv/                  # CheckV quality assessment
├── 07_annotation/              # Pharokka annotations
├── 08_taxonomy/                # vConTACT2 taxonomy results
└── pipeline_reports/           # Nextflow execution reports
```

Key files to examine:

1. **Kraken2 reports** (`02_kraken2/{sample}_kraken2.report`):
   - Shows the taxonomic classification of sequences
   - Visualize with Krona: `ktImportTaxonomy -q 2 -t 3 {sample}_kraken2.output -o {sample}_krona.html`

2. **High-quality viral contigs** (`06_checkv/{sample}_high_quality.fasta`):
   - Contains viral contigs filtered for quality
   - These are the most reliable viral sequences for further analysis

3. **Viral annotation** (`07_annotation/{sample}_pharokka/{sample}_pharokka.gff`):
   - GFF file with gene predictions and annotations
   - Use tools like Artemis for visualization: `art {sample}_pharokka.gff`

4. **Viral taxonomy** (`08_taxonomy/{sample}_vcontact/genome_by_genome_overview.csv`):
   - Shows clustering information for viral contigs
   - Use with Cytoscape for network visualization

### Analyzing the Results

To generate a comprehensive analysis report from the pipeline results:

```bash
./run_pipeline.sh --analyze
```

This will create an HTML report (`final_report/MetaVirome_Analysis_Report.html`) with visualizations and statistics about your viral communities.

## Case Studies

### Environmental Metagenomics

For environmental samples (e.g., marine, soil):

1. **Remove appropriate host**:
   - For marine samples, consider using common marine organisms as "hosts"
   - For soil, consider plant genomes relevant to your sampling site

2. **Adjust parameters**:
   ```bash
   ./run_pipeline.sh \
     --input "environmental_samples/*_R{1,2}.fastq.gz" \
     --host reference/plant_genome.fa \
     --kraken databases/kraken2/kraken2_viral_db \
     --diamond databases/diamond/viral_proteins \
     --profile singularity \
     --cpus 32 \
     --memory 128.GB
   ```

3. **Focus on specific viral types**:
   - For phages in environmental samples, examine VirSorter2 results in `results/05_virsorter2/{sample}_virsorter/final-viral-score.tsv`
   - Look for "dsDNAphage" classifications

### Clinical Virome Analysis

For clinical samples:

1. **Use human host genome**:
   ```bash
   ./run_pipeline.sh \
     --input "clinical_samples/*_R{1,2}.fastq.gz" \
     --host reference/GRCh38.fa \
     --kraken databases/kraken2/kraken2_viral_db \
     --diamond databases/diamond/viral_proteins \
     --profile docker
   ```

2. **Focus on known pathogenic viruses**:
   - Examine Kraken2 reports for known human pathogens
   - Cross-reference with DIAMOND results for confirmation

3. **Examine potential novel viruses**:
   - Look for high-quality contigs that don't match known references
   - Focus on those with high CheckV completeness scores

### Phage Discovery

For phage discovery studies:

1. **Use bacterial host genome**:
   ```bash
   ./run_pipeline.sh \
     --input "phage_samples/*_R{1,2}.fastq.gz" \
     --host reference/bacteria.fa \
     --kraken databases/kraken2/kraken2_viral_db \
     --diamond databases/diamond/viral_proteins
   ```

2. **Focus on Pharokka results**:
   - Examine `results/07_annotation/{sample}_pharokka/` for phage-specific annotations
   - Look for key phage genes (terminase, capsid, tail, etc.)

3. **Examine genome completeness**:
   - Use CheckV results to assess genome completeness
   - Complete phage genomes will often show circular topology

## Troubleshooting

### Common Issues and Solutions

1. **Out of memory errors**:
   - Increase available memory: `--memory 128.GB`
   - For metaSPAdes, use a machine with at least 64GB RAM for complex samples

2. **Excessive runtime**:
   - Use containers (Docker/Singularity) for better performance
   - Split analysis into batches of samples
   - For initial testing, subsample your data

3. **Database issues**:
   - Ensure database paths are absolute, not relative
   - Verify database files are readable by the process user
   - For Kraken2, ensure you have both taxonomy and library files

4. **Missing viral contigs**:
   - Check host removal rate - if too aggressive, you might lose viral reads
   - Adjust VirSorter2 parameters for higher sensitivity

### Pipeline Failure Debugging

If the pipeline fails:

1. **Check work directory**:
   - Examine `.command.err` files in the failed task directory
   - Look for specific error messages

2. **Restart with resume**:
   ```bash
   ./run_pipeline.sh --input "data/*_R{1,2}.fastq.gz" ... -resume
   ```

3. **Increase resources**:
   - Adjust CPU, memory, and time limits in `nextflow.config`

## FAQs

**Q: How much disk space do I need?**
A: Plan for at least 100GB for databases and 1-2TB for a full analysis of multiple samples.

**Q: How long does the pipeline take to run?**
A: Runtime varies by sample complexity and available resources. For a typical metagenomic sample (~10GB FASTQ):
- Host removal: 1-2 hours
- Kraken2 + DIAMOND: 1-3 hours
- Assembly: 2-12 hours
- VirSorter2: 1-3 hours
- CheckV + Pharokka: 1-2 hours
- vConTACT2: 1-4 hours
Total: 6-24 hours per sample with 16 CPUs

**Q: Can I run just part of the pipeline?**
A: Yes, you can modify `metavirome.nf` to include only specific processes, though this requires Nextflow knowledge.

**Q: How many viruses can I expect to find?**
A: This varies dramatically by sample type:
- Environmental: 10-1000+ viral contigs
- Clinical: 1-50 viral contigs
- Enriched phage samples: 10-300 viral contigs

**Q: How do I cite this pipeline?**
A: Please cite both the pipeline itself and the individual tools:
- Kraken2: Wood, D.E., et al. (2019)
- DIAMOND: Buchfink, B., et al. (2021)
- metaSPAdes: Nurk, S., et al. (2017)
- VirSorter2: Guo, J., et al. (2021)
- CheckV: Nayfach, S., et al. (2021)
- Pharokka: Michniewski, S., et al. (2022)
- vConTACT2: Jang, H.B., et al. (2019)
