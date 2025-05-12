# Metagenomics Analysis Pipeline

This Snakemake pipeline implements a comprehensive workflow for metagenomic analysis, following best practices for bacterial detection and novel genome discovery.

## Pipeline Overview

### Pre-processing
1. **Quality Control & Trimming** - Fastp
2. **Host Read Removal** - Bowtie2

### Bacterial Detection Track
1. **Initial Classification** - Kraken2 + Bracken
2. **Functional Profiling** - HUMAnN3

### Novel Bacterial Discovery
1. **Metagenome Assembly** - MEGAHIT
2. **Genome Binning** - MetaBAT2 → DAS Tool (refinement)
3. **Bin Quality Assessment** - CheckM
4. **Taxonomic Placement** - GTDB-Tk
5. **Genome Annotation** - Prokka → DRAM

## Directory Structure

```
.
├── Snakefile
├── config.yaml
├── run_pipeline.sh
├── metagenomics_output/
│   ├── 01_fastp/
│   ├── 02_hostfree/
│   ├── 03_kraken/
│   ├── 04_bracken/
│   ├── 05_humann/
│   ├── 06_megahit/
│   ├── 07_binning/
│   ├── 08_dastool/
│   ├── 09_checkm/
│   ├── 10_gtdbtk/
│   ├── 11_prokka/
│   ├── 12_dram/
│   └── logs/
```

## Quick Start

### 1. Setup Environment

First, make sure you have Snakemake and all the required tools installed. Using conda/mamba is recommended:

```bash
# Create conda environment for the pipeline
conda create -n metagenomics
conda activate metagenomics

# Install snakemake
conda install -c conda-forge -c bioconda snakemake

# Install required tools
conda install -c conda-forge -c bioconda fastp bowtie2 kraken2 bracken humann megahit metabat2 das_tool checkm-genome gtdbtk prokka
```

### 2. Configure the Pipeline

Edit the `config.yaml` file to set up:
- Input sample paths
- Reference genome for host filtering
- Database paths for the various tools
- Resource allocation (threads and memory)

### 3. Prepare Database Files

Several databases need to be set up before running the pipeline:

```bash
# Kraken2 database
kraken2-build --standard --threads 16 --db /path/to/kraken2_db

# HUMAnN3 databases
humann_databases --download chocophlan full /path/to/humann_nucleotide_db
humann_databases --download uniref uniref90_diamond /path/to/humann_protein_db

# GTDB-Tk database
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz -C /path/to/gtdbtk_data
```

### 4. Run the Pipeline

Make the run script executable and run the pipeline:

```bash
# Make script executable
chmod +x run_pipeline.sh

# Run in local mode with 4 cores
./run_pipeline.sh 4

# Or run on a SLURM cluster with 20 concurrent jobs
./run_pipeline.sh 20
```

## Pipeline Steps in Detail

### 1. Quality Control (Fastp)

Fastp performs quality trimming, adapter removal, and filtering of low-quality reads. It generates HTML and JSON reports summarizing the quality metrics.

### 2. Host Read Removal (Bowtie2)

Reads are mapped to the host genome to remove host contamination. Non-mapping reads (microbial reads) are retained for downstream analysis.

### 3. Taxonomic Classification (Kraken2 + Bracken)

Kraken2 performs high-speed taxonomic classification by assigning taxonomy to individual reads. Bracken refines these abundance estimates to improve accuracy at the species and genus levels.

### 4. Functional Profiling (HUMAnN3)

HUMAnN3 (HMP Unified Metabolic Analysis Network) identifies the presence/absence and abundance of microbial pathways in the metagenome, providing insights into the functional potential of the community.

### 5. Metagenome Assembly (MEGAHIT)

MEGAHIT assembles reads into contigs, which are longer sequences that represent fragments of microbial genomes. This memory-efficient assembler is well-suited for metagenomic data.

### 6. Genome Binning (MetaBAT2 → DAS Tool)

MetaBAT2 groups contigs that likely originate from the same genome into bins based on sequence composition and coverage. DAS Tool refines these bins to improve bin quality and completeness.

### 7. Bin Quality Assessment (CheckM)

CheckM evaluates the quality of each genome bin by assessing completeness and contamination based on lineage-specific marker genes.

### 8. Taxonomic Placement (GTDB-Tk)

GTDB-Tk assigns taxonomy to genome bins according to the Genome Taxonomy Database, providing consistent and up-to-date bacterial and archaeal taxonomy.

### 9. Genome Annotation (Prokka → DRAM)

Prokka provides basic structural and functional annotation of the genome bins. DRAM (Distilled and Refined Annotation of Metabolism) extends this annotation with detailed metabolic insights.

## Output Files

Key output files include:

- **Taxonomic profiles**: `04_bracken/{sample}.bracken.species.txt`
- **Functional profiles**: `05_humann/{sample}_pathabundance.tsv`
- **Assembled contigs**: `06_megahit/{sample}_final.contigs.fa`
- **Genome bins**: `08_dastool/{sample}_DASTool_bins/`
- **Bin quality report**: `09_checkm/quality_report.tsv`
- **Taxonomic assignments**: `10_gtdbtk/gtdbtk.summary.tsv`
- **Metabolic analysis**: `12_dram/distill/product.html`

## Customization

You can customize the pipeline by:

1. Editing the Snakefile to add/remove steps
2. Modifying tool parameters in the rule definitions
3. Adjusting resource allocation in the config.yaml file
4. Adding pre/post-processing steps as needed

## Troubleshooting

1. **Missing dependencies**: Ensure all tools are properly installed and in your PATH
2. **Database issues**: Verify that all databases are downloaded and properly configured
3. **Resource limitations**: Adjust memory and CPU settings in config.yaml
4. **File not found errors**: Check input file paths and permissions

## Citation

If you use this pipeline, please cite the following tools:

- Fastp: Chen S, et al. (2018). Bioinformatics.
- Bowtie2: Langmead B, et al. (2012). Nature Methods.
- Kraken2: Wood DE, et al. (2019). Genome Biology.
- Bracken: Lu J, et al. (2017). PeerJ Computer Science.
- HUMAnN3: Beghini F, et al. (2021). eLife.
- MEGAHIT: Li D, et al. (2015). Bioinformatics.
- MetaBAT2: Kang DD, et al. (2019). PeerJ.
- DAS Tool: Sieber CMK, et al. (2018). Nature Microbiology.
- CheckM: Parks DH, et al. (2015). Genome Research.
- GTDB-Tk: Chaumeil PA, et al. (2019). Bioinformatics.
- Prokka: Seemann T. (2014). Bioinformatics.
- DRAM: Shaffer M, et al. (2020). Nucleic Acids Research.
