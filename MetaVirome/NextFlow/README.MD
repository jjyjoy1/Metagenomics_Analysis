# MetaVirome Pipeline

A Nextflow pipeline for metagenomics analysis focused on virus detection, with host DNA removal.

## Overview

This pipeline implements a comprehensive workflow for viral metagenomic analysis, including:

1. **Preprocessing**: Removal of host DNA
2. **Known Virus Detection**: Fast screening with Kraken2
3. **Sensitive Viral Detection**: Protein-level detection with DIAMOND
4. **Viral Assembly**: Using metaSPAdes with metagenome-specific settings
5. **Viral Contig Identification**: Using VirSorter2 for accurate viral contig prediction
6. **Viral Quality Control**: Using CheckV to assess completeness and quality
7. **Viral Annotation**: Using Pharokka for phage annotation
8. **Novelty Assessment**: Using vConTACT2 for network-based taxonomy

## Requirements

- Nextflow (version 20.04.0 or later)
- Docker or Singularity (recommended for reproducibility)
- The following tools (can be installed in a conda environment):
  - Bowtie2
  - Kraken2 (with viral database)
  - DIAMOND (with viral protein database)
  - metaSPAdes
  - VirSorter2
  - CheckV
  - Pharokka
  - Prodigal
  - vConTACT2
  - seqtk

## Setup

### Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

### Install Required Tools

Using conda (recommended):

```bash
# Create conda environment
conda create -n metavirome
conda activate metavirome

# Install tools
conda install -c bioconda bowtie2 kraken2 diamond spades virsorter2 checkv prodigal seqtk
pip install pharokka
pip install vcontact2
```

### Prepare Databases

1. **Kraken2 Viral Database**:
   ```bash
   kraken2-build --download-library viral --db kraken2_viral_db
   kraken2-build --build --db kraken2_viral_db
   ```

2. **DIAMOND Viral Protein Database**:
   ```bash
   # Download NCBI RefSeq viral proteins
   wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.*.protein.faa.gz
   gunzip viral.*.protein.faa.gz
   cat viral.*.protein.faa > viral_proteins.faa
   
   # Format DIAMOND database
   diamond makedb --in viral_proteins.faa --db viral_proteins
   ```

## Running the Pipeline

### Basic Usage

```bash
nextflow run metavirome.nf \
  --reads "path/to/reads/*_R{1,2}.fastq.gz" \
  --host_genome "path/to/host.fasta" \
  --kraken2_db "path/to/kraken2_viral_db" \
  --diamond_db "path/to/diamond_viral_proteins"
```

### Options

```
Required arguments:
  --reads          Path to input reads (must be enclosed in quotes and include the glob pattern)
  --host_genome    Path to host genome for host DNA removal
  --kraken2_db     Path to Kraken2 viral database
  --diamond_db     Path to DIAMOND viral protein database

Optional arguments:
  --outdir         Output directory (default: ./results)
  --single_end     Specify if reads are single-end (default: false)
  --min_contig_len Minimum contig length for viral analysis (default: 1000)
  
Resource allocation:
  --max_memory     Maximum memory to use (default: 128.GB)
  --max_cpus       Maximum CPUs to use (default: 16)
  --max_time       Maximum time to run (default: 240.h)
```

### Example Command

```bash
nextflow run metavirome.nf \
  --reads "/data/samples/*_R{1,2}.fastq.gz" \
  --host_genome "/data/references/human_genome.fasta" \
  --kraken2_db "/data/databases/kraken2_viral_db" \
  --diamond_db "/data/databases/viral_proteins" \
  --outdir "/data/results/metavirome" \
  --max_cpus 24
```

## Output Structure

The pipeline produces a structured output directory:

```
results/
├── 01_host_removed/            # Host-filtered reads
├── 02_kraken2/                 # Kraken2 classification results
├── 03_diamond/                 # DIAMOND blastx results
├── 04_assembly/                # metaSPAdes assembly output
├── 05_virsorter2/              # VirSorter2 viral predictions
├── 06_checkv/                  # CheckV quality assessment
├── 07_annotation/              # Pharokka annotations
└── 08_taxonomy/                # vConTACT2 taxonomy results
```

## Performance Considerations

- Host removal is often memory-intensive, especially with large host genomes
- metaSPAdes assembly requires significant memory (recommend 64GB+ for complex metagenomes)
- VirSorter2 and vConTACT2 can be computationally intensive

## Customization

The pipeline can be customized by:

1. Modifying resource allocation parameters
2. Adding additional databases or tools
3. Adjusting filtering thresholds for viral contigs

You can also configure Nextflow to use container technologies like Docker or Singularity by creating a `nextflow.config` file.

## Citation

If you use this pipeline, please cite the tools it incorporates:

- Kraken2: Wood, D.E., et al. (2019)
- DIAMOND: Buchfink, B., et al. (2021)
- metaSPAdes: Nurk, S., et al. (2017)
- VirSorter2: Guo, J., et al. (2021)
- CheckV: Nayfach, S., et al. (2021)
- Pharokka: Michniewski, S., et al. (2022)
- vConTACT2: Jang, H.B., et al. (2019)
