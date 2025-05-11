# Novel Pathogen Detection Pipeline

This Snakemake pipeline is designed for comprehensive analysis of NGS data to detect and characterize novel bacterial and viral pathogens. The pipeline integrates multiple methods for assembly, taxonomic classification, viral detection, and anomaly detection using embedding models.

## Overview

The pipeline is organized into several modules:

1. **Assembly & Binning**: Quality control, host sequence removal, assembly, read mapping, and binning
2. **Taxonomic Classification**: Gene prediction, taxonomic classification using CAT/BAT, GTDB-Tk, and ANI analysis
3. **Viral Detection**: Detection, verification, and classification of viral sequences
4. **Embedding Models**: Generation of protein and DNA embeddings for anomaly detection
5. **Integration**: Anomaly detection and comprehensive reporting

## Requirements

- Snakemake (version 6.0+)
- Python (3.8+) with BioPython, pandas, numpy
- Bioinformatics tools:
  - FastQC
  - Bowtie2
  - MEGAHIT
  - BWA
  - Samtools
  - MetaBAT2
  - GraphBin
  - Prodigal
  - CAT/BAT
  - MMseqs2
  - GTDB-Tk
  - FastANI
  - VirSorter2
  - VIBRANT
  - DeepVirFinder
  - CheckV
  - vConTACT2
- Embedding models:
  - Protein models: ESM2, ProtT5, ESM1b
  - DNA models: DNABERT, Nucleotide Transformer

## Installation

The pipeline assumes that all required tools are installed and available in your environment. You may use conda, modules, or other environment management systems to make these tools available.

## Configuration

Edit the `config.yaml` file to set:
- Sample IDs
- Output directory
- Reference host genome
- Database paths
- Parameters for embedding models and anomaly detection

## Input Data

Place your input FASTQ files in the `data` directory following the naming convention:
- `sample_R1.fastq.gz`: Forward reads
- `sample_R2.fastq.gz`: Reverse reads

Where `sample` corresponds to the sample IDs listed in your config file.

## Running the Pipeline

1. Edit the `config.yaml` file with your settings
2. Run the pipeline with:

```bash
snakemake --cores [number_of_cores] --use-conda
```

## Output

The pipeline generates results in the specified output directory with the following structure:

```
results/
├── qc/                       # Quality control results
├── cleaned_reads/            # Host-removed reads
├── assembly/                 # Assembled contigs
├── mapping/                  # Read mapping results
├── bins/                     # Initial binning results
├── graphbin/                 # Refined bins
├── genes/                    # Gene predictions
├── cat/                      # CAT taxonomic classifications
├── bat/                      # BAT bin classifications
├── mmseqs2/                  # MMseqs2 classification results
├── gtdbtk/                   # GTDB-Tk phylogenetic placements
├── ani/                      # ANI analysis results
├── virsorter2/               # VirSorter2 results
├── vibrant/                  # VIBRANT results
├── deepvirfinder/            # DeepVirFinder results
├── checkv/                   # CheckV results
├── vcontact2/                # vConTACT2 results
├── embeddings/               # Protein and DNA embeddings
├── anomaly_detection/        # Anomaly detection results
└── reports/                  # Comprehensive reports
```

## Helper Scripts

The pipeline relies on several Python helper scripts:

- `run_protein_embedding.py`: Generates protein embeddings
- `run_dna_embedding.py`: Generates DNA embeddings
- `run_anomaly_detection.py`: Runs anomaly detection
- `generate_comprehensive_report.py`: Creates final reports

You need to implement these scripts according to your specific requirements.

## Customization

You can customize the pipeline by:
- Adding new tools or steps
- Modifying parameters in the config file
- Extending the anomaly detection methods
- Customizing the report generation


## Main Snakemake Pipeline File
The main Snakefile contains all the workflow rules organized into the sections you requested:

###Assembly & Binning:

Quality control with FastQC
Host sequence removal with Bowtie2
Assembly with MEGAHIT
Read mapping with BWA and Samtools
Binning with MetaBAT2
Bin refinement with GraphBin


###Taxonomic Classification:

Gene prediction with Prodigal
Contig classification with CAT
Bin classification with BAT
MMseqs2 for protein classification
GTDB-Tk for bacterial phylogenetic placement
ANI analysis for species delineation


###Viral Detection:

VirSorter2 for viral sequence detection
VIBRANT for viral detection and annotation
DeepVirFinder for additional viral detection
CheckV for viral quality assessment
vConTACT2 for viral genome clustering and taxonomy


###Embedding Models:

Protein embedding generation (ESM2, ProtT5, ESM1b)
DNA sequence embedding generation (DNABERT, Nucleotide Transformer)


###Integration:

Anomaly detection using multiple methods
Comprehensive report generation



##Configuration File
The config.yaml file contains all configurable parameters:

Sample IDs
Output directory
Thread count
Minimum contig size
Database paths
Embedding model selection
Anomaly detection methods

