# QIIME2 16S rRNA V3-V4 Metagenomics Pipeline

This Snakemake-based pipeline provides a comprehensive workflow for analyzing 16S rRNA gene sequencing data (V3-V4 region) using QIIME2. The pipeline includes quality control, taxonomic assignment, diversity analyses, and differential abundance testing.

## Features

- **Quality Control**: DADA2 and Deblur for denoising and ASV generation
- **Taxonomic Assignment**: Naive Bayes classifier with SILVA/Greengenes databases
- **Diversity Analyses**: Alpha diversity (Shannon, Simpson, Faith's PD) and Beta diversity (Bray-Curtis, UniFrac)
- **Differential Abundance**: ANCOM, DESeq2, and LEfSe methods
- **Snakemake Workflow**: Efficient, reproducible, and scalable analysis pipeline
- **Manifest File Input**: Simple sample management with a CSV manifest file

## Requirements

- QIIME2 (2023.5 or newer)
- Snakemake (7.0 or newer)
- R with the following packages:
  - biomformat
  - DESeq2
  - phyloseq
  - tidyverse
- Python dependencies:
  - pandas
  - numpy
  - scikit-learn
  - BioPython
- LEfSe (for differential abundance analysis)

## Setup

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/qiime2-16s-pipeline.git
   cd qiime2-16s-pipeline
   ```

2. Create a conda environment with the required dependencies:
   ```
   conda env create -f environment.yml
   conda activate qiime2-pipeline
   ```

3. Prepare your manifest and metadata files:
   - Edit `manifest.csv` to include your sample information
   - Edit `metadata.tsv` to include your sample metadata

4. Customize the configuration:
   - Edit `config.yaml` to match your experimental setup

## Pipeline Overview

The pipeline performs the following steps:

1. **Data Import**: Imports paired-end sequencing data using a manifest file
2. **Reference Database Download**: Retrieves and prepares SILVA/Greengenes reference databases
3. **Quality Control**: Processes reads using DADA2 or Deblur to generate ASVs
4. **Taxonomic Classification**: Assigns taxonomy using a trained Naive Bayes classifier
5. **Phylogenetic Analysis**: Creates a phylogenetic tree for diversity analyses
6. **Diversity Analyses**: Calculates alpha and beta diversity metrics
7. **Differential Abundance**: Identifies differentially abundant taxa between groups

## Usage

1. Activate the conda environment:
   ```
   conda activate qiime2-pipeline
   ```

2. Run the pipeline:
   ```
   snakemake --cores 8 --use-conda
   ```

3. To generate a workflow report:
   ```
   snakemake --report workflow_report.html
   ```

4. To visualize the workflow:
   ```
   snakemake --dag | dot -Tpng > dag.png
   ```

## Custom Naive Bayes Classifier

The pipeline includes a standalone scikit-learn implementation of the Naive Bayes classifier for taxonomic assignment. This can be used independently of QIIME2 if needed:

```
# Train a classifier
python scripts/sklearn_classifier.py train \
  --reference-seqs silva-138-99-v3v4-seqs.fasta \
  --reference-tax silva-138-99-tax.tsv \
  --output-classifier silva-v3v4-classifier.pkl

# Classify sequences
python scripts/sklearn_classifier.py classify \
  --classifier silva-v3v4-classifier.pkl \
  --input-seqs rep-seqs.fasta \
  --output-tax taxonomy_assignments.tsv
```

## Outputs

The pipeline generates organized outputs in the following directory structure:

```
results/
├── 1_qc/
│   ├── dada2/
│   └── deblur/
├── 2_taxonomy/
├── 3_diversity/
└── 4_differential_abundance/
```

## Customization

The pipeline is highly customizable through the `config.yaml` file. Key parameters include:

- Primer sequences for the V3-V4 region
- DADA2 and Deblur parameters
- Diversity analysis parameters
- Differential abundance settings

## License

This pipeline is available under the MIT License.

## Citation

If you use this pipeline in your research, please cite the following tools:

- QIIME2: Bolyen E, et al. (2019) Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology.
- DADA2: Callahan BJ, et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods.
- Snakemake: Köster J & Rahmann S (2012) Snakemake - A scalable bioinformatics workflow engine. Bioinformatics.
