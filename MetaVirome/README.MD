# Metagenomics Viral Detection Pipeline

This Snakemake workflow implements a comprehensive pipeline for virus detection in metagenomic data, combining tools for both known and novel virus discovery with a focus on completeness and accuracy.

## Overview

This pipeline provides:

1. **Preprocessing**
   - Quality control with FastP
   - Host DNA removal with Bowtie2

2. **Known Virus Detection**
   - Kraken2 with viral database for fast classification

3. **Sensitive Viral Detection**
   - DIAMOND alignment against viral protein database

4. **Novel Virus Discovery**
   - metaSPAdes assembly with `--meta` flag
   - VirSorter2 for viral contig identification
   - CheckV for viral quality assessment
   - Pharokka and VGAS for viral annotation
   - vConTACT2 for viral taxonomy classification

## Pipeline Structure

```
workflow/
├── Snakefile              # Main workflow file
├── config.yaml            # Configuration file
├── envs/                  # Conda environment files
│   ├── fastp.yaml
│   ├── bowtie2.yaml
│   └── ...
├── scripts/               # Helper scripts
│   ├── extract_viral_contigs.py
│   ├── extract_quality_viruses.py
│   ├── create_gene2genome.py
│   └── generate_report.py
└── results/               # Generated results folder
    ├── preprocessing/
    ├── kraken2/
    ├── assembly/
    ├── diamond/
    ├── virsorter2/
    ├── checkv/
    ├── pharokka/
    ├── vgas/
    ├── vcontact2/
    └── final_report.html
```

## Prerequisites

1. **Snakemake**: This pipeline requires Snakemake to run
   ```
   conda install -c bioconda -c conda-forge snakemake
   ```

2. **Required Databases**: Make sure to download and prepare the following databases:
   - Host genome reference (for host removal)
   - Kraken2 viral database
   - DIAMOND viral protein database
   - CheckV database
   - Pharokka database

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/metagenomics-viral-pipeline.git
   cd metagenomics-viral-pipeline
   ```

2. Create the required directory structure:
   ```
   mkdir -p raw_data logs results/{preprocessing,qc,kraken2,assembly,diamond,virsorter2,checkv,pharokka,vgas,vcontact2}
   ```

3. Install Conda environments:
   ```
   # Create all conda environments
   snakemake --use-conda --conda-create-envs-only --cores 1
   ```

## Configuration

Edit the `config.yaml` file to specify:

1. **Samples**: List of sample names to process (filenames should be in format `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`)
2. **Database paths**: Update paths to all required databases
3. **Resource configuration**: Adjust threads and memory as needed

Example:
```yaml
# Configuration file for metagenomics viral detection pipeline
samples:
  - "sample1"
  - "sample2"
  - "sample3"

# Reference databases paths
host_genome: "/path/to/host/genome/bowtie2_index/host"
kraken2_db: "/path/to/kraken2/viral_db"
diamond_viral_db: "/path/to/diamond/viral_protein_db/viral_proteins"
checkv_db: "/path/to/checkv/db"
pharokka_db: "/path/to/pharokka/db"

# Resource configurations
max_threads: 32
max_memory: "128G"
```

## Running the Pipeline

1. Place your raw sequencing files in the `raw_data/` directory with the naming convention:
   - `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`
   - `sample2_R1.fastq.gz`, `sample2_R2.fastq.gz`
   - etc.

2. Run the pipeline with Conda environments:
   ```
   snakemake --use-conda --cores [number_of_cores] --resources mem_mb=[memory_in_MB]
   ```

3. For a dry run to check the workflow:
   ```
   snakemake -n
   ```

4. To generate a workflow graph:
   ```
   snakemake --dag | dot -Tpng > dag.png
   ```

## Output

The pipeline generates:

1. **Viral Detection Results**:
   - `results/kraken2/{sample}_kraken2_report.txt`: Taxonomic classification of known viruses
   - `results/diamond/{sample}_viral_hits.txt`: Sensitive viral protein matches
   - `results/virsorter2/{sample}/final-viral-score.tsv`: Predicted viral contigs

2. **Viral Genomes and Annotation**:
   - `results/checkv/{sample}/high_quality_viruses.fasta`: High-quality viral genomes
   - `results/checkv/{sample}/medium_quality_viruses.fasta`: Medium-quality viral genomes
   - `results/pharokka/{sample}/pharokka_output.gbk`: Phage annotations
   - `results/vgas/{sample}/annotation_results.txt`: General viral annotations

3. **Viral Taxonomy**:
   - `results/vcontact2/{sample}/genome_by_genome_overview.csv`: Viral clustering results

4. **Summary Report**:
   - `results/final_report.html`: Comprehensive HTML report with visualizations

## Advanced Usage

1. **Resume a failed run**:
   ```
   snakemake --use-conda --cores [cores] --restart-times 3 --keep-going
   ```

2. **Run specific steps**:
   ```
   snakemake --use-conda --cores [cores] results/kraken2/sample1_kraken2_report.txt
   ```

3. **Execute on a cluster** (e.g., SLURM):
   ```
   snakemake --use-conda --cluster "sbatch --nodes=1 --ntasks=1 --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.time} --output=logs/slurm_%j.out" --jobs 100
   ```

## Customization

1. **Modify parameters**: Edit rule parameters in the Snakefile
2. **Add new tools**: Create new rules in the Snakefile and add corresponding Conda environment files
3. **Change output format**: Modify the report generation script (`scripts/generate_report.py`)

## Troubleshooting

1. **Memory issues**: Reduce `threads` and `mem_mb` parameters in resource-intensive rules
2. **Missing databases**: Ensure all database paths in `config.yaml` are correct
3. **Tool errors**: Check corresponding log files

## Pipeline Architecture

The pipeline follows a modular design with each tool encapsulated in its own rule. This architecture provides several benefits:

1. **Modularity**: Each step can be run independently
2. **Reproducibility**: Conda environments ensure consistent tool versions
3. **Scalability**: Resources are allocated per rule based on requirements
4. **Fault tolerance**: Failed steps can be resumed without rerunning successful steps

### Workflow Diagram

```
                                ┌─────────────────┐
                                │   Raw Reads     │
                                └────────┬────────┘
                                         │
                                         ▼
                                ┌─────────────────┐
                                │ Quality Control │
                                │     (FastP)     │
                                └────────┬────────┘
                                         │
                                         ▼
                                ┌─────────────────┐
                                │   Host Removal  │
                                │    (Bowtie2)    │
                                └────────┬────────┘
                                         │
                 ┌───────────────────────┼───────────────────────┐
                 │                       │                       │
                 ▼                       ▼                       ▼
        ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
        │  Known Viruses  │    │    Assembly     │    │ Sensitive Viral │
        │    (Kraken2)    │    │  (metaSPAdes)   │    │     Detection   │
        └─────────────────┘    └────────┬────────┘    │    (DIAMOND)    │
                                         │             └─────────────────┘
                                         │
                                         ▼
                                ┌─────────────────┐
                                │ Viral Detection │
                                │   (VirSorter2)  │
                                └────────┬────────┘
                                         │
                                         ▼
                                ┌─────────────────┐
                                │  Quality Check  │
                                │    (CheckV)     │
                                └────────┬────────┘
                                         │
                       ┌─────────────────┴─────────────────┐
                       │                                   │
                       ▼                                   ▼
            ┌─────────────────┐                  ┌─────────────────┐
            │ High-Quality    │                  │ Medium-Quality  │
            │     Viruses     │                  │     Viruses     │
            └────────┬────────┘                  └────────┬────────┘
                     │                                    │
                     ▼                                    ▼
            ┌─────────────────┐                  ┌─────────────────┐
            │ Phage Annotation│                  │ Viral Annotation│
            │    (Pharokka)   │                  │     (VGAS)      │
            └────────┬────────┘                  └────────┬────────┘
                     │                                    │
                     └────────────────┬──────────────────┘
                                      │
                                      ▼
                            ┌─────────────────┐
                            │ Viral Taxonomy  │
                            │   (vConTACT2)   │
                            └────────┬────────┘
                                     │
                                     ▼
                            ┌─────────────────┐
                            │   Final Report  │
                            └─────────────────┘
```

## Performance Optimization

This pipeline is optimized for both performance and resource efficiency:

1. **Parallelization**: 
   - Each rule uses multi-threading where supported
   - Multiple samples can be processed in parallel
   - Set `--cores` parameter to control overall parallelization

2. **Resource Management**:
   - Memory requirements are specified per rule
   - Temporary files are marked with `temp()` to free disk space
   - Intermediate files are preserved for debugging when needed

3. **Checkpointing**:
   - Snakemake automatically tracks rule dependencies
   - Failed runs can be resumed from the last successful step
   - Use `--keepgoing` to continue past failed rules

4. **Efficiency Considerations**:
   - Host removal before assembly significantly reduces computational burden
   - VirSorter2 runs before CheckV to filter for viral sequences
   - Viral quality assessment determines appropriate annotation tools

## Benchmarking

Approximate resource usage and runtime for a typical dataset (~10GB paired-end reads):

| Step | Runtime | CPU Usage | Memory Usage |
|------|---------|-----------|--------------|
| FastP QC | ~30 min | 8 cores | 8 GB |
| Host Removal | ~1-2 hours | 16 cores | 16 GB |
| Kraken2 | ~30 min | 16 cores | 64 GB |
| metaSPAdes | ~4-12 hours | 32 cores | 128 GB |
| DIAMOND | ~1-2 hours | 16 cores | 32 GB |
| VirSorter2 | ~2-4 hours | 16 cores | 32 GB |
| CheckV | ~1 hour | 16 cores | 32 GB |
| Pharokka | ~30 min | 16 cores | 32 GB |
| VGAS | ~30 min | 16 cores | 32 GB |
| vConTACT2 | ~1-2 hours | 16 cores | 64 GB |
| Report Generation | ~5 min | 1 core | 4 GB |

These estimates vary depending on dataset complexity and hardware. Adjust resource parameters in Snakefile as needed.