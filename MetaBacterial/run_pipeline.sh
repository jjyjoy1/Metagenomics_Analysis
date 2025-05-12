#!/bin/bash

# Script to run the metagenomics pipeline
# Usage: ./run_pipeline.sh [jobs] [restart_times]

# Default parameters
JOBS=${1:-1}  # Number of concurrent jobs, default 1
RESTART=${2:-3}  # Number of times to restart a job on failure, default 3

# Create log directory
mkdir -p logs

# Run snakemake
echo "Starting pipeline with $JOBS concurrent jobs and $RESTART restart attempts"

snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --jobs $JOBS \
    --use-conda \
    --restart-times $RESTART \
    --keep-going \
    --rerun-incomplete \
    --latency-wait 60 \
    --cluster "sbatch --nodes=1 --ntasks=1 --cpus-per-task={threads} --mem={resources.mem_mb}M --time=24:00:00 --output=logs/slurm-%j.out" \
    "$@"

# If you're not running on a cluster with SLURM, comment out the --cluster line
# and uncomment the following line:
# snakemake --snakefile Snakefile --configfile config.yaml --jobs $JOBS --use-conda --restart-times $RESTART --keep-going --rerun-incomplete "$@"

echo "Pipeline finished"
