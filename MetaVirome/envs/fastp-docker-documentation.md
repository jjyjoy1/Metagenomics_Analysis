# fastp Docker Container: Documentation

This document provides complete instructions for building and using a Docker container for the fastp tool, compatible with both Linux and Apple Silicon (M-series chips) platforms.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Container Information](#container-information)
- [Building the Docker Image](#building-the-docker-image)
  - [For Apple Silicon (M-series chips)](#for-apple-silicon-m-series-chips)
  - [For Linux/Intel platforms](#for-linuxintel-platforms)
- [Using the fastp Container](#using-the-fastp-container)
  - [Basic Commands](#basic-commands)
  - [Common Use Cases](#common-use-cases)
  - [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Additional Resources](#additional-resources)

## Prerequisites

Before you begin, ensure you have the following installed:
- Docker
  - On Mac: [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop/)
  - On Linux: Docker Engine (via package manager)
- Git (optional, for version control)

## Container Information

This Docker container includes:
- fastp version 0.23.4 - A tool designed to provide fast all-in-one preprocessing for FastQ files
- pigz version 2.6 - Parallel implementation of gzip for modern multi-processor systems

## Building the Docker Image

### For Apple Silicon (M-series chips)

1. Create a new directory for your Docker project:
   ```bash
   mkdir fastp-docker
   cd fastp-docker
   ```

2. Create an `environment.yml` file with the following content:
   ```yaml
   name: fastp
   channels:
     - bioconda
     - conda-forge
     - defaults
   dependencies:
     - fastp=0.23.4
     - pigz=2.6
   ```

3. Create a `Dockerfile` with the following content:
   ```dockerfile
   # Use ARM64 base image
   FROM --platform=linux/arm64 debian:bullseye-slim

   # Install dependencies
   RUN apt-get update && \
       apt-get install -y wget curl bzip2 ca-certificates git procps && \
       apt-get clean && \
       rm -rf /var/lib/apt/lists/*

   # Install miniconda for ARM64
   RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh && \
       bash miniconda.sh -b -p /opt/conda && \
       rm miniconda.sh

   # Set PATH
   ENV PATH=/opt/conda/bin:$PATH

   # Set working directory
   WORKDIR /app

   # Copy environment file
   COPY environment.yml .

   # Create environment
   RUN conda config --set channel_priority strict && \
       conda env create -f environment.yml && \
       conda clean -afy

   # Set path to the environment
   ENV PATH=/opt/conda/envs/fastp/bin:$PATH

   # Set working directory for data processing
   WORKDIR /data

   # Set entrypoint
   ENTRYPOINT ["fastp"]
   ```

4. Build the Docker image:
   ```bash
   docker build --platform linux/arm64 -t fastp:0.23.4 .
   ```

### For Linux/Intel platforms

1. Create a new directory for your Docker project (same as above)

2. Create the same `environment.yml` file as above

3. Create a `Dockerfile` with the following content:
   ```dockerfile
   # Use a minimal base image with conda pre-installed
   FROM continuumio/miniconda3:latest

   # Set working directory
   WORKDIR /app

   # Copy conda environment file
   COPY environment.yml .

   # Create conda environment based on the environment.yml file
   RUN conda env create -f environment.yml

   # Add conda environment to PATH
   SHELL ["/bin/bash", "-c"]
   RUN echo "source activate fastp" > ~/.bashrc
   ENV PATH=/opt/conda/envs/fastp/bin:$PATH

   # Set default command to run fastp
   ENTRYPOINT ["fastp"]
   ```

4. Build the Docker image:
   ```bash
   docker build -t fastp:0.23.4 .
   ```

## Using the fastp Container

### Basic Commands

1. Display help information:
   ```bash
   docker run --rm fastp:0.23.4 --help
   ```

2. Check the version:
   ```bash
   docker run --rm fastp:0.23.4 --version
   ```

### Common Use Cases

1. **Process a single-end FASTQ file**:
   ```bash
   docker run --rm -v $(pwd):/data fastp:0.23.4 \
     -i /data/input.fastq.gz \
     -o /data/output.fastq.gz
   ```

2. **Process paired-end FASTQ files**:
   ```bash
   docker run --rm -v $(pwd):/data fastp:0.23.4 \
     -i /data/input_R1.fastq.gz \
     -I /data/input_R2.fastq.gz \
     -o /data/output_R1.fastq.gz \
     -O /data/output_R2.fastq.gz
   ```

3. **Generate HTML and JSON reports**:
   ```bash
   docker run --rm -v $(pwd):/data fastp:0.23.4 \
     -i /data/input_R1.fastq.gz \
     -I /data/input_R2.fastq.gz \
     -o /data/output_R1.fastq.gz \
     -O /data/output_R2.fastq.gz \
     --html /data/report.html \
     --json /data/report.json
   ```

### Advanced Usage

1. **Quality filtering and trimming**:
   ```bash
   docker run --rm -v $(pwd):/data fastp:0.23.4 \
     -i /data/input_R1.fastq.gz \
     -I /data/input_R2.fastq.gz \
     -o /data/output_R1.fastq.gz \
     -O /data/output_R2.fastq.gz \
     --qualified_quality_phred 20 \
     --unqualified_percent_limit 40 \
     --length_required 50
   ```

2. **Adapter trimming**:
   ```bash
   docker run --rm -v $(pwd):/data fastp:0.23.4 \
     -i /data/input_R1.fastq.gz \
     -I /data/input_R2.fastq.gz \
     -o /data/output_R1.fastq.gz \
     -O /data/output_R2.fastq.gz \
     --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
     --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
   ```

3. **Multi-threading for better performance**:
   ```bash
   docker run --rm -v $(pwd):/data fastp:0.23.4 \
     -i /data/input_R1.fastq.gz \
     -I /data/input_R2.fastq.gz \
     -o /data/output_R1.fastq.gz \
     -O /data/output_R2.fastq.gz \
     --thread 8
   ```

4. **Detailed QC with base correction**:
   ```bash
   docker run --rm -v $(pwd):/data fastp:0.23.4 \
     -i /data/input_R1.fastq.gz \
     -I /data/input_R2.fastq.gz \
     -o /data/output_R1.fastq.gz \
     -O /data/output_R2.fastq.gz \
     --cut_front \
     --cut_tail \
     --cut_window_size 4 \
     --cut_mean_quality 20 \
     --correction \
     --html /data/report.html \
     --json /data/report.json
   ```

## Troubleshooting

### "Cannot connect to the Docker daemon" error
- Ensure Docker Desktop is running
- For Mac users, check that Docker Desktop is properly initialized
- Try running `docker info` to verify the Docker daemon is running

### Architecture errors
- On Apple Silicon (M-series) Macs, always use `--platform linux/arm64` when building
- If getting "exec format error" when running a container, you may have built for the wrong architecture

### Pull errors (401 Unauthorized)
- Try logging in with `docker login`
- Check your network connection
- Verify your Docker Hub account has access to the required images

### Build failures
- Check the Docker build log for specific errors
- Ensure all dependencies are available
- For conda environment issues, try using a different approach (direct binary installation)

## Additional Resources

- [fastp GitHub Repository](https://github.com/OpenGene/fastp)
- [Docker Documentation](https://docs.docker.com/)
- [Bioconda Documentation](https://bioconda.github.io/)
- [Conda Documentation](https://docs.conda.io/en/latest/)

---

Created for use with fastp v0.23.4 and pigz v2.6. Last updated: May 2025.
