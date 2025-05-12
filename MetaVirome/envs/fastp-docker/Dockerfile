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


