FROM --platform=linux/arm64 arm64v8/debian:bullseye-slim

# Install dependencies
RUN apt-get update && \
    apt-get install -y wget build-essential zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install pigz
RUN apt-get update && apt-get install -y pigz

# Install fastp from source
WORKDIR /tmp
RUN wget https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.4.tar.gz && \
    tar -xzf v0.23.4.tar.gz && \
    cd fastp-0.23.4 && \
    make && \
    make install && \
    cd .. && \
    rm -rf fastp-0.23.4 v0.23.4.tar.gz

# Set working directory
WORKDIR /data

# Set entrypoint
ENTRYPOINT ["fastp"]

