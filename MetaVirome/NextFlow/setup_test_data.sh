#!/bin/bash
# Script to set up test environment for MetaVirome Pipeline

set -e

# Create test data directory structure
mkdir -p test-data/reads
mkdir -p test-data/reference
mkdir -p test-data/databases/kraken2_mini_db
mkdir -p test-data/databases/diamond_mini_db

echo "Setting up test environment..."

# Download small test dataset (SRA)
echo "Downloading test reads..."
wget -O test-data/reads/test_R1.fastq.gz https://zenodo.org/record/3951383/files/noro_sample1_R1.fastq.gz
wget -O test-data/reads/test_R2.fastq.gz https://zenodo.org/record/3951383/files/noro_sample1_R2.fastq.gz

# Download small host genome (E. coli as an example)
echo "Downloading host genome..."
wget -O test-data/reference/host_genome.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip test-data/reference/host_genome.fasta.gz

# Create mini Kraken2 database
echo "Setting up mini Kraken2 database..."
mkdir -p tmp_viral_db
cd tmp_viral_db
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/839/085/GCF_000839085.1_ViralProj14032/GCF_000839085.1_ViralProj14032_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/854/985/GCF_000854985.1_ViralProj15201/GCF_000854985.1_ViralProj15201_genomic.fna.gz
gunzip *.gz
cat *.fna > viral_mini.fna
kraken2-build --add-to-library viral_mini.fna --db ../test-data/databases/kraken2_mini_db
kraken2-build --build --threads 2 --db ../test-data/databases/kraken2_mini_db
kraken2-build --clean --db ../test-data/databases/kraken2_mini_db
cd ..
rm -rf tmp_viral_db

# Create mini DIAMOND database
echo "Setting up mini DIAMOND database..."
mkdir -p tmp_viral_prot
cd tmp_viral_prot
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/839/085/GCF_000839085.1_ViralProj14032/GCF_000839085.1_ViralProj14032_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/854/985/GCF_000854985.1_ViralProj15201/GCF_000854985.1_ViralProj15201_protein.faa.gz
gunzip *.gz
cat *.faa > viral_mini.faa
diamond makedb --in viral_mini.faa --db ../test-data/databases/diamond_mini_db/viral_proteins
cd ..
rm -rf tmp_viral_prot

echo "Test environment setup complete!"
echo "To run the test pipeline:"
echo "nextflow run metavirome.nf -profile test"


