# scripts/create_gene2genome.py
# Script to create gene-to-genome mapping file for vConTACT2

import argparse
from Bio import SeqIO
import pandas as pd
import re
import sys
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('create_gene2genome')

def create_gene_to_genome_mapping(proteins_file, output_file):
    """
    Create gene-to-genome mapping file required by vConTACT2
    
    Parameters:
    -----------
    proteins_file : str
        Path to the predicted proteins FASTA file from Prodigal
    output_file : str
        Path to the output CSV file for gene-to-genome mapping
    """
    logger.info(f"Reading protein sequences from {proteins_file}")
    try:
        # Parse protein file and extract gene and genome IDs
        gene_to_genome = []
        
        for record in SeqIO.parse(proteins_file, "fasta"):
            # Parse Prodigal header to get gene info
            # Example header: >contig123_1 # 1 # 300 # 1 # ID=1_1;partial=00;start_type=ATG;...
            gene_id = record.id
            genome_id = gene_id.split('_')[0]
            
            gene_to_genome.append({
                'protein_id': gene_id,
                'contig_id': genome_id
            })
            
        logger.info(f"Found {len(gene_to_genome)} genes from {len(set([x['contig_id'] for x in gene_to_genome]))} genomes")
        
        # Create DataFrame and write to CSV
        df = pd.DataFrame(gene_to_genome)
        df.to_csv(output_file, index=False)
        logger.info(f"Gene-to-genome mapping written to {output_file}")
        
    except Exception as e:
        logger.error(f"Error creating gene-to-genome mapping: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create gene-to-genome mapping for vConTACT2')
    parser.add_argument('--proteins', required=True, help='Path to the predicted proteins FASTA file')
    parser.add_argument('--output', required=True, help='Path to the output CSV file')
    
    if 'snakemake' in globals():
        # Running from Snakemake
        proteins_file = snakemake.input.protein_file
        output_file = snakemake.output.gene_2_genome
    else:
        # Running from command line
        args = parser.parse_args()
        proteins_file = args.proteins
        output_file = args.output
    
    create_gene_to_genome_mapping(proteins_file, output_file)


