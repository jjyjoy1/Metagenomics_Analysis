# scripts/extract_viral_contigs.py
# Script to extract viral contigs based on VirSorter2 results

import pandas as pd
from Bio import SeqIO
import sys
import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('extract_viral_contigs')

def extract_viral_contigs(contigs_file, vs2_results_file, output_file, score_threshold=0.5):
    """
    Extract viral contigs from assembly based on VirSorter2 results
    
    Parameters:
    -----------
    contigs_file : str
        Path to the input contigs FASTA file
    vs2_results_file : str
        Path to the VirSorter2 results file (final-viral-score.tsv)
    output_file : str
        Path to the output FASTA file for viral contigs
    score_threshold : float, default=0.5
        Minimum score to consider a contig as viral
    """
    # Read VirSorter2 results
    logger.info(f"Reading VirSorter2 results from {vs2_results_file}")
    try:
        vs2_results = pd.read_csv(vs2_results_file, sep='\t')
        logger.info(f"Found {len(vs2_results)} contigs in VirSorter2 results")
    except Exception as e:
        logger.error(f"Error reading VirSorter2 results: {e}")
        sys.exit(1)
    
    # Filter by score threshold
    viral_contigs = vs2_results[vs2_results['max_score'] >= score_threshold]
    logger.info(f"Found {len(viral_contigs)} contigs with score >= {score_threshold}")
    
    # Create a set of viral contig IDs for faster lookup
    viral_contig_ids = set(viral_contigs['seqname'])
    
    # Extract viral contigs from the FASTA file
    logger.info(f"Extracting viral contigs from {contigs_file}")
    try:
        count = 0
        with open(output_file, 'w') as out_handle:
            for record in SeqIO.parse(contigs_file, "fasta"):
                if record.id in viral_contig_ids:
                    SeqIO.write(record, out_handle, "fasta")
                    count += 1
        logger.info(f"Wrote {count} viral contigs to {output_file}")
    except Exception as e:
        logger.error(f"Error extracting viral contigs: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Get file paths from Snakemake
    contigs_file = snakemake.input.contigs
    vs2_results_file = snakemake.input.vs_results
    output_file = snakemake.output.viral_contigs
    
    # Extract viral contigs
    extract_viral_contigs(contigs_file, vs2_results_file, output_file)

