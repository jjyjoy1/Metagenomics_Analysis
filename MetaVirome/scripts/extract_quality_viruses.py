# scripts/extract_quality_viruses.py
# Script to extract high-quality and medium-quality viral contigs based on CheckV results

import pandas as pd
from Bio import SeqIO
import sys
import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('extract_quality_viruses')

def extract_quality_viruses(viral_contigs_file, checkv_quality_file, hq_output_file, mq_output_file):
    """
    Extract high-quality and medium-quality viral contigs based on CheckV quality assessment
    
    Parameters:
    -----------
    viral_contigs_file : str
        Path to the input viral contigs FASTA file
    checkv_quality_file : str
        Path to the CheckV quality summary file
    hq_output_file : str
        Path to the output FASTA file for high-quality viral contigs
    mq_output_file : str
        Path to the output FASTA file for medium-quality viral contigs
    """
    # Read CheckV quality results
    logger.info(f"Reading CheckV quality results from {checkv_quality_file}")
    try:
        checkv_results = pd.read_csv(checkv_quality_file, sep='\t')
        logger.info(f"Found {len(checkv_results)} contigs in CheckV results")
    except Exception as e:
        logger.error(f"Error reading CheckV results: {e}")
        sys.exit(1)
    
    # Filter by quality levels
    high_quality = checkv_results[checkv_results['checkv_quality'] == 'High-quality']
    medium_quality = checkv_results[checkv_results['checkv_quality'] == 'Medium-quality']
    
    logger.info(f"Found {len(high_quality)} high-quality viral contigs")
    logger.info(f"Found {len(medium_quality)} medium-quality viral contigs")
    
    # Create sets of contig IDs for faster lookup
    hq_ids = set(high_quality['contig_id'])
    mq_ids = set(medium_quality['contig_id'])
    
    # Extract contigs by quality
    logger.info(f"Extracting quality-filtered viral contigs from {viral_contigs_file}")
    try:
        # Load all viral contigs into memory
        all_contigs = {record.id: record for record in SeqIO.parse(viral_contigs_file, "fasta")}
        
        # Write high-quality viral contigs
        hq_count = 0
        with open(hq_output_file, 'w') as out_handle:
            for contig_id in hq_ids:
                if contig_id in all_contigs:
                    SeqIO.write(all_contigs[contig_id], out_handle, "fasta")
                    hq_count += 1
        logger.info(f"Wrote {hq_count} high-quality viral contigs to {hq_output_file}")
        
        # Write medium-quality viral contigs
        mq_count = 0
        with open(mq_output_file, 'w') as out_handle:
            for contig_id in mq_ids:
                if contig_id in all_contigs:
                    SeqIO.write(all_contigs[contig_id], out_handle, "fasta")
                    mq_count += 1
        logger.info(f"Wrote {mq_count} medium-quality viral contigs to {mq_output_file}")
    except Exception as e:
        logger.error(f"Error extracting quality-filtered viral contigs: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Get file paths from Snakemake
    viral_contigs_file = snakemake.input.viral_contigs
    checkv_quality_file = snakemake.input.checkv_quality
    hq_output_file = snakemake.output.hq_viruses
    mq_output_file = snakemake.output.mq_viruses
    
    # Extract quality-filtered viral contigs
    extract_quality_viruses(viral_contigs_file, checkv_quality_file, hq_output_file, mq_output_file)

