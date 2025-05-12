# scripts/compare_assemblies.py
# Script to compare assemblies from all reads vs. unclassified reads

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import tempfile
import subprocess
import logging
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('compare_assemblies')

def get_contig_stats(fasta_file):
    """Get statistics for contigs in a FASTA file"""
    contig_lengths = []
    gc_contents = []
    n_count = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        length = len(seq)
        contig_lengths.append(length)
        
        gc_content = (seq.count('G') + seq.count('C')) / length * 100
        gc_contents.append(gc_content)
        
        n_count += seq.count('N')
    
    if not contig_lengths:
        return {}
    
    # Calculate N50
    sorted_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    
    cumulative = 0
    n50 = 0
    n90 = 0
    
    for length in sorted_lengths:
        cumulative += length
        if not n50 and cumulative >= total_length * 0.5:
            n50 = length
        if not n90 and cumulative >= total_length * 0.9:
            n90 = length
    
    # Count contigs by size ranges
    size_ranges = {
        '>10kb': sum(1 for l in contig_lengths if l > 10000),
        '5-10kb': sum(1 for l in contig_lengths if 5000 <= l <= 10000),
        '1-5kb': sum(1 for l in contig_lengths if 1000 <= l < 5000),
        '<1kb': sum(1 for l in contig_lengths if l < 1000)
    }
    
    return {
        'total_contigs': len(contig_lengths),
        'total_length': total_length,
        'min_length': min(contig_lengths) if contig_lengths else 0,
        'max_length': max(contig_lengths) if contig_lengths else 0,
        'mean_length': np.mean(contig_lengths) if contig_lengths else 0,
        'median_length': np.median(contig_lengths) if contig_lengths else 0,
        'n50': n50,
        'n90': n90,
        'gc_mean': np.mean(gc_contents) if gc_contents else 0,
        'gc_std': np.std(gc_contents) if gc_contents else 0,
        'n_count': n_count,
        'size_ranges': size_ranges
    }

def blast_assemblies(all_contigs, unclassified_contigs, output_dir):
    """Compare assemblies using BLAST"""
    # Create BLAST database for all_contigs
    db_path = os.path.join(output_dir, "all_contigs_db")
    makeblastdb_cmd = [
        "makeblastdb",
        "-in", all_contigs,
        "-dbtype", "nucl",
        "-out", db_path
    ]
    subprocess.run(makeblastdb_cmd, check=True)
    
    # Run BLAST
    blast_output = os.path.join(output_dir, "blast_results.xml")
    blastn_cline = NcbiblastnCommandline(
        query=unclassified_contigs,
        db=db_path,
        outfmt=5,
        out=blast_output,
        evalue=1e-10,
        perc_identity=95,
        num_threads=snakemake.threads
    )
    stdout, stderr = blastn_cline()
    
    # Parse BLAST results
    blast_results = []
    with open(blast_output) as blast_file:
        blast_records = NCBIXML.parse(blast_file)
        for record in blast_records:
            query_id = record.query
            if record.alignments:
                for alignment in record.alignments:
                    subject_id = alignment.hit_id
                    for hsp in alignment.hsps:
                        if hsp.align_length / record.query_length >= 0.8:  # 80% coverage
                            blast_results.append({
                                'query_id': query_id,
                                'subject_id': subject_id,
                                'percent_identity': hsp.identities / hsp.align_length * 100,
                                'align_length': hsp.align_length,
                                'query_length': record.query_length,
                                'e_value': hsp.expect
                            })
                            break  # Only take best hit for each alignment
    
    return pd.DataFrame(blast_results) if blast_results else pd.DataFrame()

def calculate_assembly_overlap(all_contigs_file, unclass_contigs_file, blast_df):
    """Calculate overlap between assemblies based on BLAST results"""
    # Load contig IDs
    all_contig_ids = set(str(record.id) for record in SeqIO.parse(all_contigs_file, "fasta"))
    unclass_contig_ids = set(str(record.id) for record in SeqIO.parse(unclass_contigs_file, "fasta"))
    
    # Get matched contig IDs from BLAST results
    matched_unclass_ids = set() if blast_df.empty else set(blast_df['query_id'])
    matched_all_ids = set() if blast_df.empty else set(blast_df['subject_id'])
    
    # Calculate unique and shared contigs
    unique_all = all_contig_ids - matched_all_ids
    unique_unclass = unclass_contig_ids - matched_unclass_ids
    shared = len(matched_unclass_ids)  # Number of unclassified contigs that match all_contigs
    
    return {
        'total_all': len(all_contig_ids),
        'total_unclass': len(unclass_contig_ids),
        'unique_all': len(unique_all),
        'unique_unclass': len(unique_unclass),
        'shared': shared
    }

def plot_length_distribution(all_contigs_file, unclass_contigs_file, output_file):
    """Plot length distribution of contigs"""
    # Get contig lengths
    all_lengths = [len(record.seq) for record in SeqIO.parse(all_contigs_file, "fasta")]
    unclass_lengths = [len(record.seq) for record in SeqIO.parse(unclass_contigs_file, "fasta")]
    
    # Create histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    
    bins = np.logspace(np.log10(500), np.log10(max(max(all_lengths), max(unclass_lengths)) * 1.1), 50)
    
    ax.hist(all_lengths, bins=bins, alpha=0.7, label='All Reads Assembly')
    ax.hist(unclass_lengths, bins=bins, alpha=0.7, label='Unclassified Reads Assembly')
    
    ax.set_xscale('log')
    ax.set_xlabel('Contig Length (bp, log scale)')
    ax.set_ylabel('Number of Contigs')
    ax.set_title('Contig Length Distribution')
    ax.legend()
    ax.grid(True, which="both", ls="-", alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_venn_diagram(overlap_stats, output_file):
    """Create Venn diagram of assembly overlap"""
    plt.figure(figsize=(8, 6))
    v = venn2(
        subsets=(
            overlap_stats['unique_all'],
            overlap_stats['unique_unclass'],
            overlap_stats['shared']
        ),
        set_labels=('All Reads Assembly', 'Unclassified Reads Assembly')
    )
    
    # Add counts to labels
    v.get_label_by_id('10').set_text(f'Unique to All\n{overlap_stats["unique_all"]} contigs')
    v.get_label_by_id('01').set_text(f'Unique to Unclassified\n{overlap_stats["unique_unclass"]} contigs')
    v.get_label_by_id('11').set_text(f'Shared\n{overlap_stats["shared"]} contigs')
    
    plt.title('Contig Overlap Between Assembly Methods')
    plt.savefig(output_file)
    plt.close()

def create_comparison_table(all_stats, unclass_stats, overlap_stats, output_file):
    """Create comparison table between assemblies"""
    comparison = pd.DataFrame({
        'Metric': [
            'Total Contigs',
            'Assembly Size (bp)',
            'Longest Contig (bp)',
            'Shortest Contig (bp)',
            'Mean Contig Length (bp)',
            'Median Contig Length (bp)',
            'N50 (bp)',
            'N90 (bp)',
            'Mean GC Content (%)',
            'Contigs >10kb',
            'Contigs 5-10kb',
            'Contigs 1-5kb',
            'Contigs <1kb',
            'Unique Contigs',
            'Shared Contigs'
        ],
        'All_Reads': [
            all_stats['total_contigs'],
            all_stats['total_length'],
            all_stats['max_length'],
            all_stats['min_length'],
            round(all_stats['mean_length']),
            round(all_stats['median_length']),
            all_stats['n50'],
            all_stats['n90'],
            round(all_stats['gc_mean'], 2),
            all_stats['size_ranges']['>10kb'],
            all_stats['size_ranges']['5-10kb'],
            all_stats['size_ranges']['1-5kb'],
            all_stats['size_ranges']['<1kb'],
            overlap_stats['unique_all'],
            overlap_stats['shared']
        ],
        'Unclassified_Reads': [
            unclass_stats['total_contigs'],
            unclass_stats['total_length'],
            unclass_stats['max_length'],
            unclass_stats['min_length'],
            round(unclass_stats['mean_length']),
            round(unclass_stats['median_length']),
            unclass_stats['n50'],
            unclass_stats['n90'],
            round(unclass_stats['gc_mean'], 2),
            unclass_stats['size_ranges']['>10kb'],
            unclass_stats['size_ranges']['5-10kb'],
            unclass_stats['size_ranges']['1-5kb'],
            unclass_stats['size_ranges']['<1kb'],
            overlap_stats['unique_unclass'],
            overlap_stats['shared']
        ]
    })
    
    # Calculate difference and percent change
    comparison['Difference'] = comparison['All_Reads'] - comparison['Unclassified_Reads']
    
    # Add percent change where it makes sense
    numeric_metrics = [
        'Total Contigs', 'Assembly Size (bp)', 'Longest Contig (bp)', 
        'Mean Contig Length (bp)', 'Median Contig Length (bp)', 'N50 (bp)',
        'N90 (bp)', 'Contigs >10kb', 'Contigs 5-10kb', 'Contigs 1-5kb',
        'Contigs <1kb', 'Unique Contigs'
    ]
    
    comparison['Percent_Change'] = None
    for metric in numeric_metrics:
        if comparison.loc[comparison['Metric'] == metric, 'Unclassified_Reads'].values[0] != 0:
            pct_change = (comparison.loc[comparison['Metric'] == metric, 'Difference'].values[0] / 
                         comparison.loc[comparison['Metric'] == metric, 'Unclassified_Reads'].values[0] * 100)
            comparison.loc[comparison['Metric'] == metric, 'Percent_Change'] = f"{pct_change:.1f}%"
    
    comparison.to_csv(output_file, sep='\t', index=False)
    return comparison

def main():
    """Main function to compare assemblies"""
    # Get input file paths from Snakemake
    all_contigs_file = snakemake.input.all_contigs
    unclass_contigs_file = snakemake.input.unclassified_contigs
    
    # Output files
    comparison_stats_file = snakemake.output.comparison_stats
    venn_diagram_file = snakemake.output.venn_diagram
    length_dist_file = snakemake.output.length_dist
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(comparison_stats_file)
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info("Starting assembly comparison")
    
    # Get contig statistics
    logger.info("Calculating assembly statistics")
    all_stats = get_contig_stats(all_contigs_file)
    unclass_stats = get_contig_stats(unclass_contigs_file)
    
    if not all_stats or not unclass_stats:
        logger.error("One or both assemblies are empty. Cannot compare.")
        # Create empty output files
        pd.DataFrame().to_csv(comparison_stats_file, sep='\t', index=False)
        plt.figure().savefig(venn_diagram_file)
        plt.figure().savefig(length_dist_file)
        return
    
    # Compare assemblies with BLAST
    logger.info("Running BLAST to compare assemblies")
    blast_df = blast_assemblies(all_contigs_file, unclass_contigs_file, output_dir)
    
    # Calculate overlap
    logger.info("Calculating assembly overlap")
    overlap_stats = calculate_assembly_overlap(all_contigs_file, unclass_contigs_file, blast_df)
    
    # Create comparison table
    logger.info("Creating comparison table")
    create_comparison_table(all_stats, unclass_stats, overlap_stats, comparison_stats_file)
    
    # Plot Venn diagram
    logger.info("Creating Venn diagram")
    plot_venn_diagram(overlap_stats, venn_diagram_file)
    
    # Plot length distribution
    logger.info("Creating length distribution plot")
    plot_length_distribution(all_contigs_file, unclass_contigs_file, length_dist_file)
    
    logger.info("Assembly comparison completed")

if __name__ == "__main__":
    main()


