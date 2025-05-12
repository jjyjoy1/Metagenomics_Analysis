#!/usr/bin/env python3
"""
Generate a summary report for the metagenomic analysis pipeline.
This script collects results from various tools and creates a comprehensive HTML report.

Usage:
    python generate_report.py -c config.yaml -o report.html
"""

import argparse
import os
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template
import numpy as np
from datetime import datetime
import glob
import json

def parse_args():
    parser = argparse.ArgumentParser(description='Generate a summary report for the metagenomic analysis pipeline')
    parser.add_argument('-c', '--config', required=True, help='Path to the config.yaml file')
    parser.add_argument('-o', '--output', default='report.html', help='Output HTML report file')
    return parser.parse_args()

def load_config(config_file):
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config

def get_fastp_stats(config):
    """Extract quality stats from fastp JSON files"""
    fastp_dir = os.path.join(config['outdir'], '01_fastp')
    stats = {}
    
    for sample in config['samples']:
        json_file = os.path.join(fastp_dir, f"{sample}.fastp.json")
        if os.path.exists(json_file):
            with open(json_file, 'r') as f:
                data = json.load(f)
                
            stats[sample] = {
                'raw_reads': data['summary']['before_filtering']['total_reads'],
                'clean_reads': data['summary']['after_filtering']['total_reads'],
                'raw_bases': data['summary']['before_filtering']['total_bases'],
                'clean_bases': data['summary']['after_filtering']['total_bases'],
                'q30_rate': data['summary']['after_filtering']['q30_rate'],
                'gc_content': data['summary']['after_filtering']['gc_content']
            }
    
    return pd.DataFrame(stats).T if stats else pd.DataFrame()

def get_host_removal_stats(config):
    """Calculate host removal statistics"""
    hostfree_dir = os.path.join(config['outdir'], '02_hostfree')
    fastp_dir = os.path.join(config['outdir'], '01_fastp')
    stats = {}
    
    for sample in config['samples']:
        trimmed_r1 = os.path.join(fastp_dir, f"{sample}.trimmed.R1.fastq.gz")
        trimmed_r2 = os.path.join(fastp_dir, f"{sample}.trimmed.R2.fastq.gz")
        hostfree_r1 = os.path.join(hostfree_dir, f"{sample}.hostfree.R1.fastq.gz")
        hostfree_r2 = os.path.join(hostfree_dir, f"{sample}.hostfree.R2.fastq.gz")
        
        # Get file sizes to estimate read counts (faster than actually counting)
        if os.path.exists(trimmed_r1) and os.path.exists(hostfree_r1):
            trimmed_size = os.path.getsize(trimmed_r1) + os.path.getsize(trimmed_r2)
            hostfree_size = os.path.getsize(hostfree_r1) + os.path.getsize(hostfree_r2)
            
            # Approximate percentage of non-host reads
            non_host_pct = (hostfree_size / trimmed_size) * 100 if trimmed_size > 0 else 0
            host_pct = 100 - non_host_pct
            
            stats[sample] = {
                'host_pct': host_pct,
                'non_host_pct': non_host_pct
            }
    
    return pd.DataFrame(stats).T if stats else pd.DataFrame()

def get_kraken_top_taxa(config, level='S', top_n=10):
    """Extract top taxa from Bracken results"""
    bracken_dir = os.path.join(config['outdir'], '04_bracken')
    all_data = []
    
    level_name = {'D': 'Domain', 'P': 'Phylum', 'C': 'Class', 
                 'O': 'Order', 'F': 'Family', 'G': 'Genus', 'S': 'Species'}
    
    for sample in config['samples']:
        if level == 'S':
            bracken_file = os.path.join(bracken_dir, f"{sample}.bracken.species.txt")
        elif level == 'G':
            bracken_file = os.path.join(bracken_dir, f"{sample}.bracken.genus.txt")
        else:
            continue
            
        if os.path.exists(bracken_file):
            try:
                df = pd.read_csv(bracken_file, sep='\t')
                if not df.empty:
                    # Normalize to percentage
                    df['fraction'] = df['fraction_total_reads'] * 100
                    # Get top N taxa
                    top_taxa = df.sort_values('fraction', ascending=False).head(top_n)
                    
                    for _, row in top_taxa.iterrows():
                        name = row['name'].split()
                        if level == 'S' and len(name) >= 2:
                            # Format species name
                            taxa_name = f"{name[0]} {name[1]}"
                        else:
                            taxa_name = row['name']
                            
                        all_data.append({
                            'sample': sample,
                            'taxa': taxa_name,
                            'fraction': row['fraction']
                        })
            except Exception as e:
                print(f"Error processing {bracken_file}: {e}")
                
    return pd.DataFrame(all_data) if all_data else pd.DataFrame()

def get_assembly_stats(config):
    """Get assembly statistics"""
    assembly_dir = os.path.join(config['outdir'], '06_megahit')
    stats = {}
    
    for sample in config['samples']:
        contigs_file = os.path.join(assembly_dir, f"{sample}_final.contigs.fa")
        
        if os.path.exists(contigs_file):
            # Count contigs and total length
            total_len = 0
            contig_lengths = []
            
            with open(contigs_file, 'r') as f:
                current_len = 0
                for line in f:
                    if line.startswith('>'):
                        if current_len > 0:
                            contig_lengths.append(current_len)
                            total_len += current_len
                            current_len = 0
                    else:
                        current_len += len(line.strip())
                        
                # Add the last contig
                if current_len > 0:
                    contig_lengths.append(current_len)
                    total_len += current_len
                    
            if contig_lengths:
                contig_lengths.sort(reverse=True)
                n50 = calculate_n50(contig_lengths)
                
                stats[sample] = {
                    'contig_count': len(contig_lengths),
                    'total_length': total_len,
                    'n50': n50,
                    'longest_contig': max(contig_lengths) if contig_lengths else 0,
                    'avg_contig_len': total_len / len(contig_lengths) if len(contig_lengths) > 0 else 0
                }
    
    return pd.DataFrame(stats).T if stats else pd.DataFrame()

def calculate_n50(contig_lengths):
    """Calculate N50 from a list of contig lengths"""
    total_len = sum(contig_lengths)
    target = total_len * 0.5
    running_sum = 0
    
    for length in contig_lengths:
        running_sum += length
        if running_sum >= target:
            return length
    
    return 0

def get_binning_stats(config):
    """Get binning statistics from CheckM and GTDB-Tk"""
    checkm_file = os.path.join(config['outdir'], '09_checkm', 'quality_report.tsv')
    gtdbtk_file = os.path.join(config['outdir'], '10_gtdbtk', 'gtdbtk.summary.tsv')
    
    # Process CheckM data
    checkm_data = {}
    if os.path.exists(checkm_file):
        try:
            checkm_df = pd.read_csv(checkm_file, sep='\t')
            for _, row in checkm_df.iterrows():
                bin_id = row['Bin Id']
                checkm_data[bin_id] = {
                    'Completeness': row['Completeness'],
                    'Contamination': row['Contamination'],
                    'Strain_heterogeneity': row.get('Strain heterogeneity', 0),
                    'Genome_size': row.get('Genome size', 0),
                    'GC': row.get('GC', 0),
                    'N50': row.get('N50 (scaffolds)', 0)
                }
        except Exception as e:
            print(f"Error processing CheckM file: {e}")
    
    # Process GTDB-Tk data
    gtdbtk_data = {}
    if os.path.exists(gtdbtk_file):
        try:
            gtdbtk_df = pd.read_csv(gtdbtk_file, sep='\t')
            for _, row in gtdbtk_df.iterrows():
                bin_id = row['user_genome']
                
                # Extract taxonomy
                if 'classification' in row:
                    tax = row['classification'].split(';')
                    domain = tax[0][3:] if len(tax) > 0 and tax[0].startswith('d__') else ''
                    phylum = tax[1][3:] if len(tax) > 1 and tax[1].startswith('p__') else ''
                    genus = tax[5][3:] if len(tax) > 5 and tax[5].startswith('g__') else ''
                    species = tax[6][3:] if len(tax) > 6 and tax[6].startswith('s__') else ''
                else:
                    domain, phylum, genus, species = '', '', '', ''
                
                gtdbtk_data[bin_id] = {
                    'Domain': domain,
                    'Phylum': phylum,
                    'Genus': genus,
                    'Species': species,
                    'Taxonomy': row.get('classification', '')
                }
        except Exception as e:
            print(f"Error processing GTDB-Tk file: {e}")
    
    # Combine data
    combined_data = []
    all_bins = set(list(checkm_data.keys()) + list(gtdbtk_data.keys()))
    
    for bin_id in all_bins:
        bin_data = {'bin_id': bin_id}
        
        # Add CheckM data
        if bin_id in checkm_data:
            bin_data.update(checkm_data[bin_id])
        
        # Add GTDB-Tk data
        if bin_id in gtdbtk_data:
            bin_data.update(gtdbtk_data[bin_id])
            
        combined_data.append(bin_data)
    
    return pd.DataFrame(combined_data) if combined_data else pd.DataFrame()

def get_humann_pathways(config, top_n=10):
    """Get top metabolic pathways from HUMAnN3 results"""
    humann_dir = os.path.join(config['outdir'], '05_humann')
    all_pathways = []
    
    for sample in config['samples']:
        pathway_file = os.path.join(humann_dir, f"{sample}_pathabundance.tsv")
        
        if os.path.exists(pathway_file):
            try:
                df = pd.read_csv(pathway_file, sep='\t', skiprows=[0])
                
                # Filter for community-level pathways
                community_paths = df[~df['# Pathway'].str.contains('|')].copy()
                
                # Remove "UNMAPPED" and "UNINTEGRATED"
                community_paths = community_paths[~community_paths['# Pathway'].isin(['UNMAPPED', 'UNINTEGRATED'])]
                
                # Normalize values and sort
                if not community_paths.empty:
                    path_col = community_paths.columns[0]
                    value_col = community_paths.columns[1]
                    
                    # Get top pathways
                    top_paths = community_paths.sort_values(value_col, ascending=False).head(top_n)
                    
                    for _, row in top_paths.iterrows():
                        pathway = row[path_col]
                        # Simplify pathway name
                        simple_name = pathway.split(':')[-1] if ':' in pathway else pathway
                        
                        all_pathways.append({
                            'sample': sample,
                            'pathway': simple_name,
                            'abundance': row[value_col]
                        })
            except Exception as e:
                print(f"Error processing {pathway_file}: {e}")
    
    return pd.DataFrame(all_pathways) if all_pathways else pd.DataFrame()

def generate_plots(data_dict, output_dir='plots'):
    """Generate plots for the report"""
    os.makedirs(output_dir, exist_ok=True)
    plot_files = {}
    
    # Set the style
    sns.set(style="whitegrid")
    plt.rcParams.update({'font.size': 12})
    
    # 1. Read QC plot
    if 'fastp_stats' in data_dict and not data_dict['fastp_stats'].empty:
        df = data_dict['fastp_stats']
        
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(df.index))
        width = 0.35
        
        ax.bar(x - width/2, df['raw_reads']/1000000, width, label='Raw Reads')
        ax.bar(x + width/2, df['clean_reads']/1000000, width, label='Clean Reads')
        
        ax.set_xlabel('Sample')
        ax.set_ylabel('Millions of Reads')
        ax.set_title('Read Counts Before and After QC')
        ax.set_xticks(x)
        ax.set_xticklabels(df.index, rotation=45, ha='right')
        ax.legend()
        
        fig.tight_layout()
        plot_file = os.path.join(output_dir, 'qc_stats.png')
        plt.savefig(plot_file, dpi=300)
        plt.close()
        plot_files['qc_stats'] = plot_file
    
    # 2. Host removal plot
    if 'host_stats' in data_dict and not data_dict['host_stats'].empty:
        df = data_dict['host_stats']
        
        fig, ax = plt.subplots(figsize=(10, 6))
        df.plot(kind='bar', stacked=True, ax=ax)
        
        ax.set_xlabel('Sample')
        ax.set_ylabel('Percentage')
        ax.set_title('Host vs. Non-host Read Percentage')
        ax.legend(title='Read Type')
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'host_removal.png')
        plt.savefig(plot_file, dpi=300)
        plt.close()
        plot_files['host_removal'] = plot_file
    
    # 3. Taxonomic composition plot
    if 'taxonomy_data' in data_dict and not data_dict['taxonomy_data'].empty:
        df = data_dict['taxonomy_data']
        
        # Create a pivot table for plotting
        pivot_df = df.pivot(index='taxa', columns='sample', values='fraction').fillna(0)
        
        fig, ax = plt.subplots(figsize=(12, 8))
        pivot_df.plot(kind='barh', stacked=True, ax=ax, cmap='tab20')
        
        ax.set_xlabel('Relative Abundance (%)')
        ax.set_title('Top Species Composition')
        ax.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'taxonomy.png')
        plt.savefig(plot_file, dpi=300)
        plt.close()
        plot_files['taxonomy'] = plot_file
    
    # 4. Assembly stats plot
    if 'assembly_stats' in data_dict and not data_dict['assembly_stats'].empty:
        df = data_dict['assembly_stats']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # Contig count
        df['contig_count'].plot(kind='bar', ax=ax1, color='skyblue')
        ax1.set_xlabel('Sample')
        ax1.set_ylabel('Number of Contigs')
        ax1.set_title('Contig Count by Sample')
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
        
        # N50 values
        df['n50'].plot(kind='bar', ax=ax2, color='lightgreen')
        ax2.set_xlabel('Sample')
        ax2.set_ylabel('N50 (bp)')
        ax2.set_title('N50 by Sample')
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'assembly_stats.png')
        plt.savefig(plot_file, dpi=300)
        plt.close()
        plot_files['assembly_stats'] = plot_file
    
    # 5. Bin quality plot
    if 'binning_stats' in data_dict and not data_dict['binning_stats'].empty:
        df = data_dict['binning_stats']
        
        if 'Completeness' in df.columns and 'Contamination' in df.columns:
            fig, ax = plt.subplots(figsize=(10, 8))
            
            scatter = ax.scatter(
                df['Contamination'], 
                df['Completeness'],
                c=df['Completeness'] - 5*df['Contamination'],  # Quality score coloring
                cmap='viridis',
                s=80,
                alpha=0.7
            )
            
            # Add bin labels
            for i, bin_id in enumerate(df['bin_id']):
                ax.annotate(bin_id, (df['Contamination'].iloc[i], df['Completeness'].iloc[i]),
                           fontsize=8, ha='right', va='bottom')
            
            # Add quality thresholds
            ax.axhline(y=90, color='r', linestyle='--', alpha=0.3)
            ax.axhline(y=50, color='orange', linestyle='--', alpha=0.3)
            ax.axvline(x=5, color='r', linestyle='--', alpha=0.3)
            
            ax.set_xlabel('Contamination (%)')
            ax.set_ylabel('Completeness (%)')
            ax.set_title('MAG Quality Assessment')
            
            # Add colorbar
            cbar = plt.colorbar(scatter)
            cbar.set_label('Quality Score (Completeness - 5*Contamination)')
            
            # Add quality legend
            ax.text(0.02, 0.98, "High-quality: >90% complete, <5% contamination", 
                   transform=ax.transAxes, fontsize=9, va='top', color='darkgreen')
            ax.text(0.02, 0.94, "Medium-quality: >50% complete, <10% contamination", 
                   transform=ax.transAxes, fontsize=9, va='top', color='darkorange')
            
            plt.tight_layout()
            plot_file = os.path.join(output_dir, 'bin_quality.png')
            plt.savefig(plot_file, dpi=300)
            plt.close()
            plot_files['bin_quality'] = plot_file
    
    # 6. Metabolic pathways heatmap
    if 'pathway_data' in data_dict and not data_dict['pathway_data'].empty:
        df = data_dict['pathway_data']
        
        # Create a pivot table for the heatmap
        pivot_df = df.pivot(index='pathway', columns='sample', values='abundance')
        
        # Normalize by sample (column)
        norm_df = pivot_df.apply(lambda x: x / x.sum() * 100, axis=0)
        
        fig, ax = plt.subplots(figsize=(12, 10))
        heatmap = sns.heatmap(norm_df, annot=False, cmap='YlGnBu', fmt='.1f', 
                            linewidths=0.5, ax=ax)
        
        ax.set_title('Top Metabolic Pathways (Normalized Abundance)')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'pathways.png')
        plt.savefig(plot_file, dpi=300)
        plt.close()
        plot_files['pathways'] = plot_file
    
    return plot_files

def generate_html_report(data_dict, plot_files, output_file):
    """Generate an HTML report with all the collected data"""
    # HTML template
    template_str = '''
    <!DOCTYPE html>
    <html>
    <head>
        <title>Metagenomic Analysis Report</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                line-height: 1.6;
                margin: 20px;
                color: #333;
            }
            h1, h2, h3 {
                color: #2c3e50;
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
            }
            .section {
                margin-bottom: 30px;
                border: 1px solid #ddd;
                border-radius: 5px;
                padding: 20px;
                background-color: #f9f9f9;
            }
            table {
                border-collapse: collapse;
                width: 100%;
                margin-bottom: 20px;
            }
            th, td {
                text-align: left;
                padding: 8px;
                border: 1px solid #ddd;
            }
            th {
                background-color: #4CAF50;
                color: white;
            }
            tr:nth-child(even) {
                background-color: #f2f2f2;
            }
            .plot-container {
                text-align: center;
                margin: 20px 0;
            }
            .plot-container img {
                max-width: 100%;
                height: auto;
                border: 1px solid #ddd;
                border-radius: 5px;
                box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            }
            .header {
                background-color: #2c3e50;
                color: white;
                padding: 20px;
                border-radius: 5px;
                margin-bottom: 20px;
            }
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>Metagenomic Analysis Report</h1>
                <p>Generated on {{ date }}</p>
            </div>
            
            <div class="section">
                <h2>Summary</h2>
                <p>This report summarizes the metagenomic analysis results.</p>
                <ul>
                    <li><strong>Number of samples:</strong> {{ sample_count }}</li>
                    <li><strong>Pipeline version:</strong> 1.0.0</li>
                </ul>
            </div>
            
            {% if fastp_stats is defined and fastp_stats|length > 0 %}
            <div class="section">
                <h2>Quality Control Statistics</h2>
                
                {% if 'qc_stats' in plot_files %}
                <div class="plot-container">
                    <img src="{{ plot_files.qc_stats }}" alt="QC Statistics">
                </div>
                {% endif %}
                
                <table>
                    <tr>
                        <th>Sample</th>
                        <th>Raw Reads</th>
                        <th>Clean Reads</th>
                        <th>Retention (%)</th>
                        <th>Q30 Rate (%)</th>
                        <th>GC Content (%)</th>
                    </tr>
                    {% for idx, row in fastp_stats.iterrows() %}
                    <tr>
                        <td>{{ idx }}</td>
                        <td>{{ '{:,}'.format(row['raw_reads']) }}</td>
                        <td>{{ '{:,}'.format(row['clean_reads']) }}</td>
                        <td>{{ '{:.2f}'.format(row['clean_reads'] / row['raw_reads'] * 100) }}</td>
                        <td>{{ '{:.2f}'.format(row['q30_rate'] * 100) }}</td>
                        <td>{{ '{:.2f}'.format(row['gc_content'] * 100) }}</td>
                    </tr>
                    {% endfor %}
                </table>
            </div>
            {% endif %}
            
            {% if host_stats is defined and host_stats|length > 0 %}
            <div class="section">
                <h2>Host Removal Statistics</h2>
                
                {% if 'host_removal' in plot_files %}
                <div class="plot-container">
                    <img src="{{ plot_files.host_removal }}" alt="Host Removal">
                </div>
                {% endif %}
                
                <table>
                    <tr>
                        <th>Sample</th>
                        <th>Host Reads (%)</th>
                        <th>Non-host Reads (%)</th>
                    </tr>
                    {% for idx, row in host_stats.iterrows() %}
                    <tr>
                        <td>{{ idx }}</td>
                        <td>{{ '{:.2f}'.format(row['host_pct']) }}</td>
                        <td>{{ '{:.2f}'.format(row['non_host_pct']) }}</td>
                    </tr>
                    {% endfor %}
                </table>
            </div>
            {% endif %}
            
            {% if taxonomy_data is defined and taxonomy_data|length > 0 %}
            <div class="section">
                <h2>Taxonomic Composition</h2>
                
                {% if 'taxonomy' in plot_files %}
                <div class="plot-container">
                    <img src="{{ plot_files.taxonomy }}" alt="Taxonomic Composition">
                </div>
                {% endif %}
                
                <h3>Top Species by Sample</h3>
                {% for sample in taxonomy_data['sample'].unique() %}
                <h4>{{ sample }}</h4>
                <table>
                    <tr>
                        <th>Species</th>
                        <th>Relative Abundance (%)</th>
                    </tr>
                    {% for _, row in taxonomy_data[taxonomy_data['sample'] == sample].iterrows() %}
                    <tr>
                        <td>{{ row['taxa'] }}</td>
                        <td>{{ '{:.2f}'.format(row['fraction']) }}</td>
                    </tr>
                    {% endfor %}
                </table>
                {% endfor %}
            </div>
            {% endif %}
            
            {% if assembly_stats is defined and assembly_stats|length > 0 %}
            <div class="section">
                <h2>Assembly Statistics</h2>
                
                {% if 'assembly_stats' in plot_files %}
                <div class="plot-container">
                    <img src="{{ plot_files.assembly_stats }}" alt="Assembly Statistics">
                </div>
                {% endif %}
                
                <table>
                    <tr>
                        <th>Sample</th>
                        <th>Contigs</th>
                        <th>Total Length (Mbp)</th>
                        <th>N50 (bp)</th>
                        <th>Longest Contig (bp)</th>
                        <th>Average Length (bp)</th>
                    </tr>
                    {% for idx, row in assembly_stats.iterrows() %}
                    <tr>
                        <td>{{ idx }}</td>
                        <td>{{ '{:,}'.format(row['contig_count']) }}</td>
                        <td>{{ '{:.2f}'.format(row['total_length'] / 1000000) }}</td>
                        <td>{{ '{:,}'.format(row['n50']) }}</td>
                        <td>{{ '{:,}'.format(row['longest_contig']) }}</td>
                        <td>{{ '{:,}'.format(row['avg_contig_len']|int) }}</td>
                    </tr>
                    {% endfor %}
                </table>
            </div>
            {% endif %}
            
            {% if binning_stats is defined and binning_stats|length > 0 %}
            <div class="section">
                <h2>Genome Bin Quality</h2>
                
                {% if 'bin_quality' in plot_files %}
                <div class="plot-container">
                    <img src="{{ plot_files.bin_quality }}" alt="Bin Quality">
                </div>
                {% endif %}
                
                <table>
                    <tr>
                        <th>Bin ID</th>
                        <th>Completeness (%)</th>
                        <th>Contamination (%)</th>
                        <th>Quality Score</th>
                        <th>Genome Size (Mbp)</th>
                        <th>GC (%)</th>
                        <th>Taxonomy</th>
                    </tr>
                    {% for _, row in binning_stats.iterrows() %}
                    <tr>
                        <td>{{ row['bin_id'] }}</td>
                        <td>{{ '{:.2f}'.format(row['Completeness']) }}</td>
                        <td>{{ '{:.2f}'.format(row['Contamination']) }}</td>
                        <td>{{ '{:.2f}'.format(row['Completeness'] - 5*row['Contamination']) }}</td>
                        <td>{{ '{:.2f}'.format(row['Genome_size'] / 1000000) if 'Genome_size' in row else 'N/A' }}</td>
                        <td>{{ '{:.2f}'.format(row['GC']) if 'GC' in row else 'N/A' }}</td>
                        <td>{{ row['Genus'] + ' ' + row['Species'] if 'Genus' in row and 'Species' in row and row['Species'] else row.get('Taxonomy', 'Unknown') }}</td>
                    </tr>
                    {% endfor %}
                </table>
            </div>
            {% endif %}
            
            {% if pathway_data is defined and pathway_data|length > 0 %}
            <div class="section">
                <h2>Functional Profiling</h2>
                
                {% if 'pathways' in plot_files %}
                <div class="plot-container">
                    <img src="{{ plot_files.pathways }}" alt="Metabolic Pathways">
                </div>
                {% endif %}
                
                <h3>Top Metabolic Pathways by Sample</h3>
                {% for sample in pathway_data['sample'].unique() %}
                <h4>{{ sample }}</h4>
                <table>
                    <tr>
                        <th>Pathway</th>
                        <th>Abundance</th>
                    </tr>
                    {% for _, row in pathway_data[pathway_data['sample'] == sample].iterrows() %}
                    <tr>
                        <td>{{ row['pathway'] }}</td>
                        <td>{{ '{:.2f}'.format(row['abundance']) }}</td>
                    </tr>
                    {% endfor %}
                </table>
                {% endfor %}
            </div>
            {% endif %}
            
            <div class="section">
                <h2>Methods</h2>
                <p>This analysis was performed using the following bioinformatics tools:</p>
                <ul>
                    <li><strong>Quality Control:</strong> Fastp v0.23.2</li>
                    <li><strong>Host Removal:</strong> Bowtie2 v2.4.5</li>
                    <li><strong>Taxonomic Classification:</strong> Kraken2 v2.1.2 + Bracken v2.7</li>
                    <li><strong>Functional Profiling:</strong> HUMAnN3 v3.6</li>
                    <li><strong>Assembly:</strong> MEGAHIT v1.2.9</li>
                    <li><strong>Binning:</strong> MetaBAT2 v2.15 + DAS Tool v1.1.5</li>
                    <li><strong>Bin Quality:</strong> CheckM v1.1.3</li>
                    <li><strong>Taxonomic Placement:</strong> GTDB-Tk v2.1.0</li>
                    <li><strong>Annotation:</strong> Prokka v1.14.6 + DRAM v1.4.0</li>
                </ul>
            </div>
        </div>
    </body>
    </html>
    '''
    
    template = Template(template_str)
    
    # Render HTML
    html = template.render(
        date=datetime.now().strftime('%Y-%m-%d %H:%M'),
        sample_count=len(data_dict.get('fastp_stats', [])),
        fastp_stats=data_dict.get('fastp_stats', pd.DataFrame()),
        host_stats=data_dict.get('host_stats', pd.DataFrame()),
        taxonomy_data=data_dict.get('taxonomy_data', pd.DataFrame()),
        assembly_stats=data_dict.get('assembly_stats', pd.DataFrame()),
        binning_stats=data_dict.get('binning_stats', pd.DataFrame()),
        pathway_data=data_dict.get('pathway_data', pd.DataFrame()),
        plot_files=plot_files
    )
    
    # Write HTML to file
    with open(output_file, 'w') as f:
        f.write(html)
    
    return output_file

def main():
    args = parse_args()
    config = load_config(args.config)
    
    # Collect data
    fastp_stats = get_fastp_stats(config)
    host_stats = get_host_removal_stats(config)
    taxonomy_data = get_kraken_top_taxa(config)
    assembly_stats = get_assembly_stats(config)
    binning_stats = get_binning_stats(config)
    pathway_data = get_humann_pathways(config)
    
    # Organize data
    data_dict = {
        'fastp_stats': fastp_stats,
        'host_stats': host_stats,
        'taxonomy_data': taxonomy_data,
        'assembly_stats': assembly_stats,
        'binning_stats': binning_stats,
        'pathway_data': pathway_data
    }
    
    # Generate plots
    plot_files = generate_plots(data_dict)
    
    # Generate HTML report
    report_file = generate_html_report(data_dict, plot_files, args.output)
    
    print(f"Report generated: {report_file}")

if __name__ == "__main__":
    main()

