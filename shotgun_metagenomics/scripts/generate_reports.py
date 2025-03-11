#!/usr/bin/env python3
"""
Generate HTML reports for taxonomy, AMR, strain analysis, and outbreak tracking.
"""

import pandas as pd
import numpy as np
import json
import os
import glob
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Environment, FileSystemLoader
from snakemake.utils import logger

def parse_kraken_reports(report_files):
    """Parse Kraken2 reports and extract taxonomic information"""
    taxonomy_data = {}
    
    for report_file in report_files:
        sample_id = os.path.basename(report_file).split('.')[0]
        taxonomy_data[sample_id] = {}
        
        try:
            # Read Kraken2 report
            with open(report_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 6:
                        continue
                    
                    percentage = float(parts[0])
                    tax_level = parts[3]
                    tax_id = parts[4]
                    name = parts[5].strip()
                    
                    # Only keep species level and above with percentage > 0.1%
                    if tax_level in ['S', 'G', 'F', 'O', 'C', 'P', 'K'] and percentage > 0.1:
                        taxonomy_data[sample_id][name] = {
                            'percentage': percentage,
                            'tax_level': tax_level,
                            'tax_id': tax_id
                        }
        
        except Exception as e:
            logger.warning(f"Error processing {report_file}: {e}")
    
    return taxonomy_data

def parse_metaphlan_results(metaphlan_file):
    """Parse merged MetaPhlAn results"""
    try:
        df = pd.read_csv(metaphlan_file, sep='\t', index_col=0)
        return df
    except Exception as e:
        logger.warning(f"Error processing {metaphlan_file}: {e}")
        return pd.DataFrame()

def parse_resfinder_results(resfinder_files):
    """Parse ResFinder results from multiple samples"""
    amr_data = {}
    
    for json_file in resfinder_files:
        sample_id = os.path.basename(json_file).split('.')[0]
        amr_data[sample_id] = {}
        
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Extract antimicrobial resistance genes
            if 'seq_regions' in data:
                for gene_id, gene_data in data['seq_regions'].items():
                    gene_name = gene_data.get('name', gene_id)
                    resistance = gene_data.get('resistance', [])
                    identity = gene_data.get('identity', 0)
                    coverage = gene_data.get('coverage', 0)
                    
                    amr_data[sample_id][gene_name] = {
                        'resistance': resistance,
                        'identity': identity,
                        'coverage': coverage
                    }
        
        except Exception as e:
            logger.warning(f"Error processing {json_file}: {e}")
    
    return amr_data

def parse_deeparg_results(deeparg_file):
    """Parse combined DeepARG results"""
    try:
        df = pd.read_csv(deeparg_file, sep='\t')
        return df
    except Exception as e:
        logger.warning(f"Error processing {deeparg_file}: {e}")
        return pd.DataFrame()

def parse_strain_results(strain_file):
    """Parse strain analysis and outbreak tracking results"""
    try:
        df = pd.read_csv(strain_file, sep='\t')
        return df
    except Exception as e:
        logger.warning(f"Error processing {strain_file}: {e}")
        return pd.DataFrame()

def create_taxonomy_plots(kraken_data, metaphlan_df, output_dir):
    """Create taxonomy visualization plots"""
    os.makedirs(output_dir, exist_ok=True)
    plot_files = []
    
    # Create stacked bar plot for top 10 species across samples from Kraken2 data
    if kraken_data:
        # Combine all species across samples
        all_species = set()
        for sample in kraken_data.values():
            all_species.update(sample.keys())
        
        # Filter for species level (usually starts with no spaces, then one space)
        species_only = [sp for sp in all_species if sp.lstrip()[0].isupper() and sp.count(' ') == 1]
        
        # Create a DataFrame with samples as columns and species as rows
        data = []
        for species in species_only:
            row = {'Species': species}
            for sample_id, sample_data in kraken_data.items():
                row[sample_id] = sample_data.get(species, {}).get('percentage', 0)
            data.append(row)
        
        species_df = pd.DataFrame(data)
        species_df = species_df.set_index('Species')
        
        # Sort by total abundance and get top 10
        species_df['total'] = species_df.sum(axis=1)
        top_species = species_df.sort_values('total', ascending=False).head(10)
        top_species = top_species.drop('total', axis=1)
        
        # Create stacked bar plot
        plt.figure(figsize=(12, 8))
        top_species.T.plot(kind='bar', stacked=True, figsize=(12, 8))
        plt.title('Top 10 Species Abundance Across Samples (Kraken2)')
        plt.xlabel('Sample')
        plt.ylabel('Relative Abundance (%)')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Save the plot
        kraken_plot = os.path.join(output_dir, 'kraken2_top_species.png')
        plt.savefig(kraken_plot)
        plt.close()
        plot_files.append(('Kraken2 Top Species', os.path.basename(kraken_plot)))
    
    # Create heatmap from MetaPhlAn data
    if not metaphlan_df.empty:
        # Filter for species level
        species_rows = [idx for idx in metaphlan_df.index if len(idx.split('|')) == 7]
        species_df = metaphlan_df.loc[species_rows]
        
        # Extract species names
        species_df['Species'] = species_df.index.map(lambda x: x.split('|')[-1].replace('s__', ''))
        species_df = species_df.set_index('Species')
        
        # Sort by mean abundance and get top 15
        species_df['mean'] = species_df.mean(axis=1)
        top_species = species_df.sort_values('mean', ascending=False).head(15)
        top_species = top_species.drop('mean', axis=1)
        
        # Create heatmap
        plt.figure(figsize=(12, 10))
        sns.heatmap(top_species, annot=False, cmap='viridis', linewidths=.5)
        plt.title('Top 15 Species Abundance Heatmap (MetaPhlAn)')
        plt.tight_layout()
        
        # Save the plot
        metaphlan_plot = os.path.join(output_dir, 'metaphlan_heatmap.png')
        plt.savefig(metaphlan_plot)
        plt.close()
        plot_files.append(('MetaPhlAn Heatmap', os.path.basename(metaphlan_plot)))
    
    return plot_files

def create_amr_plots(resfinder_data, deeparg_df, output_dir):
    """Create AMR visualization plots"""
    os.makedirs(output_dir, exist_ok=True)
    plot_files = []
    
    # Create bar plot for ResFinder results
    if resfinder_data:
        # Count the number of resistance genes per antimicrobial class
        amr_classes = {}
        
        for sample_id, genes in resfinder_data.items():
            sample_classes = {}
            
            for gene_name, gene_data in genes.items():
                for resistance in gene_data['resistance']:
                    if resistance not in sample_classes:
                        sample_classes[resistance] = 0
                    sample_classes[resistance] += 1
            
            for amr_class, count in sample_classes.items():
                if amr_class not in amr_classes:
                    amr_classes[amr_class] = {}
                amr_classes[amr_class][sample_id] = count
        
        # Convert to DataFrame
        amr_df = pd.DataFrame(amr_classes).fillna(0)
        
        # Plot
        plt.figure(figsize=(14, 8))
        amr_df.plot(kind='bar', figsize=(14, 8))
        plt.title('Antimicrobial Resistance Genes by Class (ResFinder)')
        plt.xlabel('Sample')
        plt.ylabel('Number of Resistance Genes')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Save the plot
        resfinder_plot = os.path.join(output_dir, 'resfinder_amr_classes.png')
        plt.savefig(resfinder_plot)
        plt.close()
        plot_files.append(('ResFinder AMR Classes', os.path.basename(resfinder_plot)))
    
    # Create heatmap for DeepARG results
    if not deeparg_df.empty:
        # Get ARG categories and sample columns
        arg_categories = deeparg_df['ARG-category'].unique()
        sample_columns = [col for col in deeparg_df.columns if col not in ['ARG', 'ARG-category', 'ARG-type']]
        
        # Create a DataFrame with categories as rows and samples as columns
        category_data = {}
        
        for category in arg_categories:
            category_df = deeparg_df[deeparg_df['ARG-category'] == category]
            category_sum = category_df[sample_columns].sum()
            category_data[category] = category_sum
        
        category_df = pd.DataFrame(category_data).T
        
        # Create heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(category_df, annot=True, fmt=".2f", cmap='YlOrRd', linewidths=.5)
        plt.title('ARG Categories Abundance Heatmap (DeepARG)')
        plt.tight_layout()
        
        # Save the plot
        deeparg_plot = os.path.join(output_dir, 'deeparg_heatmap.png')
        plt.savefig(deeparg_plot)
        plt.close()
        plot_files.append(('DeepARG Heatmap', os.path.basename(deeparg_plot)))
    
    return plot_files

def create_strain_plots(strain_df, output_dir):
    """Create strain analysis and outbreak tracking plots"""
    os.makedirs(output_dir, exist_ok=True)
    plot_files = []
    
    if not strain_df.empty:
        # Create a network visualization of outbreak clusters
        # Count samples per cluster
        cluster_counts = strain_df.groupby(['species', 'cluster_id']).size().reset_index(name='sample_count')
        
        # Plot
        plt.figure(figsize=(10, 6))
        sns.barplot(data=cluster_counts, x='species', y='sample_count', hue='cluster_id')
        plt.title('Potential Outbreak Clusters')
        plt.xlabel('Species')
        plt.ylabel('Number of Samples')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        # Save the plot
        cluster_plot = os.path.join(output_dir, 'outbreak_clusters.png')
        plt.savefig(cluster_plot)
        plt.close()
        plot_files.append(('Outbreak Clusters', os.path.basename(cluster_plot)))
        
        # Create similarity heatmap
        pivot_df = strain_df.pivot_table(
            index='sample_id', 
            columns='species', 
            values='strain_similarity',
            aggfunc='first'
        ).fillna(0)
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(pivot_df, annot=True, fmt=".2f", cmap='viridis', linewidths=.5)
        plt.title('Strain Similarity Heatmap')
        plt.tight_layout()
        
        # Save the plot
        similarity_plot = os.path.join(output_dir, 'strain_similarity.png')
        plt.savefig(similarity_plot)
        plt.close()
        plot_files.append(('Strain Similarity', os.path.basename(similarity_plot)))
    
    return plot_files

def generate_html_reports(taxonomy_plots, amr_plots, strain_plots, output_files):
    """Generate HTML reports using the plots"""
    # Create a simple HTML template
    template = """<!DOCTYPE html>
<html>
<head>
    <title>{{ title }}</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .header { background-color: #f2f2f2; padding: 10px; border-radius: 5px; }
        .plot { margin: 20px 0; }
        .plot img { max-width: 100%; border: 1px solid #ddd; }
        .footer { margin-top: 30px; font-size: 0.8em; color: #666; }
    </style>
</head>
<body>
    <div class="header">
        <h1>{{ title }}</h1>
        <p>Generated on {{ date }}</p>
    </div>
    
    <div class="content">
        <p>{{ description }}</p>
        
        {% for plot_title, plot_file in plots %}
        <div class="plot">
            <h2>{{ plot_title }}</h2>
            <img src="{{ plot_dir }}/{{ plot_file }}" alt="{{ plot_title }}">
        </div>
        {% endfor %}
    </div>
    
    <div class="footer">
        <p>Report generated by the Shotgun Metagenomics Analysis Pipeline</p>
    </div>
</body>
</html>
"""
    
    from datetime import datetime
    today = datetime.now().strftime('%Y-%m-%d')
    
    # Generate taxonomy report
    with open(output_files['taxonomy_report'], 'w') as f:
        html = template.replace('{{ title }}', 'Taxonomic Profiling Report')
        html = html.replace('{{ date }}', today)
        html = html.replace('{{ description }}', 'This report shows the taxonomic composition of the analyzed samples.')
        
        plots_html = ''
        for plot_title, plot_file in taxonomy_plots:
            plots_html += f'<div class="plot">\n'
            plots_html += f'    <h2>{plot_title}</h2>\n'
            plots_html += f'    <img src="../plots/{plot_file}" alt="{plot_title}">\n'
            plots_html += f'</div>\n'
        
        html = html.replace('{% for plot_title, plot_file in plots %}\n        <div class="plot">\n            <h2>{{ plot_title }}</h2>\n            <img src="{{ plot_dir }}/{{ plot_file }}" alt="{{ plot_title }}">\n        </div>\n        {% endfor %}', plots_html)
        
        f.write(html)
    
    # Generate AMR report
    with open(output_files['amr_report'], 'w') as f:
        html = template.replace('{{ title }}', 'Antimicrobial Resistance Report')
        html = html.replace('{{ date }}', today)
        html = html.replace('{{ description }}', 'This report shows the antimicrobial resistance genes and their abundance in the analyzed samples.')
        
        plots_html = ''
        for plot_title, plot_file in amr_plots:
            plots_html += f'<div class="plot">\n'
            plots_html += f'    <h2>{plot_title}</h2>\n'
            plots_html += f'    <img src="../plots/{plot_file}" alt="{plot_title}">\n'
            plots_html += f'</div>\n'
        
        html = html.replace('{% for plot_title, plot_file in plots %}\n        <div class="plot">\n            <h2>{{ plot_title }}</h2>\n            <img src="{{ plot_dir }}/{{ plot_file }}" alt="{{ plot_title }}">\n        </div>\n        {% endfor %}', plots_html)
        
        f.write(html)
    
    # Generate strain report
    with open(output_files['strain_report'], 'w') as f:
        html = template.replace('{{ title }}', 'Strain Analysis Report')
        html = html.replace('{{ date }}', today)
        html = html.replace('{{ description }}', 'This report shows the strain-level analysis of the samples, including strain identification and sequence typing.')
        
        plots_html = ''
        for plot_title, plot_file in strain_plots:
            plots_html += f'<div class="plot">\n'
            plots_html += f'    <h2>{plot_title}</h2>\n'
            plots_html += f'    <img src="../plots/{plot_file}" alt="{plot_title}">\n'
            plots_html += f'</div>\n'
        
        html = html.replace('{% for plot_title, plot_file in plots %}\n        <div class="plot">\n            <h2>{{ plot_title }}</h2>\n            <img src="{{ plot_dir }}/{{ plot_file }}" alt="{{ plot_title }}">\n        </div>\n        {% endfor %}', plots_html)
        
        f.write(html)
    
    # Generate outbreak report (clone of strain report with focus on outbreaks)
    with open(output_files['outbreak_report'], 'w') as f:
        html = template.replace('{{ title }}', 'Outbreak Tracking Report')
        html = html.replace('{{ date }}', today)
        html = html.replace('{{ description }}', 'This report shows potential outbreak clusters identified in the samples.')
        
        plots_html = ''
        for plot_title, plot_file in strain_plots:
            if 'Outbreak' in plot_title:
                plots_html += f'<div class="plot">\n'
                plots_html += f'    <h2>{plot_title}</h2>\n'
                plots_html += f'    <img src="../plots/{plot_file}" alt="{plot_title}">\n'
                plots_html += f'</div>\n'
        
        html = html.replace('{% for plot_title, plot_file in plots %}\n        <div class="plot">\n            <h2>{{ plot_title }}</h2>\n            <img src="{{ plot_dir }}/{{ plot_file }}" alt="{{ plot_title }}">\n        </div>\n        {% endfor %}', plots_html)
        
        f.write(html)

def main():
    # Get input files from Snakemake
    kraken_files = snakemake.input.kraken
    metaphlan_file = snakemake.input.metaphlan
    resfinder_files = snakemake.input.resfinder
    deeparg_file = snakemake.input.deeparg
    strain_file = snakemake.input.strain
    
    # Get output files
    output_files = {
        'taxonomy_report': snakemake.output.taxonomy_report,
        'amr_report': snakemake.output.amr_report,
        'strain_report': snakemake.output.strain_report,
        'outbreak_report': snakemake.output.outbreak_report
    }
    
    # Create plots directory
    plots_dir = os.path.join(os.path.dirname(output_files['taxonomy_report']), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Parse input files
    kraken_data = parse_kraken_reports(kraken_files)
    metaphlan_df = parse_metaphlan_results(metaphlan_file)
    resfinder_data = parse_resfinder_results(resfinder_files)
    deeparg_df = parse_deeparg_results(deeparg_file)
    strain_df = parse_strain_results(strain_file)
    
    # Create visualization plots
    taxonomy_plots = create_taxonomy_plots(kraken_data, metaphlan_df, plots_dir)
    amr_plots = create_amr_plots(resfinder_data, deeparg_df, plots_dir)
    strain_plots = create_strain_plots(strain_df, plots_dir)
    
    # Generate HTML reports
    generate_html_reports(taxonomy_plots, amr_plots, strain_plots, output_files)
    
    logger.info("Generated all reports successfully")

if __name__ == "__main__":
    main()

