#!/usr/bin/env python
"""
Generate comprehensive report integrating results from anomaly detection,
taxonomic classification, and viral analysis to identify novel pathogens.
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import logging
import glob
import json
from datetime import datetime
import warnings
from Bio import SeqIO
warnings.filterwarnings("ignore")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('comprehensive_report')

def parse_args():
    parser = argparse.ArgumentParser(description='Generate comprehensive report')
    
    # Input files and directories
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--anomalies_protein', help='Protein anomaly detection results')
    parser.add_argument('--anomalies_dna', help='DNA anomaly detection results')
    parser.add_argument('--visualization_dir', help='Directory with visualization plots')
    
    parser.add_argument('--cat_taxonomy', help='CAT taxonomy file')
    parser.add_argument('--gtdbtk_dir', help='GTDB-Tk output directory')
    parser.add_argument('--ani_dir', help='ANI results directory')
    
    parser.add_argument('--checkv_dir', help='CheckV output directory')
    parser.add_argument('--vcontact2_dir', help='vConTACT2 output directory')
    parser.add_argument('--high_confidence_viral', help='High confidence viral contigs file')
    
    # Output files
    parser.add_argument('--output_report', required=True, help='Output HTML report file')
    parser.add_argument('--output_novel_pathogens', required=True, help='Output TSV file for novel pathogens')
    
    return parser.parse_args()

def read_anomaly_results(protein_file, dna_file):
    """Read and process anomaly detection results."""
    protein_anomalies = None
    dna_anomalies = None
    
    if protein_file and os.path.exists(protein_file):
        protein_anomalies = pd.read_csv(protein_file, sep='\t')
        logger.info(f"Loaded protein anomalies: {len(protein_anomalies)} entries")
    
    if dna_file and os.path.exists(dna_file):
        dna_anomalies = pd.read_csv(dna_file, sep='\t')
        logger.info(f"Loaded DNA anomalies: {len(dna_anomalies)} entries")
    
    return protein_anomalies, dna_anomalies

def read_cat_taxonomy(cat_file):
    """Read and process CAT taxonomy results."""
    if not cat_file or not os.path.exists(cat_file):
        return None
    
    try:
        cat_taxonomy = pd.read_csv(cat_file, sep='\t')
        
        # Rename columns if necessary
        if '# contig' in cat_taxonomy.columns:
            cat_taxonomy = cat_taxonomy.rename(columns={'# contig': 'contig_id'})
        
        logger.info(f"Loaded CAT taxonomy: {len(cat_taxonomy)} entries")
        return cat_taxonomy
    except Exception as e:
        logger.warning(f"Error reading CAT taxonomy: {e}")
        return None

def read_gtdbtk_results(gtdbtk_dir):
    """Read and process GTDB-Tk classification results."""
    if not gtdbtk_dir or not os.path.exists(gtdbtk_dir):
        return None, None
    
    bac_file = os.path.join(gtdbtk_dir, 'classify', 'gtdbtk.bac120.summary.tsv')
    ar_file = os.path.join(gtdbtk_dir, 'classify', 'gtdbtk.ar53.summary.tsv')
    
    bac_results = None
    ar_results = None
    
    if os.path.exists(bac_file):
        try:
            bac_results = pd.read_csv(bac_file, sep='\t')
            logger.info(f"Loaded GTDB-Tk bacterial classifications: {len(bac_results)} entries")
        except Exception as e:
            logger.warning(f"Error reading GTDB-Tk bacterial results: {e}")
    
    if os.path.exists(ar_file):
        try:
            ar_results = pd.read_csv(ar_file, sep='\t')
            logger.info(f"Loaded GTDB-Tk archaeal classifications: {len(ar_results)} entries")
        except Exception as e:
            logger.warning(f"Error reading GTDB-Tk archaeal results: {e}")
    
    return bac_results, ar_results

def read_ani_results(ani_dir):
    """Read and process ANI results."""
    if not ani_dir or not os.path.exists(ani_dir):
        return None
    
    summary_file = os.path.join(ani_dir, 'ani_summary.tsv')
    
    if os.path.exists(summary_file):
        try:
            ani_results = pd.read_csv(summary_file, sep='\t')
            logger.info(f"Loaded ANI results: {len(ani_results)} entries")
            return ani_results
        except Exception as e:
            logger.warning(f"Error reading ANI results: {e}")
    
    return None

def read_checkv_results(checkv_dir, high_conf_file):
    """Read and process CheckV results."""
    checkv_results = None
    high_conf_viral_contigs = set()
    
    if checkv_dir and os.path.exists(checkv_dir):
        quality_file = os.path.join(checkv_dir, 'quality_summary.tsv')
        if os.path.exists(quality_file):
            try:
                checkv_results = pd.read_csv(quality_file, sep='\t')
                logger.info(f"Loaded CheckV results: {len(checkv_results)} entries")
            except Exception as e:
                logger.warning(f"Error reading CheckV results: {e}")
    
    if high_conf_file and os.path.exists(high_conf_file):
        try:
            for record in SeqIO.parse(high_conf_file, "fasta"):
                high_conf_viral_contigs.add(record.id)
            logger.info(f"Loaded {len(high_conf_viral_contigs)} high-confidence viral contigs")
        except Exception as e:
            logger.warning(f"Error reading high-confidence viral contigs: {e}")
    
    return checkv_results, high_conf_viral_contigs

def read_vcontact2_results(vcontact2_dir):
    """Read and process vConTACT2 results."""
    if not vcontact2_dir or not os.path.exists(vcontact2_dir):
        return None
    
    genome_clusters_file = os.path.join(vcontact2_dir, 'genome_by_genome_overview.csv')
    
    if os.path.exists(genome_clusters_file):
        try:
            vcontact2_results = pd.read_csv(genome_clusters_file)
            logger.info(f"Loaded vConTACT2 results: {len(vcontact2_results)} entries")
            return vcontact2_results
        except Exception as e:
            logger.warning(f"Error reading vConTACT2 results: {e}")
    
    return None

def identify_potential_novel_pathogens(protein_anomalies, dna_anomalies, cat_taxonomy, 
                                      gtdbtk_bac, gtdbtk_ar, ani_results, 
                                      checkv_results, high_conf_viral, vcontact2_results):
    """Identify potential novel pathogens based on anomaly detection and taxonomic classification."""
    # Initialize results structure
    novel_pathogens = []
    
    # Process anomalies from protein data
    if protein_anomalies is not None:
        # Group by sequence ID to combine results from different methods
        protein_grouped = protein_anomalies.groupby('id')
        
        for protein_id, group in protein_grouped:
            # Extract contig ID from protein ID (assuming format: contig_123_gene_1)
            contig_id = protein_id.split('_')[0]
            
            # Count methods that flagged this as anomaly
            anomaly_methods = group[group['anomaly'] == True]['method'].tolist()
            
            # Skip if not flagged by enough methods (at least 2)
            if len(anomaly_methods) < 2:
                continue
            
            # Get taxonomy from CAT if available
            taxonomy = "Unknown"
            if cat_taxonomy is not None:
                cat_match = cat_taxonomy[cat_taxonomy['contig_id'] == contig_id]
                if not cat_match.empty:
                    taxonomy = cat_match['classification'].iloc[0]
            
            # Check if this is a viral protein
            is_viral = False
            viral_quality = "N/A"
            if checkv_results is not None and contig_id in checkv_results['contig_id'].values:
                is_viral = True
                viral_idx = checkv_results[checkv_results['contig_id'] == contig_id].index[0]
                viral_quality = checkv_results.loc[viral_idx, 'checkv_quality']
            
            # Get highest anomaly score
            max_score = group['score'].max()
            
            # Add to novel pathogens list
            novel_pathogens.append({
                'id': protein_id,
                'contig_id': contig_id,
                'type': 'Protein',
                'anomaly_methods': ','.join(anomaly_methods),
                'anomaly_score': max_score,
                'taxonomy': taxonomy,
                'is_viral': is_viral,
                'viral_quality': viral_quality,
                'ani_match': 'N/A',
                'ani_value': 'N/A',
                'novel_assessment': 'Potential novel protein'
            })
    
    # Process anomalies from DNA data (contigs)
    if dna_anomalies is not None:
        # Group by sequence ID to combine results from different methods
        dna_grouped = dna_anomalies.groupby('id')
        
        for contig_id, group in dna_grouped:
            # Count methods that flagged this as anomaly
            anomaly_methods = group[group['anomaly'] == True]['method'].tolist()
            
            # Skip if not flagged by enough methods (at least 2)
            if len(anomaly_methods) < 2:
                continue
            
            # Get taxonomy from CAT if available
            taxonomy = "Unknown"
            if cat_taxonomy is not None:
                cat_match = cat_taxonomy[cat_taxonomy['contig_id'] == contig_id]
                if not cat_match.empty:
                    taxonomy = cat_match['classification'].iloc[0]
            
            # Check if this is a viral contig
            is_viral = False
            viral_quality = "N/A"
            if checkv_results is not None and contig_id in checkv_results['contig_id'].values:
                is_viral = True
                viral_idx = checkv_results[checkv_results['contig_id'] == contig_id].index[0]
                viral_quality = checkv_results.loc[viral_idx, 'checkv_quality']
            
            # Check for GTDB-Tk and ANI matches
            gtdb_classification = "N/A"
            ani_match = "N/A"
            ani_value = "N/A"
            
            # Look for matches in GTDB-Tk results
            if gtdbtk_bac is not None and contig_id in gtdbtk_bac['user_genome'].values:
                gtdb_idx = gtdbtk_bac[gtdbtk_bac['user_genome'] == contig_id].index[0]
                gtdb_classification = gtdbtk_bac.loc[gtdb_idx, 'classification']
            
            elif gtdbtk_ar is not None and contig_id in gtdbtk_ar['user_genome'].values:
                gtdb_idx = gtdbtk_ar[gtdbtk_ar['user_genome'] == contig_id].index[0]
                gtdb_classification = gtdbtk_ar.loc[gtdb_idx, 'classification']
            
            # Look for matches in ANI results
            if ani_results is not None and contig_id in ani_results['Bin_ID'].values:
                ani_idx = ani_results[ani_results['Bin_ID'] == contig_id].index[0]
                ani_match = ani_results.loc[ani_idx, 'Reference_Genome']
                ani_value = ani_results.loc[ani_idx, 'ANI']
            
            # Get highest anomaly score
            max_score = group['score'].max()
            
            # Determine if it's a potential novel pathogen
            novel_assessment = 'Potential novel contig'
            
            if is_viral:
                if contig_id in high_conf_viral:
                    novel_assessment = 'Novel viral pathogen (high confidence)'
                else:
                    novel_assessment = 'Potential novel viral pathogen'
                
                # Add vConTACT2 information if available
                if vcontact2_results is not None:
                    # Format might vary, try different possible ID formats
                    for col in ['Genome', 'genome']:
                        if col in vcontact2_results.columns and any(vcontact2_results[col].str.contains(contig_id, na=False)):
                            vc2_idx = vcontact2_results[vcontact2_results[col].str.contains(contig_id, na=False)].index[0]
                            
                            for cluster_col in ['VC', 'VC Status', 'VC Subcluster']:
                                if cluster_col in vcontact2_results.columns:
                                    cluster_value = vcontact2_results.loc[vc2_idx, cluster_col]
                                    if pd.notna(cluster_value) and cluster_value != '':
                                        novel_assessment += f' (vConTACT2 cluster: {cluster_value})'
                                        break
            
            elif taxonomy != 'Unknown' and 'bacteria' in taxonomy.lower():
                # For bacterial contigs
                if ani_value != 'N/A' and ani_value != 'NA':
                    try:
                        ani_float = float(ani_value)
                        if ani_float < 95:
                            novel_assessment = 'Novel bacterial species'
                        elif ani_float < 80:
                            novel_assessment = 'Novel bacterial genus'
                        else:
                            novel_assessment = 'Known bacterial species'
                    except (ValueError, TypeError):
                        pass
            
            # Add to novel pathogens list
            novel_pathogens.append({
                'id': contig_id,
                'contig_id': contig_id,
                'type': 'Contig',
                'anomaly_methods': ','.join(anomaly_methods),
                'anomaly_score': max_score,
                'taxonomy': taxonomy,
                'gtdb_classification': gtdb_classification,
                'is_viral': is_viral,
                'viral_quality': viral_quality,
                'ani_match': ani_match,
                'ani_value': ani_value,
                'novel_assessment': novel_assessment
            })
    
    # Convert to DataFrame
    novel_df = pd.DataFrame(novel_pathogens)
    
    if not novel_df.empty:
        # Sort by anomaly score (descending)
        novel_df = novel_df.sort_values('anomaly_score', ascending=False)
    
    return novel_df

def create_summary_stats(protein_anomalies, dna_anomalies, cat_taxonomy, 
                        gtdbtk_bac, gtdbtk_ar, checkv_results, novel_pathogens):
    """Create summary statistics for the report."""
    summary = {}
    
    # Count anomalies
    if protein_anomalies is not None:
        # Count unique proteins flagged as anomalies by at least one method
        unique_protein_anomalies = protein_anomalies[protein_anomalies['anomaly'] == True]['id'].nunique()
        summary['protein_anomalies'] = unique_protein_anomalies
    else:
        summary['protein_anomalies'] = 0
    
    if dna_anomalies is not None:
        # Count unique contigs flagged as anomalies by at least one method
        unique_dna_anomalies = dna_anomalies[dna_anomalies['anomaly'] == True]['id'].nunique()
        summary['dna_anomalies'] = unique_dna_anomalies
    else:
        summary['dna_anomalies'] = 0
    
    # Count taxonomic classifications
    if cat_taxonomy is not None:
        # Count by category
        if 'classification' in cat_taxonomy.columns:
            tax_counts = cat_taxonomy['classification'].value_counts().to_dict()
            # Make keys more friendly
            summary['taxonomic_counts'] = {k.replace(';', '; '): v for k, v in tax_counts.items()}
        else:
            summary['taxonomic_counts'] = {}
    else:
        summary['taxonomic_counts'] = {}
    
    # Count GTDB-Tk classifications
    summary['bacterial_bins'] = 0
    summary['archaeal_bins'] = 0
    
    if gtdbtk_bac is not None:
        summary['bacterial_bins'] = len(gtdbtk_bac)
    
    if gtdbtk_ar is not None:
        summary['archaeal_bins'] = len(gtdbtk_ar)
    
    # Count viral contigs
    summary['viral_contigs'] = 0
    summary['high_quality_viral'] = 0
    summary['complete_viral'] = 0
    
    if checkv_results is not None:
        summary['viral_contigs'] = len(checkv_results)
        
        if 'checkv_quality' in checkv_results.columns:
            # Count high-quality and complete viral genomes
            summary['high_quality_viral'] = sum(checkv_results['checkv_quality'] == 'High-quality')
            summary['complete_viral'] = sum(checkv_results['checkv_quality'] == 'Complete')
    
    # Count novel pathogens
    if novel_pathogens is not None and not novel_pathogens.empty:
        summary['novel_pathogens'] = len(novel_pathogens)
        summary['novel_viral'] = sum(novel_pathogens['is_viral'] == True)
        summary['novel_bacterial'] = sum((novel_pathogens['novel_assessment'].str.contains('bacterial')) & 
                                       (novel_pathogens['novel_assessment'].str.contains('Novel')))
    else:
        summary['novel_pathogens'] = 0
        summary['novel_viral'] = 0
        summary['novel_bacterial'] = 0
    
    return summary

def generate_html_report(sample_id, summary_stats, novel_pathogens, visualization_dir=None):
    """Generate HTML report with all results."""
    # Create HTML header
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Novel Pathogen Detection Report - {sample_id}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                margin: 0;
                padding: 20px;
                max-width: 1200px;
                margin: 0 auto;
            }}
            h1, h2, h3, h4 {{
                color: #2c3e50;
            }}
            h1 {{
                border-bottom: 2px solid #3498db;
                padding-bottom: 10px;
            }}
            .report-header {{
                background-color: #f8f9fa;
                padding: 20px;
                border-radius: 5px;
                margin-bottom: 20px;
            }}
            .summary-box {{
                display: flex;
                flex-wrap: wrap;
                gap: 20px;
                margin-bottom: 30px;
            }}
            .stat-card {{
                background: #fff;
                border: 1px solid #ddd;
                border-radius: 5px;
                padding: 15px;
                flex: 1 1 200px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            }}
            .stat-card h3 {{
                margin-top: 0;
                color: #3498db;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 30px;
            }}
            th, td {{
                padding: 12px 15px;
                border: 1px solid #ddd;
                text-align: left;
            }}
            th {{
                background-color: #f4f4f4;
                font-weight: bold;
            }}
            tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            tr:hover {{
                background-color: #f1f1f1;
            }}
            .high-score {{
                color: #e74c3c;
                font-weight: bold;
            }}
            .medium-score {{
                color: #f39c12;
            }}
            .low-score {{
                color: #2ecc71;
            }}
            .visualizations {{
                display: flex;
                flex-wrap: wrap;
                gap: 20px;
                margin-bottom: 30px;
            }}
            .visualization-card {{
                flex: 1 1 45%;
                max-width: 600px;
                background: #fff;
                border: 1px solid #ddd;
                border-radius: 5px;
                padding: 15px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            }}
            .visualization-card img {{
                max-width: 100%;
                height: auto;
                display: block;
                margin: 0 auto;
            }}
            .tax-chart {{
                width: 100%;
                height: 400px;
                margin-bottom: 30px;
            }}
            footer {{
                text-align: center;
                margin-top: 40px;
                padding-top: 20px;
                border-top: 1px solid #eee;
                color: #777;
            }}
        </style>
    </head>
    <body>
        <div class="report-header">
            <h1>Novel Pathogen Detection Report</h1>
            <p><strong>Sample ID:</strong> {sample_id}</p>
            <p><strong>Date Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <h2>Summary Statistics</h2>
        <div class="summary-box">
            <div class="stat-card">
                <h3>Anomaly Detection</h3>
                <p><strong>Protein Anomalies:</strong> {summary_stats.get('protein_anomalies', 0)}</p>
                <p><strong>Contig Anomalies:</strong> {summary_stats.get('dna_anomalies', 0)}</p>
            </div>
            
            <div class="stat-card">
                <h3>Taxonomic Classification</h3>
                <p><strong>Bacterial Bins:</strong> {summary_stats.get('bacterial_bins', 0)}</p>
                <p><strong>Archaeal Bins:</strong> {summary_stats.get('archaeal_bins', 0)}</p>
            </div>
            
            <div class="stat-card">
                <h3>Viral Analysis</h3>
                <p><strong>Viral Contigs:</strong> {summary_stats.get('viral_contigs', 0)}</p>
                <p><strong>High-Quality Viral:</strong> {summary_stats.get('high_quality_viral', 0)}</p>
                <p><strong>Complete Viral:</strong> {summary_stats.get('complete_viral', 0)}</p>
            </div>
            
            <div class="stat-card">
                <h3>Novel Pathogens</h3>
                <p><strong>Total Potential Novel:</strong> {summary_stats.get('novel_pathogens', 0)}</p>
                <p><strong>Novel Viral Pathogens:</strong> {summary_stats.get('novel_viral', 0)}</p>
                <p><strong>Novel Bacterial Pathogens:</strong> {summary_stats.get('novel_bacterial', 0)}</p>
            </div>
        </div>
    """
    
    # Add visualization section if visualizations are available
    if visualization_dir and os.path.exists(visualization_dir):
        html += f"""
        <h2>Visualizations</h2>
        <div class="visualizations">
        """
        
        # Find visualization images
        plot_files = glob.glob(os.path.join(visualization_dir, '*.png'))
        
        # Sort them by method and type
        for plot_file in sorted(plot_files):
            plot_name = os.path.basename(plot_file)
            # Skip taxonomy plots for this section
            if 'taxonomy' in plot_name:
                continue
            
            # Extract method and dimension reduction from filename
            try:
                method, dim_red = plot_name.replace('.png', '').split('_')
                title = f"{method.replace('_', ' ').title()} - {dim_red.upper()}"
            except:
                title = plot_name
            
            # Embed image
            html += f"""
            <div class="visualization-card">
                <h3>{title}</h3>
                <img src="data:image/png;base64,{encode_image(plot_file)}" alt="{title}" />
            </div>
            """
        
        html += """
        </div>
        """
        
        # Add taxonomy visualizations if available
        taxonomy_plots = [f for f in plot_files if 'taxonomy' in f]
        if taxonomy_plots:
            html += f"""
            <h2>Taxonomy Visualizations</h2>
            <div class="visualizations">
            """
            
            for plot_file in sorted(taxonomy_plots):
                plot_name = os.path.basename(plot_file)
                
                # Extract method and dimension reduction from filename
                try:
                    parts = plot_name.replace('.png', '').split('_')
                    method = parts[0]
                    title = f"Taxonomy with {method.replace('_', ' ').title()}"
                except:
                    title = plot_name
                
                # Embed image
                html += f"""
                <div class="visualization-card">
                    <h3>{title}</h3>
                    <img src="data:image/png;base64,{encode_image(plot_file)}" alt="{title}" />
                </div>
                """
            
            html += """
            </div>
            """
    
    # Add taxonomy distribution chart (placeholder - this would be replaced with actual chart)
    if summary_stats.get('taxonomic_counts'):
        tax_data = []
        for tax, count in summary_stats['taxonomic_counts'].items():
            # Simplify taxonomy names for display
            if len(tax) > 30:
                parts = tax.split(';')
                if len(parts) > 2:
                    # Take first and last two parts
                    simplified = f"{parts[0]};...;{parts[-2]};{parts[-1]}"
                else:
                    simplified = tax
            else:
                simplified = tax
            
            tax_data.append({'name': simplified, 'count': count})
        
        # Sort by count descending
        tax_data = sorted(tax_data, key=lambda x: x['count'], reverse=True)
        
        # Take top 10
        tax_data = tax_data[:10]
        
        # Create chart data
        chart_labels = [item['name'] for item in tax_data]
        chart_values = [item['count'] for item in tax_data]
        
        # Add chart
        html += f"""
        <h2>Taxonomic Distribution</h2>
        <div class="tax-chart">
            <canvas id="taxonomyChart"></canvas>
        </div>
        
        <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
        <script>
            document.addEventListener('DOMContentLoaded', function() {{
                const ctx = document.getElementById('taxonomyChart').getContext('2d');
                const taxonomyChart = new Chart(ctx, {{
                    type: 'bar',
                    data: {{
                        labels: {json.dumps(chart_labels)},
                        datasets: [{{
                            label: 'Count',
                            data: {json.dumps(chart_values)},
                            backgroundColor: 'rgba(54, 162, 235, 0.5)',
                            borderColor: 'rgba(54, 162, 235, 1)',
                            borderWidth: 1
                        }}]
                    }},
                    options: {{
                        indexAxis: 'y',
                        plugins: {{
                            legend: {{
                                display: false
                            }},
                            title: {{
                                display: true,
                                text: 'Top 10 Taxonomic Classifications'
                            }}
                        }},
                        scales: {{
                            x: {{
                                beginAtZero: true,
                                title: {{
                                    display: true,
                                    text: 'Count'
                                }}
                            }},
                            y: {{
                                title: {{
                                    display: true,
                                    text: 'Taxonomy'
                                }}
                            }}
                        }}
                    }}
                }});
            }});
        </script>
        """
    
    # Add novel pathogens table
    if novel_pathogens is not None and not novel_pathogens.empty:
        html += """
        <h2>Potential Novel Pathogens</h2>
        <table>
            <thead>
                <tr>
                    <th>ID</th>
                    <th>Type</th>
                    <th>Anomaly Score</th>
                    <th>Taxonomy</th>
                    <th>Viral</th>
                    <th>ANI</th>
                    <th>Assessment</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for _, row in novel_pathogens.iterrows():
            # Format anomaly score with color based on value
            try:
                score = float(row['anomaly_score'])
                if score > 0.8:
                    score_class = "high-score"
                elif score > 0.5:
                    score_class = "medium-score"
                else:
                    score_class = "low-score"
                score_formatted = f'<span class="{score_class}">{score:.3f}</span>'
            except (ValueError, TypeError):
                score_formatted = row['anomaly_score']
            
            # Format ANI value
            ani_value = row['ani_value']
            if ani_value != 'N/A' and ani_value != 'NA':
                try:
                    ani_formatted = f"{float(ani_value):.2f}%"
                except (ValueError, TypeError):
                    ani_formatted = ani_value
            else:
                ani_formatted = "N/A"
            
            html += f"""
                <tr>
                    <td>{row['id']}</td>
                    <td>{row['type']}</td>
                    <td>{score_formatted}</td>
                    <td>{row['taxonomy']}</td>
                    <td>{'Yes' if row['is_viral'] else 'No'}</td>
                    <td>{ani_formatted}</td>
                    <td>{row['novel_assessment']}</td>
                </tr>
            """
        
        html += """
            </tbody>
        </table>
        """
    
    # Footer
    html += """
        <footer>
            <p>Generated by Novel Pathogen Detection Pipeline</p>
        </footer>
        
    </body>
    </html>
    """
    
    return html

def encode_image(image_path):
    """Encode image as base64 string for embedding in HTML."""
    import base64
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def main():
    args = parse_args()
    
    # Read results files
    protein_anomalies, dna_anomalies = read_anomaly_results(
        args.anomalies_protein, args.anomalies_dna
    )
    
    cat_taxonomy = read_cat_taxonomy(args.cat_taxonomy)
    
    gtdbtk_bac, gtdbtk_ar = read_gtdbtk_results(args.gtdbtk_dir)
    
    ani_results = read_ani_results(args.ani_dir)
    
    checkv_results, high_conf_viral = read_checkv_results(
        args.checkv_dir, args.high_confidence_viral
    )
    
    vcontact2_results = read_vcontact2_results(args.vcontact2_dir)
    
    # Identify potential novel pathogens
    novel_pathogens = identify_potential_novel_pathogens(
        protein_anomalies, dna_anomalies, cat_taxonomy,
        gtdbtk_bac, gtdbtk_ar, ani_results,
        checkv_results, high_conf_viral, vcontact2_results
    )
    
    # Create summary statistics
    summary_stats = create_summary_stats(
        protein_anomalies, dna_anomalies, cat_taxonomy,
        gtdbtk_bac, gtdbtk_ar, checkv_results, novel_pathogens
    )
    
    # Generate HTML report
    html_report = generate_html_report(
        args.sample_id, summary_stats, novel_pathogens, args.visualization_dir
    )
    
    # Save HTML report
    os.makedirs(os.path.dirname(args.output_report), exist_ok=True)
    with open(args.output_report, 'w') as f:
        f.write(html_report)
    logger.info(f"Saved HTML report to {args.output_report}")
    
    # Save novel pathogens list
    if novel_pathogens is not None and not novel_pathogens.empty:
        novel_pathogens.to_csv(args.output_novel_pathogens, sep='\t', index=False)
        logger.info(f"Saved novel pathogens list to {args.output_novel_pathogens}")
    else:
        # Create empty file with headers
        with open(args.output_novel_pathogens, 'w') as f:
            f.write("id\tcontig_id\ttype\tanomaly_methods\tanomaly_score\ttaxonomy\tgtdb_classification\tis_viral\tviral_quality\tani_match\tani_value\tnovel_assessment\n")
        logger.info(f"No novel pathogens identified. Created empty file at {args.output_novel_pathogens}")
    
    logger.info("Report generation completed successfully")

if __name__ == "__main__":
    main()

