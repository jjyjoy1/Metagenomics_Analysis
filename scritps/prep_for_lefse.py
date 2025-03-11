#!/usr/bin/env python3

"""
Script to prepare QIIME2 data for LEfSe analysis.
Converts feature table and taxonomy data to the format required by LEfSe.
"""

import argparse
import pandas as pd
import qiime2
import numpy as np
from biom import Table
from qiime2 import Artifact

def parse_args():
    parser = argparse.ArgumentParser(description='Prepare QIIME2 data for LEfSe analysis')
    parser.add_argument('--table', required=True, help='Path to QIIME2 feature table (.qza)')
    parser.add_argument('--taxonomy', required=True, help='Path to QIIME2 taxonomy file (.qza)')
    parser.add_argument('--metadata', required=True, help='Path to metadata file')
    parser.add_argument('--group-col', required=True, help='Column in metadata file for grouping')
    parser.add_argument('--output', required=True, help='Output file path for LEfSe input')
    return parser.parse_args()

def load_feature_table(table_path):
    """Load QIIME2 feature table artifact."""
    table_artifact = Artifact.load(table_path)
    feature_table = table_artifact.view(Table)
    df = feature_table.to_dataframe().T  # Transpose to have samples as rows
    return df

def load_taxonomy(taxonomy_path):
    """Load QIIME2 taxonomy artifact."""
    taxonomy_artifact = Artifact.load(taxonomy_path)
    taxonomy_df = taxonomy_artifact.view(pd.DataFrame)
    # Parse taxonomy strings
    taxonomy_series = taxonomy_df['Taxon'].apply(lambda x: pd.Series(x.split(';')))
    # Rename columns to taxonomic levels
    tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    taxonomy_series.columns = tax_levels[:len(taxonomy_series.columns)]
    # Clean up taxonomy names
    for col in taxonomy_series.columns:
        taxonomy_series[col] = taxonomy_series[col].str.strip()
    # Add feature IDs
    taxonomy_series.index = taxonomy_df.index
    return taxonomy_series

def load_metadata(metadata_path):
    """Load metadata file."""
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    return metadata

def prepare_lefse_input(feature_table, taxonomy, metadata, group_col):
    """Prepare data in LEfSe format."""
    # Get the group information from metadata
    groups = metadata[group_col]
    
    # Ensure samples match between feature table and metadata
    shared_samples = set(feature_table.index).intersection(set(groups.index))
    if not shared_samples:
        raise ValueError("No matching samples between feature table and metadata")
    
    feature_table = feature_table.loc[shared_samples]
    groups = groups.loc[shared_samples]
    
    # Create a merged dataframe for LEfSe input
    # First, normalize the feature table (relative abundance)
    norm_table = feature_table.div(feature_table.sum(axis=1), axis=0) * 100
    
    # Create LEfSe input
    lefse_df = pd.DataFrame(index=norm_table.index)
    lefse_df['class'] = groups
    
    # Add taxonomy and abundance information
    for feature_id in norm_table.columns:
        if feature_id in taxonomy.index:
            tax_string = '|'.join([level for level in taxonomy.loc[feature_id] if not pd.isna(level)])
            tax_string = tax_string.replace('[', '').replace(']', '')
            for sample in norm_table.index:
                lefse_df.loc[sample, tax_string] = norm_table.loc[sample, feature_id]
    
    return lefse_df

def main():
    args = parse_args()
    
    # Load data
    feature_table = load_feature_table(args.table)
    taxonomy = load_taxonomy(args.taxonomy)
    metadata = load_metadata(args.metadata)
    
    # Prepare LEfSe input
    lefse_input = prepare_lefse_input(feature_table, taxonomy, metadata, args.group_col)
    
    # Write to file in LEfSe format
    lefse_input.to_csv(args.output, sep='\t', index_label='sample_id')
    print(f"LEfSe input file created: {args.output}")

if __name__ == "__main__":
    main()

