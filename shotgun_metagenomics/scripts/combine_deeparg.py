#!/usr/bin/env python3
"""
Combine DeepARG results from multiple samples into a single abundance table.
"""

import pandas as pd
import glob
import os
import sys
from snakemake.utils import logger

def main():
    # Get input files from snakemake
    input_files = snakemake.input
    output_file = snakemake.output[0]
    
    # Create a dictionary to store the ARG abundance data
    arg_abundance = {}
    
    # Process each DeepARG output file
    for file_path in input_files:
        # Extract sample name from the filename
        sample_name = os.path.basename(file_path).split('.')[0]
        
        try:
            # Read the DeepARG output file
            df = pd.read_csv(file_path, sep='\t')
            
            # Group by ARG category and ARG type, and sum the abundances
            arg_summary = df.groupby(['ARG-category', 'ARG-type']).agg({'abundance-rpkm': 'sum'}).reset_index()
            
            # Create a unique identifier for each ARG (category + type)
            arg_summary['ARG'] = arg_summary['ARG-category'] + '|' + arg_summary['ARG-type']
            
            # Add the sample to the ARG abundance dictionary
            for _, row in arg_summary.iterrows():
                arg_id = row['ARG']
                if arg_id not in arg_abundance:
                    arg_abundance[arg_id] = {
                        'ARG-category': row['ARG-category'],
                        'ARG-type': row['ARG-type']
                    }
                
                # Add the abundance for this sample
                arg_abundance[arg_id][sample_name] = row['abundance-rpkm']
        
        except Exception as e:
            logger.warning(f"Error processing {file_path}: {e}")
    
    # Convert the dictionary to a DataFrame
    result_df = pd.DataFrame.from_dict(arg_abundance, orient='index')
    
    # Fill missing values with 0
    result_df = result_df.fillna(0)
    
    # Reset the index to make the ARG ID a column
    result_df = result_df.reset_index().rename(columns={'index': 'ARG'})
    
    # Rearrange columns to have ARG-category and ARG-type first, then samples
    first_cols = ['ARG', 'ARG-category', 'ARG-type']
    sample_cols = [col for col in result_df.columns if col not in first_cols]
    result_df = result_df[first_cols + sorted(sample_cols)]
    
    # Write the combined data to the output file
    result_df.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"Combined DeepARG results written to {output_file}")

if __name__ == "__main__":
    main()

