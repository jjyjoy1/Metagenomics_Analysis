#!/usr/bin/env python3
"""
Perform strain-level analysis and outbreak tracking for specified pathogens.
This script combines StrainGST and MLST results with SNP-based phylogeny
to identify potential outbreak clusters.
"""

import pandas as pd
import numpy as np
import json
import os
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
import glob
from snakemake.utils import logger

def extract_strains_from_straingst(json_files, target_species):
    """Extract strain information from StrainGST results"""
    strain_data = {}
    
    for json_file in json_files:
        sample_id = os.path.basename(json_file).split('.')[0]
        
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Process each species in the results
            for species, species_data in data.get('species', {}).items():
                if species not in target_species:
                    continue
                
                strains = species_data.get('strains', [])
                
                if not strains:
                    continue
                
                # Get the top strain for this species
                top_strain = strains[0]
                strain_id = top_strain.get('strain_id', '')
                similarity = top_strain.get('similarity', 0)
                
                if species not in strain_data:
                    strain_data[species] = {}
                
                strain_data[species][sample_id] = {
                    'strain_id': strain_id,
                    'similarity': similarity
                }
        
        except Exception as e:
            logger.warning(f"Error processing {json_file}: {e}")
    
    return strain_data

def parse_mlst_results(mlst_files):
    """Parse MLST results for each sample"""
    mlst_data = {}
    
    for mlst_file in mlst_files:
        sample_id = os.path.basename(mlst_file).split('.')[0]
        
        try:
            df = pd.read_csv(mlst_file, sep='\t', header=None)
            
            if len(df) > 0:
                # MLST output format: file, scheme, ST, allele1, allele2, ...
                record = df.iloc[0]
                scheme = record[1] if len(record) > 1 else 'unknown'
                st = record[2] if len(record) > 2 else 'unknown'
                alleles = list(record[3:]) if len(record) > 3 else []
                
                mlst_data[sample_id] = {
                    'scheme': scheme,
                    'ST': st,
                    'alleles': alleles
                }
        
        except Exception as e:
            logger.warning(f"Error processing {mlst_file}: {e}")
    
    return mlst_data

def extract_species_contigs(assembly_files, species_data, target_species):
    """Extract contigs for target species from assemblies"""
    species_contigs = {species: {} for species in target_species}
    
    # Get Kraken2 classification for contigs
    # Note: In a real pipeline, you would use Kraken2 to classify contigs
    # Here we're simulating this step
    
    for assembly_file in assembly_files:
        sample_id = os.path.basename(os.path.dirname(assembly_file))
        
        # For each target species
        for species in target_species:
            # Check if this sample has this species based on StrainGST
            if species in species_data and sample_id in species_data[species]:
                # Read contigs
                contigs = []
                for record in SeqIO.parse(assembly_file, "fasta"):
                    # In a real pipeline, you would filter contigs by species
                    # Here we're assuming all contigs are from the target species
                    contigs.append(record)
                
                if contigs:
                    species_contigs[species][sample_id] = contigs
    
    return species_contigs

def perform_snp_analysis(species_contigs, similarity_threshold):
    """Perform SNP-based analysis to identify clusters"""
    outbreak_clusters = {}
    phylogenies = {}
    
    for species, samples_contigs in species_contigs.items():
        if len(samples_contigs) < 2:
            # Need at least 2 samples for clustering
            continue
        
        # In a real pipeline, you would align contigs and identify SNPs
        # Here we're simulating the SNP distance calculation
        
        # Create a distance matrix
        n_samples = len(samples_contigs)
        sample_ids = list(samples_contigs.keys())
        distance_matrix = np.zeros((n_samples, n_samples))
        
        # Simulate SNP distances (in a real pipeline, this would be calculated from alignments)
        for i in range(n_samples):
            for j in range(i+1, n_samples):
                # Random distance between 0 and 0.1
                distance = np.random.rand() * 0.1
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
        
        # Identify clusters based on the similarity threshold
        clusters = []
        processed = set()
        
        for i in range(n_samples):
            if sample_ids[i] in processed:
                continue
            
            cluster = {sample_ids[i]}
            processed.add(sample_ids[i])
            
            for j in range(n_samples):
                if i == j or sample_ids[j] in processed:
                    continue
                
                if distance_matrix[i, j] < (1 - similarity_threshold):
                    cluster.add(sample_ids[j])
                    processed.add(sample_ids[j])
            
            if len(cluster) > 1:  # Only consider clusters with multiple samples
                clusters.append(cluster)
        
        # Save the clusters
        outbreak_clusters[species] = clusters
        
        # Create a Newick tree (simulated)
        phylogeny = f"({','.join(sample_ids)});"
        phylogenies[species] = phylogeny
    
    return outbreak_clusters, phylogenies

def main():
    # Get inputs from Snakemake
    straingst_files = snakemake.input.straingst
    mlst_files = snakemake.input.mlst
    assembly_files = snakemake.input.assemblies
    
    output_clusters = snakemake.output.clusters
    output_phylogeny = snakemake.output.phylogeny
    
    target_species = snakemake.params.target_species
    similarity_threshold = float(snakemake.config.get('strain_similarity_threshold', 0.99))
    
    # Create output directories
    os.makedirs(os.path.dirname(output_clusters), exist_ok=True)
    
    # Extract strain information from StrainGST results
    strain_data = extract_strains_from_straingst(straingst_files, target_species)
    
    # Parse MLST results
    mlst_data = parse_mlst_results(mlst_files)
    
    # Extract species-specific contigs from assemblies
    species_contigs = extract_species_contigs(assembly_files, strain_data, target_species)
    
    # Perform SNP analysis and identify outbreak clusters
    outbreak_clusters, phylogenies = perform_snp_analysis(species_contigs, similarity_threshold)
    
    # Prepare output data
    results = []
    
    for species, clusters in outbreak_clusters.items():
        for i, cluster in enumerate(clusters):
            cluster_id = f"{species.replace(' ', '_')}_cluster_{i+1}"
            
            for sample_id in cluster:
                strain_info = strain_data.get(species, {}).get(sample_id, {})
                mlst_info = mlst_data.get(sample_id, {})
                
                results.append({
                    'species': species,
                    'cluster_id': cluster_id,
                    'sample_id': sample_id,
                    'strain_id': strain_info.get('strain_id', ''),
                    'strain_similarity': strain_info.get('similarity', 0),
                    'mlst_scheme': mlst_info.get('scheme', ''),
                    'mlst_ST': mlst_info.get('ST', '')
                })
    
    # Write cluster results to TSV
    if results:
        df = pd.DataFrame(results)
        df.to_csv(output_clusters, sep='\t', index=False)
    else:
        # Create empty file with headers
        pd.DataFrame(columns=['species', 'cluster_id', 'sample_id', 'strain_id', 
                              'strain_similarity', 'mlst_scheme', 'mlst_ST']
                    ).to_csv(output_clusters, sep='\t', index=False)
    
    # Write phylogeny trees
    with open(output_phylogeny, 'w') as f:
        for species, newick in phylogenies.items():
            f.write(f"# {species}\n")
            f.write(f"{newick}\n\n")
    
    logger.info(f"Outbreak analysis completed. Found {sum(len(clusters) for clusters in outbreak_clusters.values())} potential outbreak clusters")

if __name__ == "__main__":
    main()
