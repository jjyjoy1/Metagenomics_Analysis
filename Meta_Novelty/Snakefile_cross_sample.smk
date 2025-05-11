#############################################################################
# Part 6: Cross-Sample Analysis
#############################################################################

# Configuration for cross-sample analysis
METADATA = config.get("metadata", "metadata.csv")
CROSS_SAMPLE_DIR = OUTPUT_DIR + "/cross_sample"

# Checkpoint to ensure all individual samples are processed
checkpoint all_samples_processed:
    input:
        # Key outputs from individual sample processing
        expand(OUTPUT_DIR + "/genes/{sample}_proteins.faa", sample=SAMPLES),
        expand(OUTPUT_DIR + "/cat/{sample}_cat.taxonomy.tsv", sample=SAMPLES),
        expand(OUTPUT_DIR + "/graphbin/{sample}_graphbin_bins.csv", sample=SAMPLES)
    output:
        touch(CROSS_SAMPLE_DIR + "/all_samples_processed.flag")
    shell:
        """
        mkdir -p {CROSS_SAMPLE_DIR}
        """

#############################################################################
# Core Matrix 1: Pan-genome Matrix Generation
#############################################################################

# Collect and cluster genes from all samples
rule collect_genes:
    input:
        checkpoint_flag = CROSS_SAMPLE_DIR + "/all_samples_processed.flag",
        gene_files = expand(OUTPUT_DIR + "/genes/{sample}_genes.fna", sample=SAMPLES)
    output:
        combined_genes = CROSS_SAMPLE_DIR + "/all_genes.fna",
        sample_map = CROSS_SAMPLE_DIR + "/gene_sample_map.tsv"
    shell:
        """
        # Combine all gene files with sample identifiers
        touch {output.combined_genes}
        touch {output.sample_map}
        
        for gene_file in {input.gene_files}; do
            sample=$(basename "$gene_file" | sed 's/_genes.fna//')
            
            # Add sample identifier to gene headers and append to combined file
            python -c "
from Bio import SeqIO
import os

sample = '{sample}'
gene_file = '{gene_file}'
output_file = '{output.combined_genes}'
sample_map = '{output.sample_map}'

with open(sample_map, 'a') as map_handle:
    for record in SeqIO.parse(gene_file, 'fasta'):
        original_id = record.id
        record.id = f'{{sample}}__{{record.id}}'
        record.description = f'{{sample}}__{{record.description}}'
        
        # Write to combined gene file
        with open(output_file, 'a') as out_handle:
            SeqIO.write(record, out_handle, 'fasta')
        
        # Write mapping information
        map_handle.write(f'{{record.id}}\\t{{sample}}\\t{{original_id}}\\n')
"
        done
        """

# Cluster genes to create pan-genome
rule cluster_genes:
    input:
        combined_genes = rules.collect_genes.output.combined_genes,
        sample_map = rules.collect_genes.output.sample_map
    output:
        clusters = CROSS_SAMPLE_DIR + "/gene_clusters.tsv",
        representative_genes = CROSS_SAMPLE_DIR + "/representative_genes.faa"
    params:
        identity = config.get("gene_cluster_identity", 0.95),
        coverage = config.get("gene_cluster_coverage", 0.9)
    threads: THREADS
    shell:
        """
        # Cluster genes using CD-HIT
        cd-hit-est -i {input.combined_genes} -o {output.representative_genes} \
            -c {params.identity} -aS {params.coverage} -M 0 -T {threads} \
            -d 0 -g 1 > {CROSS_SAMPLE_DIR}/cd-hit.log
        
        # Convert CD-HIT output to cluster file
        python -c "
import re

cluster_file = '{output.representative_genes}.clstr'
output_file = '{output.clusters}'

clusters = {{}}
with open(cluster_file, 'r') as f:
    current_cluster = None
    for line in f:
        if line.startswith('>Cluster'):
            current_cluster = line.strip().split()[1]
            clusters[current_cluster] = []
        else:
            gene_match = re.search(r'>([^.]+)', line)
            if gene_match:
                gene_id = gene_match.group(1)
                clusters[current_cluster].append(gene_id)

with open(output_file, 'w') as out:
    out.write('cluster_id\\tgene_id\\n')
    for cluster_id, genes in clusters.items():
        for gene in genes:
            out.write(f'{{cluster_id}}\\t{{gene}}\\n')
"
        """

# Generate pan-genome matrix
rule generate_pangenome_matrix:
    input:
        clusters = rules.cluster_genes.output.clusters,
        sample_map = rules.collect_genes.output.sample_map
    output:
        matrix = CROSS_SAMPLE_DIR + "/pan_genome_matrix.tsv",
        core_genes = CROSS_SAMPLE_DIR + "/core_genes.txt",
        accessory_genes = CROSS_SAMPLE_DIR + "/accessory_genes.txt",
        unique_genes = CROSS_SAMPLE_DIR + "/unique_genes.txt"
    params:
        core_threshold = config.get("core_gene_threshold", 0.95)  # Percentage of samples required for core genes
    shell:
        """
        python -c "
import pandas as pd
import numpy as np

# Load sample mapping
sample_map = pd.read_csv('{input.sample_map}', sep='\\t', 
                         names=['gene_id', 'sample', 'original_id'])
sample_list = sorted(sample_map['sample'].unique())

# Load cluster information
clusters = pd.read_csv('{input.clusters}', sep='\\t')

# Create mapping from gene to cluster
gene_to_cluster = dict(zip(clusters['gene_id'], clusters['cluster_id']))

# Prepare data for the matrix
matrix_data = []
for sample in sample_list:
    # Get genes for this sample
    sample_genes = set(sample_map[sample_map['sample'] == sample]['gene_id'])
    
    # Get clusters for these genes
    sample_clusters = [gene_to_cluster[gene] for gene in sample_genes if gene in gene_to_cluster]
    
    # Count genes per cluster for this sample (should be 0 or 1 for most cases)
    cluster_counts = pd.Series(sample_clusters).value_counts()
    
    # Create row with sample and cluster information
    row = {'sample': sample}
    
    # Add cluster columns (1 for present, 0 for absent)
    unique_clusters = sorted(clusters['cluster_id'].unique())
    for cluster in unique_clusters:
        row[f'cluster_{{cluster}}'] = 1 if cluster in cluster_counts.index else 0
    
    matrix_data.append(row)

# Create the matrix DataFrame
matrix_df = pd.DataFrame(matrix_data)
matrix_df.set_index('sample', inplace=True)

# Save the matrix
matrix_df.to_csv('{output.matrix}', sep='\\t')

# Identify core, accessory, and unique genes
sample_count = len(sample_list)
core_threshold = {params.core_threshold} * sample_count
cluster_presence = matrix_df.sum(axis=0)

core_clusters = cluster_presence[cluster_presence >= core_threshold].index.tolist()
unique_clusters = cluster_presence[cluster_presence == 1].index.tolist()
accessory_clusters = cluster_presence[(cluster_presence > 1) & (cluster_presence < core_threshold)].index.tolist()

# Strip 'cluster_' prefix
core_clusters = [c.replace('cluster_', '') for c in core_clusters]
unique_clusters = [c.replace('cluster_', '') for c in unique_clusters]
accessory_clusters = [c.replace('cluster_', '') for c in accessory_clusters]

# Save core, accessory, and unique gene lists
with open('{output.core_genes}', 'w') as f:
    for cluster in core_clusters:
        f.write(f'{{cluster}}\\n')

with open('{output.accessory_genes}', 'w') as f:
    for cluster in accessory_clusters:
        f.write(f'{{cluster}}\\n')

with open('{output.unique_genes}', 'w') as f:
    for cluster in unique_clusters:
        f.write(f'{{cluster}}\\n')

print(f'Total clusters: {{len(cluster_presence)}}')
print(f'Core gene clusters: {{len(core_clusters)}}')
print(f'Accessory gene clusters: {{len(accessory_clusters)}}')
print(f'Unique gene clusters: {{len(unique_clusters)}}')
"
        """

#############################################################################
# Core Matrix 2: Functional Profile Matrix Generation
#############################################################################

# Add eggNOG-mapper annotation for KO terms and pathways
rule eggnog_annotation:
    input:
        proteins = OUTPUT_DIR + "/genes/{sample}_proteins.faa"
    output:
        annotations = OUTPUT_DIR + "/functional/{sample}_eggnog.emapper.annotations",
        orthologs = OUTPUT_DIR + "/functional/{sample}_eggnog.emapper.seed_orthologs"
    params:
        eggnog_db = config.get("eggnog_db", "/path/to/eggnog_db"),
        output_prefix = OUTPUT_DIR + "/functional/{sample}_eggnog"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTPUT_DIR}/functional
        
        # Run eggNOG-mapper
        emapper.py -i {input.proteins} \
            --output {params.output_prefix} \
            -m diamond \
            --data_dir {params.eggnog_db} \
            --cpu {threads} \
            --go_evidence non-electronic \
            --target_orthologs all \
            --seed_ortholog_evalue 0.001 \
            --seed_ortholog_score 60 \
            --tax_scope auto \
            --annotation_tax_scope auto \
            --override
        """

# Extract and compile KO annotations across samples
rule extract_ko_annotations:
    input:
        annotations = OUTPUT_DIR + "/functional/{sample}_eggnog.emapper.annotations"
    output:
        ko_file = OUTPUT_DIR + "/functional/{sample}_ko_counts.tsv"
    shell:
        """
        # Extract KO terms and count occurrences
        python -c "
import pandas as pd
import re

# Read eggNOG annotations (skip comments at the beginning)
with open('{input.annotations}', 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            header_line = i - 1
            break

# Read the actual data
annot_df = pd.read_csv('{input.annotations}', sep='\\t', skiprows=header_line)

# KEGG KOs are in the KEGG_ko column, format: ko:K00001,ko:K00002,...
ko_counts = {{}}

if 'KEGG_ko' in annot_df.columns:
    for ko_str in annot_df['KEGG_ko'].dropna():
        kos = re.findall(r'ko:(K\\d+)', ko_str)
        for ko in kos:
            ko_counts[ko] = ko_counts.get(ko, 0) + 1

# Write KO counts to file
with open('{output.ko_file}', 'w') as f:
    f.write('KO\\tcount\\n')
    for ko, count in sorted(ko_counts.items()):
        f.write(f'{{ko}}\\t{{count}}\\n')
"
        """

# Combine KO counts across samples
rule combine_ko_counts:
    input:
        checkpoint_flag = CROSS_SAMPLE_DIR + "/all_samples_processed.flag",
        ko_files = expand(OUTPUT_DIR + "/functional/{sample}_ko_counts.tsv", sample=SAMPLES)
    output:
        ko_matrix = CROSS_SAMPLE_DIR + "/ko_abundance_matrix.tsv"
    shell:
        """
        # Combine KO counts from all samples
        python -c "
import pandas as pd
import glob
import os

# Get all KO files
ko_files = '{input.ko_files}'.split()

# Combine KO counts
ko_matrix = {{}}

for ko_file in ko_files:
    sample = os.path.basename(ko_file).replace('_ko_counts.tsv', '')
    ko_df = pd.read_csv(ko_file, sep='\\t')
    
    # Add sample as a column
    ko_matrix[sample] = pd.Series(ko_df['count'].values, index=ko_df['KO'])

# Convert to DataFrame
ko_matrix_df = pd.DataFrame(ko_matrix).fillna(0)

# Save the matrix
ko_matrix_df.to_csv('{output.ko_matrix}', sep='\\t')
"
        """

# Extract and compile pathway annotations
rule extract_pathway_annotations:
    input:
        annotations = OUTPUT_DIR + "/functional/{sample}_eggnog.emapper.annotations"
    output:
        pathway_file = OUTPUT_DIR + "/functional/{sample}_pathway_counts.tsv"
    shell:
        """
        # Extract KEGG pathway information and count occurrences
        python -c "
import pandas as pd
import re

# Read eggNOG annotations (skip comments at the beginning)
with open('{input.annotations}', 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            header_line = i - 1
            break

# Read the actual data
annot_df = pd.read_csv('{input.annotations}', sep='\\t', skiprows=header_line)

# KEGG pathways are in the KEGG_Pathway column, format: map00010,map00020,...
pathway_counts = {{}}

if 'KEGG_Pathway' in annot_df.columns:
    for pathway_str in annot_df['KEGG_Pathway'].dropna():
        pathways = re.findall(r'(map\\d+)', pathway_str)
        for pathway in pathways:
            pathway_counts[pathway] = pathway_counts.get(pathway, 0) + 1

# Write pathway counts to file
with open('{output.pathway_file}', 'w') as f:
    f.write('Pathway\\tcount\\n')
    for pathway, count in sorted(pathway_counts.items()):
        f.write(f'{{pathway}}\\t{{count}}\\n')
"
        """

# Combine pathway counts across samples
rule combine_pathway_counts:
    input:
        checkpoint_flag = CROSS_SAMPLE_DIR + "/all_samples_processed.flag",
        pathway_files = expand(OUTPUT_DIR + "/functional/{sample}_pathway_counts.tsv", sample=SAMPLES)
    output:
        pathway_matrix = CROSS_SAMPLE_DIR + "/pathway_abundance_matrix.tsv"
    shell:
        """
        # Combine pathway counts from all samples
        python -c "
import pandas as pd
import glob
import os

# Get all pathway files
pathway_files = '{input.pathway_files}'.split()

# Combine pathway counts
pathway_matrix = {{}}

for pathway_file in pathway_files:
    sample = os.path.basename(pathway_file).replace('_pathway_counts.tsv', '')
    pathway_df = pd.read_csv(pathway_file, sep='\\t')
    
    # Add sample as a column
    pathway_matrix[sample] = pd.Series(pathway_df['count'].values, index=pathway_df['Pathway'])

# Convert to DataFrame
pathway_matrix_df = pd.DataFrame(pathway_matrix).fillna(0)

# Save the matrix
pathway_matrix_df.to_csv('{output.pathway_matrix}', sep='\\t')
"
        """

#############################################################################
# Core Matrix 3: Taxonomic Abundance Matrix Generation
#############################################################################

# Process CAT taxonomy results for each sample
rule process_cat_taxonomy:
    input:
        taxonomy = OUTPUT_DIR + "/cat/{sample}_cat.taxonomy.tsv",
        coverage = OUTPUT_DIR + "/mapping/{sample}_coverage.txt"
    output:
        processed = OUTPUT_DIR + "/taxonomy/{sample}_processed_taxonomy.tsv"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/taxonomy
        
        # Process CAT taxonomy results with coverage information
        python -c "
import pandas as pd
import numpy as np

# Read taxonomy file (skip comment lines)
taxonomy_df = pd.read_csv('{input.taxonomy}', sep='\\t', comment='#')
taxonomy_df = taxonomy_df.rename(columns={{'# contig': 'contig_id'}})

# Read coverage file
coverage_df = pd.read_csv('{input.coverage}', sep=' ', names=['contig_id', 'coverage'])

# Merge taxonomy and coverage information
merged_df = pd.merge(taxonomy_df, coverage_df, on='contig_id', how='left')
merged_df['coverage'] = merged_df['coverage'].fillna(0)

# Extract taxonomic ranks (superkingdom, phylum, class, order, family, genus, species)
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
for rank in ranks:
    merged_df[rank] = merged_df['classification'].apply(
        lambda x: next((t.split(':')[1] for t in x.split(';') if t.startswith(rank)), 'Unknown')
        if isinstance(x, str) else 'Unknown'
    )

# Save processed taxonomy
merged_df.to_csv('{output.processed}', sep='\\t', index=False)
"
        """

# Combine taxonomic information across samples
rule combine_taxonomy:
    input:
        checkpoint_flag = CROSS_SAMPLE_DIR + "/all_samples_processed.flag",
        taxonomy_files = expand(OUTPUT_DIR + "/taxonomy/{sample}_processed_taxonomy.tsv", sample=SAMPLES)
    output:
        phylum_matrix = CROSS_SAMPLE_DIR + "/phylum_abundance_matrix.tsv",
        genus_matrix = CROSS_SAMPLE_DIR + "/genus_abundance_matrix.tsv",
        species_matrix = CROSS_SAMPLE_DIR + "/species_abundance_matrix.tsv"
    shell:
        """
        # Combine taxonomy information from all samples
        python -c "
import pandas as pd
import numpy as np
import os

# Get all taxonomy files
taxonomy_files = '{input.taxonomy_files}'.split()

# Initialize abundance matrices
phylum_abundance = {{}}
genus_abundance = {{}}
species_abundance = {{}}

# Process each sample
for tax_file in taxonomy_files:
    sample = os.path.basename(tax_file).replace('_processed_taxonomy.tsv', '')
    tax_df = pd.read_csv(tax_file, sep='\\t')
    
    # Calculate abundances at different taxonomic levels
    # Weight by coverage - more coverage means more abundance
    
    # Phylum level
    phylum_counts = tax_df.groupby('phylum')['coverage'].sum()
    total_coverage = phylum_counts.sum()
    if total_coverage > 0:  # Avoid division by zero
        phylum_rel_abundance = phylum_counts / total_coverage
    else:
        phylum_rel_abundance = phylum_counts  # All zeros
    
    phylum_abundance[sample] = phylum_rel_abundance
    
    # Genus level
    genus_counts = tax_df.groupby('genus')['coverage'].sum()
    if total_coverage > 0:  # Avoid division by zero
        genus_rel_abundance = genus_counts / total_coverage
    else:
        genus_rel_abundance = genus_counts  # All zeros
    
    genus_abundance[sample] = genus_rel_abundance
    
    # Species level
    species_counts = tax_df.groupby('species')['coverage'].sum()
    if total_coverage > 0:  # Avoid division by zero
        species_rel_abundance = species_counts / total_coverage
    else:
        species_rel_abundance = species_counts  # All zeros
    
    species_abundance[sample] = species_rel_abundance

# Convert to DataFrames
phylum_df = pd.DataFrame(phylum_abundance).fillna(0)
genus_df = pd.DataFrame(genus_abundance).fillna(0)
species_df = pd.DataFrame(species_abundance).fillna(0)

# Save matrices
phylum_df.to_csv('{output.phylum_matrix}', sep='\\t')
genus_df.to_csv('{output.genus_matrix}', sep='\\t')
species_df.to_csv('{output.species_matrix}', sep='\\t')
"
        """

# Process bin-level taxonomy
rule process_bin_taxonomy:
    input:
        bat_taxonomy = OUTPUT_DIR + "/bat/{sample}_bat.taxonomy.tsv",
        graphbin_csv = OUTPUT_DIR + "/graphbin/{sample}_graphbin_bins.csv",
        coverage = OUTPUT_DIR + "/mapping/{sample}_coverage.txt"
    output:
        bin_taxonomy = OUTPUT_DIR + "/taxonomy/{sample}_bin_taxonomy.tsv"
    shell:
        """
        # Process BAT taxonomy results at bin level
        python -c "
import pandas as pd
import numpy as np

# Read BAT taxonomy file
try:
    bat_df = pd.read_csv('{input.bat_taxonomy}', sep='\\t', comment='#')
    # Handle different column naming in BAT output
    if '# bin' in bat_df.columns:
        bat_df = bat_df.rename(columns={{'# bin': 'bin_id'}})
    else:
        bat_df['bin_id'] = bat_df.iloc[:, 0]  # Assume first column is bin ID
except Exception as e:
    print(f'Error reading BAT taxonomy: {{e}}')
    bat_df = pd.DataFrame(columns=['bin_id', 'classification'])

# Read GraphBin contig to bin assignments
try:
    graphbin_df = pd.read_csv('{input.graphbin_csv}', header=0)
    graphbin_df.columns = ['contig_id', 'bin_id']
except Exception as e:
    print(f'Error reading GraphBin assignments: {{e}}')
    graphbin_df = pd.DataFrame(columns=['contig_id', 'bin_id'])

# Read coverage information
try:
    coverage_df = pd.read_csv('{input.coverage}', sep=' ', names=['contig_id', 'coverage'])
except Exception as e:
    print(f'Error reading coverage file: {{e}}')
    coverage_df = pd.DataFrame(columns=['contig_id', 'coverage'])

# Calculate bin coverage as average of contig coverages
bin_coverage = {{}}
for bin_id in graphbin_df['bin_id'].unique():
    if bin_id == 'unbinned':
        continue
    
    bin_contigs = graphbin_df[graphbin_df['bin_id'] == bin_id]['contig_id']
    contig_coverages = coverage_df[coverage_df['contig_id'].isin(bin_contigs)]['coverage']
    
    if not contig_coverages.empty:
        bin_coverage[bin_id] = contig_coverages.mean()
    else:
        bin_coverage[bin_id] = 0

# Extract taxonomic ranks for bins
result_df = pd.DataFrame()
for bin_id, coverage in bin_coverage.items():
    # Find bin in BAT results
    bin_taxonomy = bat_df[bat_df['bin_id'] == bin_id]['classification'].values
    
    if len(bin_taxonomy) > 0:
        classification = bin_taxonomy[0]
    else:
        classification = 'Unknown'
    
    # Extract taxonomic ranks
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    rank_values = {{}}
    
    for rank in ranks:
        if isinstance(classification, str):
            rank_values[rank] = next((t.split(':')[1] for t in classification.split(';') 
                                     if t.startswith(rank)), 'Unknown')
        else:
            rank_values[rank] = 'Unknown'
    
    # Create row
    row = {{
        'bin_id': bin_id,
        'coverage': coverage,
        'classification': classification
    }}
    row.update(rank_values)
    
    # Append to result
    result_df = pd.concat([result_df, pd.DataFrame([row])], ignore_index=True)

# Save bin taxonomy
if not result_df.empty:
    result_df.to_csv('{output.bin_taxonomy}', sep='\\t', index=False)
else:
    # Create empty file with headers
    with open('{output.bin_taxonomy}', 'w') as f:
        f.write('bin_id\\tcoverage\\tclassification\\tsuperkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n')
"
        """

# Combine bin-level taxonomy across samples
rule combine_bin_taxonomy:
    input:
        checkpoint_flag = CROSS_SAMPLE_DIR + "/all_samples_processed.flag",
        bin_taxonomy_files = expand(OUTPUT_DIR + "/taxonomy/{sample}_bin_taxonomy.tsv", sample=SAMPLES)
    output:
        bin_matrix = CROSS_SAMPLE_DIR + "/bin_taxonomy_matrix.tsv"
    shell:
        """
        # Combine bin-level taxonomy from all samples
        python -c "
import pandas as pd
import os

# Get all bin taxonomy files
bin_files = '{input.bin_taxonomy_files}'.split()

# Combine bin taxonomies
all_bins = []

for bin_file in bin_files:
    sample = os.path.basename(bin_file).replace('_bin_taxonomy.tsv', '')
    
    try:
        bin_df = pd.read_csv(bin_file, sep='\\t')
        bin_df['sample'] = sample
        all_bins.append(bin_df)
    except Exception as e:
        print(f'Error processing {{bin_file}}: {{e}}')

if all_bins:
    # Combine all bin data
    combined_df = pd.concat(all_bins, ignore_index=True)
    
    # Create a unique identifier for each bin across samples
    combined_df['bin_identifier'] = combined_df['sample'] + '--' + combined_df['bin_id']
    
    # Set index for the final matrix
    combined_df.set_index('bin_identifier', inplace=True)
    
    # Save the combined matrix
    combined_df.to_csv('{output.bin_matrix}', sep='\\t')
else:
    # Create empty file with basic headers
    with open('{output.bin_matrix}', 'w') as f:
        f.write('bin_identifier\\tsample\\tbin_id\\tcoverage\\tclassification\\tsuperkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n')
"
        """

# Integrate cross-sample matrices with metadata
rule integrate_with_metadata:
    input:
        pangenome_matrix = CROSS_SAMPLE_DIR + "/pan_genome_matrix.tsv",
        ko_matrix = CROSS_SAMPLE_DIR + "/ko_abundance_matrix.tsv",
        pathway_matrix = CROSS_SAMPLE_DIR + "/pathway_abundance_matrix.tsv",
        phylum_matrix = CROSS_SAMPLE_DIR + "/phylum_abundance_matrix.tsv",
        genus_matrix = CROSS_SAMPLE_DIR + "/genus_abundance_matrix.tsv",
        species_matrix = CROSS_SAMPLE_DIR + "/species_abundance_matrix.tsv",
        bin_matrix = CROSS_SAMPLE_DIR + "/bin_taxonomy_matrix.tsv",
        metadata = METADATA
    output:
        integrated_metadata = CROSS_SAMPLE_DIR + "/integrated_metadata.tsv",
        matrix_stats = CROSS_SAMPLE_DIR + "/matrix_statistics.tsv"
    shell:
        """
        # Integrate all matrices with metadata for further analysis
        python -c "
import pandas as pd
import os
import numpy as np

# Load metadata
try:
    metadata = pd.read_csv('{input.metadata}')
    # Ensure sample_id column exists
    if 'sample_id' not in metadata.columns:
        # Try to determine which column contains sample IDs
        # This is a heuristic and might need adjustment
        potential_id_cols = [col for col in metadata.columns if 'id' in col.lower() or 'sample' in col.lower()]
        if potential_id_cols:
            metadata = metadata.rename(columns={{potential_id_cols[0]: 'sample_id'}})
        else:
            # Assume first column is sample ID
            metadata = metadata.rename(columns={{metadata.columns[0]: 'sample_id'}})
    
    # Ensure sample_id is the index
    metadata = metadata.set_index('sample_id')
except Exception as e:
    print(f'Error loading metadata: {{e}}')
    # Create empty metadata with samples from matrices
    metadata = pd.DataFrame(index=pd.read_csv('{input.pangenome_matrix}', sep='\\t', index_col=0).index)

# Collect matrix statistics
matrix_stats = {{
    'matrix_name': [],
    'rows': [],
    'columns': [],
    'sparsity': [],
    'non_zero_entries': []
}}

# Function to calculate matrix statistics
def add_matrix_stats(name, matrix_df):
    total_cells = matrix_df.shape[0] * matrix_df.shape[1]
    non_zero = np.count_nonzero(matrix_df.values)
    sparsity = 1 - (non_zero / total_cells) if total_cells > 0 else 1
    
    matrix_stats['matrix_name'].append(name)
    matrix_stats['rows'].append(matrix_df.shape[0])
    matrix_stats['columns'].append(matrix_df.shape[1])
    matrix_stats['sparsity'].append(sparsity)
    matrix_stats['non_zero_entries'].append(non_zero)

# Process each matrix
matrices = [
    ('Pan-genome', '{input.pangenome_matrix}'),
    ('KO Abundance', '{input.ko_matrix}'),
    ('Pathway Abundance', '{input.pathway_matrix}'),
    ('Phylum Abundance', '{input.phylum_matrix}'),
    ('Genus Abundance', '{input.genus_matrix}'),
    ('Species Abundance', '{input.species_matrix}')
]

matrix_features = {{}}

for name, path in matrices:
    try:
        df = pd.read_csv(path, sep='\\t', index_col=0)
        
        # Collect basic stats
        add_matrix_stats(name, df)
        
        # For taxonomic matrices, calculate diversity metrics
        if name in ['Phylum Abundance', 'Genus Abundance', 'Species Abundance']:
            # Shannon diversity
            for sample in df.columns:
                abundances = df[sample].values
                abundances = abundances[abundances > 0]  # Only consider non-zero
                if len(abundances) > 0:
                    shannon = -np.sum(abundances * np.log(abundances))
                    matrix_features[f'{{name}}_shannon_{{sample}}'] = shannon
                else:
                    matrix_features[f'{{name}}_shannon_{{sample}}'] = 0
        
        # For functional matrices, calculate richness
        if name in ['KO Abundance', 'Pathway Abundance']:
            for sample in df.columns:
                richness = np.sum(df[sample] > 0)
                matrix_features[f'{{name}}_richness_{{sample}}'] = richness
   
        # For pangenome, calculate core/accessory/unique proportions
        if name == 'Pan-genome':
            # Count clusters per sample
            for sample in df.index:
                present_genes = np.sum(df.loc[sample] > 0)
                total_genes = len(df.columns)
                matrix_features[f'pangenome_gene_count_{sample}'] = present_genes
                matrix_features[f'pangenome_gene_proportion_{sample}'] = present_genes / total_genes if total_genes > 0 else 0
    except Exception as e:
        print(f'Error processing {name} matrix: {e}')

# Process bin matrix separately due to different structure
try:
    bin_df = pd.read_csv('{input.bin_matrix}', sep='\\t')
    bin_counts = bin_df.groupby('sample')['bin_id'].count()
    
    for sample, count in bin_counts.items():
        matrix_features[f'bin_count_{sample}'] = count
        
    # Calculate taxonomic diversity of bins
    for sample in bin_df['sample'].unique():
        sample_bins = bin_df[bin_df['sample'] == sample]
        phylum_diversity = len(sample_bins['phylum'].unique())
        genus_diversity = len(sample_bins['genus'].unique())
        matrix_features[f'bin_phylum_diversity_{sample}'] = phylum_diversity
        matrix_features[f'bin_genus_diversity_{sample}'] = genus_diversity
except Exception as e:
    print(f'Error processing bin matrix: {e}')

# Create feature DataFrame
feature_df = pd.DataFrame()
for feature, value in matrix_features.items():
    # Extract sample from feature name
    parts = feature.split('_')
    sample = parts[-1]
    feature_type = '_'.join(parts[:-1])
    
    # Add to DataFrame
    if sample not in feature_df.index:
        feature_df.loc[sample, feature_type] = value
    else:
        feature_df.loc[sample, feature_type] = value

# Merge features with metadata
if not feature_df.empty:
    # Ensure sample IDs match
    common_samples = set(metadata.index).intersection(feature_df.index)
    
    if common_samples:
        integrated_df = pd.merge(
            metadata.loc[list(common_samples)], 
            feature_df.loc[list(common_samples)],
            left_index=True, right_index=True
        )
    else:
        # No common samples, just concatenate
        integrated_df = pd.concat([metadata, feature_df], axis=1)
else:
    integrated_df = metadata.copy()

# Save integrated metadata
integrated_df.to_csv('{output.integrated_metadata}', sep='\\t')

# Save matrix statistics
pd.DataFrame(matrix_stats).to_csv('{output.matrix_stats}', sep='\\t', index=False)
"
        """

# Create basic visualizations of the matrices
rule visualize_matrices:
    input:
        pangenome_matrix = CROSS_SAMPLE_DIR + "/pan_genome_matrix.tsv",
        ko_matrix = CROSS_SAMPLE_DIR + "/ko_abundance_matrix.tsv",
        phylum_matrix = CROSS_SAMPLE_DIR + "/phylum_abundance_matrix.tsv",
        species_matrix = CROSS_SAMPLE_DIR + "/species_abundance_matrix.tsv",
        metadata = METADATA,
        integrated_metadata = CROSS_SAMPLE_DIR + "/integrated_metadata.tsv"
    output:
        heatmap_dir = directory(CROSS_SAMPLE_DIR + "/visualizations/heatmaps"),
        pca_dir = directory(CROSS_SAMPLE_DIR + "/visualizations/pca"),
        summary_plots = directory(CROSS_SAMPLE_DIR + "/visualizations/summary")
    shell:
        """
        mkdir -p {output.heatmap_dir}
        mkdir -p {output.pca_dir}
        mkdir -p {output.summary_plots}
        
        # Create visualizations using Python with matplotlib
        python -c "
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os

# Set plot style
sns.set(style='whitegrid')
plt.rcParams['figure.figsize'] = (12, 8)

# Function to create heatmap
def create_heatmap(matrix_path, output_path, title, max_features=50):
    try:
        # Load data
        df = pd.read_csv(matrix_path, sep='\\t', index_col=0)
        
        # If too many features, select most variable
        if df.shape[1] > max_features:
            variances = df.var(axis=0)
            top_features = variances.sort_values(ascending=False).head(max_features).index
            df = df[top_features]
        
        # Create heatmap
        plt.figure(figsize=(max(10, df.shape[1]/4), max(8, df.shape[0]/4)))
        sns.heatmap(df, cmap='viridis', xticklabels=True, yticklabels=True)
        plt.title(title)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()
    except Exception as e:
        print(f'Error creating heatmap for {{matrix_path}}: {{e}}')

# Function to perform PCA and visualize
def create_pca_plot(matrix_path, metadata_path, output_path, title, max_features=1000):
    try:
        # Load data
        df = pd.read_csv(matrix_path, sep='\\t', index_col=0)
        metadata = pd.read_csv(metadata_path, sep='\\t', index_col=0)
        
        # If too many features, select most variable
        if df.shape[1] > max_features:
            variances = df.var(axis=0)
            top_features = variances.sort_values(ascending=False).head(max_features).index
            df = df[top_features]
        
        # Standardize the data
        X = StandardScaler().fit_transform(df)
        
        # Apply PCA
        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(X)
        
        # Create DataFrame with principal components
        pca_df = pd.DataFrame(data=principal_components, 
                             columns=['PC1', 'PC2'], 
                             index=df.index)
        
        # Try to add metadata coloring
        plt.figure(figsize=(10, 8))
        
        # Check for categorical columns in metadata that match samples
        common_samples = set(pca_df.index).intersection(metadata.index)
        if common_samples:
            pca_df = pca_df.loc[list(common_samples)]
            metadata = metadata.loc[list(common_samples)]
            
            # Find categorical columns with few unique values
            categorical_cols = []
            for col in metadata.columns:
                if metadata[col].dtype == 'object' and metadata[col].nunique() <= 10:
                    categorical_cols.append(col)
            
            if categorical_cols:
                # Use first categorical column for coloring
                col = categorical_cols[0]
                categories = metadata[col].astype(str)
                
                # Plot with categorical coloring
                for category in categories.unique():
                    subset = pca_df[categories == category]
                    plt.scatter(subset['PC1'], subset['PC2'], label=category, alpha=0.7)
                
                plt.legend(title=col)
            else:
                # No useful categorical column, just plot points
                plt.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.7)
        else:
            # No matching samples, just plot points
            plt.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.7)
        
        # Add axis labels with explained variance
        explained_variance = pca.explained_variance_ratio_
        plt.xlabel(f'PC1 ({{explained_variance[0]:.2%}} explained variance)')
        plt.ylabel(f'PC2 ({{explained_variance[1]:.2%}} explained variance)')
        
        plt.title(title)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()
    except Exception as e:
        print(f'Error creating PCA plot for {{matrix_path}}: {{e}}')

# Function to create summary plots
def create_summary_plots(matrix_path, output_dir, matrix_type):
    try:
        # Load data
        df = pd.read_csv(matrix_path, sep='\\t', index_col=0)
        
        if matrix_type == 'taxonomy':
            # Plot relative abundance of top 10 taxa
            abundances = df.mean(axis=1).sort_values(ascending=False)
            top_taxa = abundances.head(10)
            
            plt.figure(figsize=(12, 6))
            top_taxa.plot(kind='bar')
            plt.title(f'Top 10 Most Abundant Taxa')
            plt.ylabel('Mean Relative Abundance')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'top_taxa_abundance.png'), dpi=300)
            plt.close()
            
            # Plot stacked bar chart of sample composition
            top_taxa_df = df.loc[top_taxa.index]
            top_taxa_df = top_taxa_df.append(pd.Series(1 - top_taxa_df.sum(), name='Other'))
            
            plt.figure(figsize=(15, 8))
            top_taxa_df.T.plot(kind='bar', stacked=True, colormap='tab20')
            plt.title('Taxonomic Composition per Sample')
            plt.xlabel('Sample')
            plt.ylabel('Relative Abundance')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'taxonomic_composition.png'), dpi=300)
            plt.close()
            
        elif matrix_type == 'functional':
            # Count features per sample
            features_per_sample = (df > 0).sum(axis=0)
            
            plt.figure(figsize=(12, 6))
            features_per_sample.plot(kind='bar')
            plt.title('Functional Richness per Sample')
            plt.ylabel('Number of Features')
            plt.xlabel('Sample')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'functional_richness.png'), dpi=300)
            plt.close()
            
            # Plot top 10 features
            mean_abundances = df.mean(axis=1).sort_values(ascending=False)
            top_features = mean_abundances.head(10)
            
            plt.figure(figsize=(12, 6))
            top_features.plot(kind='bar')
            plt.title(f'Top 10 Most Abundant Functional Features')
            plt.ylabel('Mean Abundance')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'top_functional_features.png'), dpi=300)
            plt.close()
            
        elif matrix_type == 'pangenome':
            # Calculate gene presence distribution
            gene_presence = df.sum(axis=1)
            gene_count_bins = np.linspace(0, gene_presence.max(), 20)
            
            plt.figure(figsize=(12, 6))
            plt.hist(gene_presence, bins=gene_count_bins)
            plt.title('Distribution of Gene Counts per Sample')
            plt.xlabel('Number of Genes Present')
            plt.ylabel('Number of Samples')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'gene_count_distribution.png'), dpi=300)
            plt.close()
            
            # Calculate gene frequency across samples
            gene_freq = df.sum(axis=0) / df.shape[0]
            freq_bins = np.linspace(0, 1, 20)
            
            plt.figure(figsize=(12, 6))
            plt.hist(gene_freq, bins=freq_bins)
            plt.title('Distribution of Gene Frequencies Across Samples')
            plt.xlabel('Fraction of Samples Containing Gene')
            plt.ylabel('Number of Genes')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'gene_frequency_distribution.png'), dpi=300)
            plt.close()
    except Exception as e:
        print(f'Error creating summary plots for {{matrix_path}}: {{e}}')

# Create heatmaps
create_heatmap('{input.pangenome_matrix}', '{output.heatmap_dir}/pangenome_heatmap.png', 
               'Pan-genome Presence/Absence Heatmap')
create_heatmap('{input.ko_matrix}', '{output.heatmap_dir}/ko_heatmap.png', 
               'KO Abundance Heatmap')
create_heatmap('{input.phylum_matrix}', '{output.heatmap_dir}/phylum_heatmap.png', 
               'Phylum Abundance Heatmap')
create_heatmap('{input.species_matrix}', '{output.heatmap_dir}/species_heatmap.png', 
               'Species Abundance Heatmap')

# Create PCA plots
create_pca_plot('{input.pangenome_matrix}', '{input.integrated_metadata}',
                '{output.pca_dir}/pangenome_pca.png', 'Pan-genome PCA')
create_pca_plot('{input.ko_matrix}', '{input.integrated_metadata}',
                '{output.pca_dir}/ko_pca.png', 'KO Abundance PCA')
create_pca_plot('{input.phylum_matrix}', '{input.integrated_metadata}',
                '{output.pca_dir}/phylum_pca.png', 'Phylum Abundance PCA')
create_pca_plot('{input.species_matrix}', '{input.integrated_metadata}',
                '{output.pca_dir}/species_pca.png', 'Species Abundance PCA')

# Create summary plots
create_summary_plots('{input.pangenome_matrix}', '{output.summary_plots}', 'pangenome')
create_summary_plots('{input.ko_matrix}', '{output.summary_plots}', 'functional')
create_summary_plots('{input.phylum_matrix}', '{output.summary_plots}', 'taxonomy')
create_summary_plots('{input.species_matrix}', '{output.summary_plots}', 'taxonomy')
"
        """

# All cross-sample analyses rule
rule all_cross_sample:
    input:
        # Pan-genome matrices
        CROSS_SAMPLE_DIR + "/pan_genome_matrix.tsv",
        CROSS_SAMPLE_DIR + "/core_genes.txt",
        CROSS_SAMPLE_DIR + "/accessory_genes.txt",
        CROSS_SAMPLE_DIR + "/unique_genes.txt",
        
        # Functional profiles
        CROSS_SAMPLE_DIR + "/ko_abundance_matrix.tsv",
        CROSS_SAMPLE_DIR + "/pathway_abundance_matrix.tsv",
        
        # Taxonomic abundance
        CROSS_SAMPLE_DIR + "/phylum_abundance_matrix.tsv",
        CROSS_SAMPLE_DIR + "/genus_abundance_matrix.tsv",
        CROSS_SAMPLE_DIR + "/species_abundance_matrix.tsv",
        CROSS_SAMPLE_DIR + "/bin_taxonomy_matrix.tsv",
        
        # Integration and visualization
        CROSS_SAMPLE_DIR + "/integrated_metadata.tsv",
        CROSS_SAMPLE_DIR + "/matrix_statistics.tsv",
        CROSS_SAMPLE_DIR + "/visualizations/heatmaps",
        CROSS_SAMPLE_DIR + "/visualizations/pca",
        CROSS_SAMPLE_DIR + "/visualizations/summary"
    output:
        touch(CROSS_SAMPLE_DIR + "/cross_sample_analysis_complete.flag")




