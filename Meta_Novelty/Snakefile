# Snakemake Pipeline for Novel Pathogen Detection
# Based on the provided Nextflow template

# Configuration
configfile: "config.yaml"

# Default parameters
SAMPLES = config.get("samples", [])
OUTPUT_DIR = config.get("output_dir", "results")
MIN_CONTIG_SIZE = config.get("min_contig_size", 1000)
THREADS = config.get("threads", 8)

# Define expected output files for the pipeline
rule all:
    input:
        # QC
        expand(OUTPUT_DIR + "/qc/{sample}_fastqc.html", sample=SAMPLES),
        # Assembly & Binning
        expand(OUTPUT_DIR + "/assembly/{sample}_contigs.fasta", sample=SAMPLES),
        expand(OUTPUT_DIR + "/graphbin/{sample}_graphbin_bins.csv", sample=SAMPLES),
        # Gene Prediction
        expand(OUTPUT_DIR + "/genes/{sample}_proteins.faa", sample=SAMPLES),
        # Taxonomic Classification
        expand(OUTPUT_DIR + "/cat/{sample}_cat.taxonomy.tsv", sample=SAMPLES),
        expand(OUTPUT_DIR + "/gtdbtk/{sample}_gtdbtk/classify/gtdbtk.bac120.summary.tsv", sample=SAMPLES),
        expand(OUTPUT_DIR + "/ani/{sample}_ani_results/ani_summary.tsv", sample=SAMPLES),
        # Viral Detection
        expand(OUTPUT_DIR + "/checkv/{sample}_high_confidence_viral_contigs.fasta", sample=SAMPLES),
        expand(OUTPUT_DIR + "/vcontact2/{sample}_vcontact2/genome_by_genome_overview.csv", sample=SAMPLES),
        # Embeddings
        expand(OUTPUT_DIR + "/embeddings/{sample}_protein_embeddings.npy", sample=SAMPLES),
        expand(OUTPUT_DIR + "/embeddings/{sample}_dna_embeddings.npy", sample=SAMPLES),
        # Integration & Analysis
        expand(OUTPUT_DIR + "/anomaly_detection/{sample}_anomalies_protein.tsv", sample=SAMPLES),
        expand(OUTPUT_DIR + "/reports/{sample}_comprehensive_report.html", sample=SAMPLES)

#############################################################################
# Part 1: Quality Control, Assembly & Binning
#############################################################################

# Quality control for raw reads
rule fastqc:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
    output:
        html = OUTPUT_DIR + "/qc/{sample}_fastqc.html",
        zip = OUTPUT_DIR + "/qc/{sample}_fastqc.zip"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTPUT_DIR}/qc
        fastqc -o {OUTPUT_DIR}/qc {input.r1} {input.r2}
        mv {OUTPUT_DIR}/qc/*_fastqc.html {output.html}
        mv {OUTPUT_DIR}/qc/*_fastqc.zip {output.zip}
        """

# Remove host sequences
rule host_removal:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
    output:
        r1_clean = OUTPUT_DIR + "/cleaned_reads/{sample}_clean_R1.fastq.gz",
        r2_clean = OUTPUT_DIR + "/cleaned_reads/{sample}_clean_R2.fastq.gz"
    params:
        ref_host = config.get("ref_host", "reference/hg38.fa")
    threads: THREADS
    shell:
        """
        mkdir -p {OUTPUT_DIR}/cleaned_reads
        bowtie2 -x {params.ref_host} -1 {input.r1} -2 {input.r2} \
            --un-conc-gz {OUTPUT_DIR}/cleaned_reads/{wildcards.sample}_clean_R%.fastq.gz \
            -S /dev/null -p {threads}
        """

# Assemble non-host reads
rule megahit_assembly:
    input:
        r1 = rules.host_removal.output.r1_clean,
        r2 = rules.host_removal.output.r2_clean
    output:
        contigs = OUTPUT_DIR + "/assembly/{sample}_contigs.fasta",
        graph = OUTPUT_DIR + "/assembly/{sample}_assembly_graph.gfa"
    params:
        min_contig_size = MIN_CONTIG_SIZE,
        out_dir = OUTPUT_DIR + "/assembly/{sample}_megahit"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTPUT_DIR}/assembly
        megahit -1 {input.r1} -2 {input.r2} -o {params.out_dir} --num-cpu-threads {threads}
        cp {params.out_dir}/final.contigs.fa {output.contigs}
        
        # Create assembly graph if it doesn't exist
        python -c "
import os
from Bio import SeqIO

# Create a minimal GFA if MEGAHIT doesn't provide one
if not os.path.exists('{params.out_dir}/assembly_graph.gfa'):
    with open('{output.graph}', 'w') as f:
        f.write('H\\tVN:Z:1.0\\n')
        for record in SeqIO.parse('{output.contigs}', 'fasta'):
            f.write(f'S\\t{{record.id}}\\t{{record.seq}}\\tLN:i:{{len(record.seq)}}\\n')
else:
    os.system('cp {params.out_dir}/assembly_graph.gfa {output.graph}')
"

        # Filter contigs by minimum size
        python -c "
from Bio import SeqIO
records = [r for r in SeqIO.parse('{output.contigs}', 'fasta') if len(r.seq) >= {params.min_contig_size}]
SeqIO.write(records, '{output.contigs}', 'fasta')
"
        """

# Map reads back to contigs for coverage information
rule read_mapping:
    input:
        contigs = rules.megahit_assembly.output.contigs,
        r1 = rules.host_removal.output.r1_clean,
        r2 = rules.host_removal.output.r2_clean
    output:
        bam = OUTPUT_DIR + "/mapping/{sample}.bam",
        bai = OUTPUT_DIR + "/mapping/{sample}.bam.bai",
        coverage = OUTPUT_DIR + "/mapping/{sample}_coverage.txt"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTPUT_DIR}/mapping
        
        # Index contigs
        bwa index {input.contigs}
        
        # Map reads to contigs
        bwa mem -t {threads} {input.contigs} {input.r1} {input.r2} | \
            samtools view -bS - > {output.bam}.temp
        
        # Sort and index BAM file
        samtools sort -@ {threads} {output.bam}.temp -o {output.bam}
        samtools index {output.bam}
        
        # Calculate coverage
        samtools depth -a {output.bam} | \
            awk '{{sum[$1]+=$3; cnt[$1]++}} END {{for (contig in sum) print contig, sum[contig]/cnt[contig]}}' > {output.coverage}
        
        rm -f {output.bam}.temp
        """

# Bin contigs into potential genomes
rule metabat2_binning:
    input:
        contigs = rules.megahit_assembly.output.contigs,
        coverage = rules.read_mapping.output.coverage
    output:
        bin_dir = directory(OUTPUT_DIR + "/bins/{sample}_bins"),
        summary = OUTPUT_DIR + "/bins/{sample}_bins.summary.tsv"
    params:
        min_contig_size = MIN_CONTIG_SIZE
    threads: THREADS
    shell:
        """
        mkdir -p {output.bin_dir}
        metabat2 -i {input.contigs} -a {input.coverage} -o {output.bin_dir}/bin \
            -m {params.min_contig_size} -t {threads} --unbinned
        
        # Create a summary of the bins
        echo -e "Bin_Id\\tContigs\\tTotal_Length\\tN50" > {output.summary}
        for bin in {output.bin_dir}/bin*.fa; do
            bin_id=$(basename ${{bin%.fa}})
            contig_count=$(grep -c ">" $bin)
            total_length=$(grep -v ">" $bin | tr -d '\\n' | wc -c)
            
            # Calculate N50
            python -c "
from Bio import SeqIO
lengths = sorted([len(rec.seq) for rec in SeqIO.parse('$bin', 'fasta')], reverse=True)
total = sum(lengths)
n50 = 0
running_sum = 0
for l in lengths:
    running_sum += l
    if running_sum >= total/2:
        n50 = l
        break
print(n50)
" > n50.txt
            n50=$(cat n50.txt)
            
            echo -e "$bin_id\\t$contig_count\\t$total_length\\t$n50" >> {output.summary}
        done
        """

# Improve binning using assembly graph
rule graphbin:
    input:
        contigs = rules.megahit_assembly.output.contigs,
        graph = rules.megahit_assembly.output.graph,
        bins = rules.metabat2_binning.output.bin_dir,
        bin_summary = rules.metabat2_binning.output.summary
    output:
        graphbin_dir = directory(OUTPUT_DIR + "/graphbin/{sample}_graphbin"),
        graphbin_csv = OUTPUT_DIR + "/graphbin/{sample}_graphbin_bins.csv"
    shell:
        """
        mkdir -p {output.graphbin_dir}
        
        # Prepare the initial binning result in the format required by GraphBin
        python -c "
import os
import glob

bins_path = '{input.bins}'
bin_files = glob.glob(os.path.join(bins_path, 'bin*.fa'))

with open('initial_contig_bins.csv', 'w') as outfile:
    outfile.write('contig_id,bin_id\\n')
    for bin_file in bin_files:
        bin_name = os.path.basename(bin_file).replace('.fa', '')
        with open(bin_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig_id = line.strip()[1:].split()[0]
                    outfile.write(f'{{contig_id}},{{bin_name}}\\n')
"
        
        # Run GraphBin
        graphbin --assembler megahit --graph {input.graph} --contigs {input.contigs} \
            --binned initial_contig_bins.csv --output {output.graphbin_csv}
        
        # Reorganize bins based on GraphBin output
        python -c "
import os
import shutil
from Bio import SeqIO

# Read GraphBin output
graphbin_result = {{}}
with open('{output.graphbin_csv}', 'r') as f:
    next(f)  # Skip header
    for line in f:
        contig_id, bin_id = line.strip().split(',')
        if bin_id != 'unbinned':
            graphbin_result[contig_id] = bin_id

# Create new bins
os.makedirs('{output.graphbin_dir}', exist_ok=True)
bin_contents = {{}}

for record in SeqIO.parse('{input.contigs}', 'fasta'):
    contig_id = record.id
    if contig_id in graphbin_result:
        bin_id = graphbin_result[contig_id]
        if bin_id not in bin_contents:
            bin_contents[bin_id] = []
        bin_contents[bin_id].append(record)

# Write new bin files
for bin_id, records in bin_contents.items():
    SeqIO.write(records, os.path.join('{output.graphbin_dir}', f'{{bin_id}}.fa'), 'fasta')
"
        """

#############################################################################
# Part 2: Gene Prediction & Taxonomic Classification
#############################################################################

# Predict genes from contigs
rule gene_prediction:
    input:
        contigs = rules.megahit_assembly.output.contigs
    output:
        gff = OUTPUT_DIR + "/genes/{sample}_genes.gff",
        proteins = OUTPUT_DIR + "/genes/{sample}_proteins.faa",
        genes = OUTPUT_DIR + "/genes/{sample}_genes.fna"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/genes
        prodigal -i {input.contigs} -o {output.gff} -a {output.proteins} -d {output.genes} -p meta -q
        """

# CAT taxonomic classification
rule cat_classification:
    input:
        contigs = rules.megahit_assembly.output.contigs,
        proteins = rules.gene_prediction.output.proteins
    output:
        cat_dir = directory(OUTPUT_DIR + "/cat/{sample}_cat"),
        taxonomy = OUTPUT_DIR + "/cat/{sample}_cat.taxonomy.tsv"
    threads: THREADS
    shell:
        """
        mkdir -p {output.cat_dir}
        
        # Run CAT on contigs
        CAT contigs -c {input.contigs} -p {input.proteins} -o {output.cat_dir}/{wildcards.sample}_cat \
            -t {threads}
        
        # Generate taxonomy file
        CAT add_names -i {output.cat_dir}/{wildcards.sample}_cat.contig2classification.txt \
            -o {output.taxonomy}
        """

# BAT for bin classification
rule bat_classification:
    input:
        graphbin_dir = rules.graphbin.output.graphbin_dir,
        graphbin_csv = rules.graphbin.output.graphbin_csv,
        proteins = rules.gene_prediction.output.proteins
    output:
        bat_dir = directory(OUTPUT_DIR + "/bat/{sample}_bat"),
        taxonomy = OUTPUT_DIR + "/bat/{sample}_bat.taxonomy.tsv"
    threads: THREADS
    shell:
        """
        mkdir -p {output.bat_dir}
        
        # Run BAT on bins
        BAT bins -i {input.graphbin_dir} -p {input.proteins} -o {output.bat_dir}/{wildcards.sample}_bat \
            -t {threads}
        
        # Generate taxonomy file
        CAT add_names -i {output.bat_dir}/{wildcards.sample}_bat.bin2classification.txt \
            -o {output.taxonomy}
        """

# MMseqs2 clustering and classification
rule mmseqs2_classification:
    input:
        contigs = rules.megahit_assembly.output.contigs,
        proteins = rules.gene_prediction.output.proteins
    output:
        mmseqs2_dir = directory(OUTPUT_DIR + "/mmseqs2/{sample}_mmseqs2"),
        taxonomy = OUTPUT_DIR + "/mmseqs2/{sample}_mmseqs2_taxonomy.tsv"
    params:
        mmseqs2_db = config.get("mmseqs2_db", "/path/to/mmseqs2_db/mmseqs2_db")
    threads: THREADS
    shell:
        """
        mkdir -p {output.mmseqs2_dir}
        
        # Create MMseqs2 database
        mmseqs createdb {input.proteins} {output.mmseqs2_dir}/query_db
        
        # Run MMseqs2 taxonomy assignment
        mmseqs taxonomy {output.mmseqs2_dir}/query_db {params.mmseqs2_db} \
            {output.mmseqs2_dir}/taxonomy_results tmp \
            --threads {threads} --tax-lineage true
        
        # Get taxonomy report
        mmseqs createtsv {output.mmseqs2_dir}/query_db {output.mmseqs2_dir}/taxonomy_results \
            {output.taxonomy} --threads {threads}
        """

# GTDB-Tk for bacterial phylogenetic placement
rule gtdbtk:
    input:
        graphbin_dir = rules.graphbin.output.graphbin_dir
    output:
        gtdbtk_dir = directory(OUTPUT_DIR + "/gtdbtk/{sample}_gtdbtk"),
        bac_summary = OUTPUT_DIR + "/gtdbtk/{sample}_gtdbtk/classify/gtdbtk.bac120.summary.tsv",
        ar_summary = OUTPUT_DIR + "/gtdbtk/{sample}_gtdbtk/classify/gtdbtk.ar53.summary.tsv"
    params:
        gtdb_db = config.get("gtdb_db", "/path/to/gtdb_db")
    threads: THREADS
    shell:
        """
        mkdir -p {output.gtdbtk_dir}
        export GTDBTK_DATA_PATH={params.gtdb_db}
        
        # Run GTDB-Tk classify workflow
        gtdbtk classify_wf --genome_dir {input.graphbin_dir} \
            --out_dir {output.gtdbtk_dir} \
            --cpus {threads} \
            --extension fa
            
        # Create empty summary files if they don't exist
        if [ ! -f {output.bac_summary} ]; then
            mkdir -p $(dirname {output.bac_summary})
            echo -e "user_genome\\tclassification\\tfastani_reference\\tfastani_reference_radius\\tfastani_taxonomy\\tfastani_ani\\tfastani_af\\tother_related_references(genome_id,species_name,radius,ANI,AF)\\tpplacer_taxonomy\\ttyping_status\\tother_related_classifications\\tfastani_reference_is_rep_radius\\tpplacer_placement_method\\taa_percent\\tred_value\\twarnings\\tnote" > {output.bac_summary}
        fi
        
        if [ ! -f {output.ar_summary} ]; then
            mkdir -p $(dirname {output.ar_summary})
            echo -e "user_genome\\tclassification\\tfastani_reference\\tfastani_reference_radius\\tfastani_taxonomy\\tfastani_ani\\tfastani_af\\tother_related_references(genome_id,species_name,radius,ANI,AF)\\tpplacer_taxonomy\\ttyping_status\\tother_related_classifications\\tfastani_reference_is_rep_radius\\tpplacer_placement_method\\taa_percent\\tred_value\\twarnings\\tnote" > {output.ar_summary}
        fi
        """

# Calculate Average Nucleotide Identity (ANI)
rule ani:
    input:
        graphbin_dir = rules.graphbin.output.graphbin_dir,
        gtdbtk_dir = rules.gtdbtk.output.gtdbtk_dir
    output:
        ani_dir = directory(OUTPUT_DIR + "/ani/{sample}_ani_results"),
        summary = OUTPUT_DIR + "/ani/{sample}_ani_results/ani_summary.tsv"
    params:
        gtdb_genomes = config.get("gtdb_genomes", "/path/to/gtdb_genomes")
    threads: THREADS
    shell:
        """
        mkdir -p {output.ani_dir}
        
        # Extract reference genomes based on GTDB-Tk results
        python -c "
import os
import pandas as pd
import glob

# Check if classification files exist
bac_file = '{input.gtdbtk_dir}/classify/gtdbtk.bac120.summary.tsv'
ar_file = '{input.gtdbtk_dir}/classify/gtdbtk.ar53.summary.tsv'

ref_genomes = {{}}

if os.path.exists(bac_file) and os.path.getsize(bac_file) > 0:
    try:
        bac_df = pd.read_csv(bac_file, sep='\\t')
        for idx, row in bac_df.iterrows():
            bin_id = row['user_genome']
            if 'fastani_reference' in row and pd.notna(row['fastani_reference']):
                closest_ref = row['fastani_reference']
                ref_genomes[bin_id] = closest_ref
    except Exception as e:
        print(f'Error processing bacterial file: {{e}}')

if os.path.exists(ar_file) and os.path.getsize(ar_file) > 0:
    try:
        ar_df = pd.read_csv(ar_file, sep='\\t')
        for idx, row in ar_df.iterrows():
            bin_id = row['user_genome']
            if 'fastani_reference' in row and pd.notna(row['fastani_reference']):
                closest_ref = row['fastani_reference']
                ref_genomes[bin_id] = closest_ref
    except Exception as e:
        print(f'Error processing archaeal file: {{e}}')

with open('reference_genomes.txt', 'w') as f:
    for bin_id, ref in ref_genomes.items():
        f.write(f'{{bin_id}}\\t{{ref}}\\n')
"
        
        # Run FastANI for each bin against its closest reference
        echo -e "Bin_ID\\tReference_Genome\\tANI\\tQuery_Coverage\\tReference_Coverage" > {output.summary}
        
        if [ -s reference_genomes.txt ]; then
            while read bin_id ref_genome; do
                ref_path="{params.gtdb_genomes}/$ref_genome.fna"
                
                # Create a placeholder if reference doesn't exist
                if [ ! -f "$ref_path" ]; then
                    echo -e "$bin_id\\t$ref_genome\\tNA\\tNA\\tNA" >> {output.summary}
                    continue
                fi
                
                # Run FastANI
                fastANI -q {input.graphbin_dir}/$bin_id.fa -r $ref_path \
                    -o {output.ani_dir}/$bin_id"_vs_"$ref_genome.ani
                
                # Process results
                if [ -s {output.ani_dir}/$bin_id"_vs_"$ref_genome.ani ]; then
                    ani_value=$(cut -f3 {output.ani_dir}/$bin_id"_vs_"$ref_genome.ani)
                    query_cov=$(cut -f4 {output.ani_dir}/$bin_id"_vs_"$ref_genome.ani)
                    ref_cov=$(cut -f5 {output.ani_dir}/$bin_id"_vs_"$ref_genome.ani)
                    echo -e "$bin_id\\t$ref_genome\\t$ani_value\\t$query_cov\\t$ref_cov" >> {output.summary}
                else
                    echo -e "$bin_id\\t$ref_genome\\tNA\\tNA\\tNA" >> {output.summary}
                fi
            done < reference_genomes.txt
        else
            echo "No reference genomes found in GTDB-Tk results" >> {output.ani_dir}/README.txt
        fi
        """

#############################################################################
# Part 3: Viral Detection
#############################################################################

# Viral detection using VirSorter2
rule virsorter2:
    input:
        contigs = rules.megahit_assembly.output.contigs
    output:
        vs2_dir = directory(OUTPUT_DIR + "/virsorter2/{sample}_virsorter2"),
        viral_contigs = OUTPUT_DIR + "/virsorter2/{sample}_viral_contigs_vs2.fasta"
    params:
        min_contig_size = MIN_CONTIG_SIZE
    threads: THREADS
    shell:
        """
        mkdir -p {output.vs2_dir}
        
        # Run VirSorter2
        virsorter run --seqfile {input.contigs} --working-dir {output.vs2_dir} \
            --include-groups dsDNAphage,ssDNA,RNA --min-length {params.min_contig_size} \
            -j {threads} --keep-original-seq
        
        # Extract viral sequences
        python -c "
import os
import pandas as pd
from Bio import SeqIO

# Read the VirSorter2 results
results_file = '{output.vs2_dir}/final-viral-score.tsv'
if os.path.exists(results_file):
    df = pd.read_csv(results_file, sep='\\t')
    viral_contigs = df[df['max_score'] >= 0.9]['seqname'].tolist()
else:
    viral_contigs = []

# Extract the viral sequences
with open('{input.contigs}', 'r') as input_handle:
    seqs = SeqIO.to_dict(SeqIO.parse(input_handle, 'fasta'))
    
viral_seqs = [seqs[contig] for contig in viral_contigs if contig in seqs]

with open('{output.viral_contigs}', 'w') as output_handle:
    SeqIO.write(viral_seqs, output_handle, 'fasta')
"
        
        # Create empty file if no viral contigs are found
        if [ ! -s {output.viral_contigs} ]; then
            touch {output.viral_contigs}
        fi
        """

# Viral detection using VIBRANT
rule vibrant:
    input:
        contigs = rules.megahit_assembly.output.contigs
    output:
        vibrant_dir = directory(OUTPUT_DIR + "/vibrant/{sample}_vibrant"),
        viral_contigs = OUTPUT_DIR + "/vibrant/{sample}_viral_contigs_vibrant.fasta"
    params:
        min_contig_size = MIN_CONTIG_SIZE
    threads: THREADS
    shell:
        """
        mkdir -p {output.vibrant_dir}
        
        # Run VIBRANT
        VIBRANT_run.py -i {input.contigs} -folder {output.vibrant_dir} -t {threads} \
            -l {params.min_contig_size}
        
        # Extract viral sequences
        if [ -f "{output.vibrant_dir}/VIBRANT_phages_{wildcards.sample}_contigs.fasta/VIBRANT_phages_{wildcards.sample}_contigs.fasta.fasta" ]; then
            cp "{output.vibrant_dir}/VIBRANT_phages_{wildcards.sample}_contigs.fasta/VIBRANT_phages_{wildcards.sample}_contigs.fasta.fasta" "{output.viral_contigs}"
        else
            touch "{output.viral_contigs}"
        fi
        """

# Additional viral detection with DeepVirFinder
rule deep_vir_finder:
    input:
        contigs = rules.megahit_assembly.output.contigs
    output:
        dvf_dir = directory(OUTPUT_DIR + "/deepvirfinder/{sample}_dvf"),
        viral_contigs = OUTPUT_DIR + "/deepvirfinder/{sample}_viral_contigs_dvf.fasta"
    params:
        min_contig_size = MIN_CONTIG_SIZE,
        dvf_script = config.get("dvf_script", "/path/to/DeepVirFinder/dvf.py")
    threads: THREADS
    shell:
        """
        mkdir -p {output.dvf_dir}
        
        # Run DeepVirFinder
        python {params.dvf_script} -i {input.contigs} -o {output.dvf_dir} -c {threads} -l {params.min_contig_size}
        
        # Extract viral sequences with score > 0.9 and p-value < 0.05
        python -c "
import pandas as pd
from Bio import SeqIO

# Read DeepVirFinder results
result_file = '{output.dvf_dir}/{wildcards.sample}_contigs.fasta.dvfpred.txt'
try:
    df = pd.read_csv(result_file, sep='\\t')
    viral_contigs = df[(df['score'] > 0.9) & (df['pvalue'] < 0.05)]['name'].tolist()
except Exception as e:
    print(f'Error reading DVF results: {{e}}')
    viral_contigs = []

# Extract viral sequences
with open('{input.contigs}', 'r') as input_handle:
    seqs = SeqIO.to_dict(SeqIO.parse(input_handle, 'fasta'))
    
viral_seqs = [seqs[contig] for contig in viral_contigs if contig in seqs]

with open('{output.viral_contigs}', 'w') as output_handle:
    SeqIO.write(viral_seqs, output_handle, 'fasta')
"
        
        # Create empty file if no viral contigs are found
        if [ ! -s {output.viral_contigs} ]; then
            touch {output.viral_contigs}
        fi
        """

# Combine viral predictions and assess quality with CheckV
rule checkv:
    input:
        vs2_dir = rules.virsorter2.output.vs2_dir,
        vs2_fasta = rules.virsorter2.output.viral_contigs,
        vibrant_dir = rules.vibrant.output.vibrant_dir,
        vibrant_fasta = rules.vibrant.output.viral_contigs,
        dvf_dir = rules.deep_vir_finder.output.dvf_dir,
        dvf_fasta = rules.deep_vir_finder.output.viral_contigs
    output:
        checkv_dir = directory(OUTPUT_DIR + "/checkv/{sample}_checkv"),
        combined_viral = OUTPUT_DIR + "/checkv/{sample}_combined_viral_contigs.fasta",
        high_conf_viral = OUTPUT_DIR + "/checkv/{sample}_high_confidence_viral_contigs.fasta"
    params:
        checkv_db = config.get("checkv_db", "/path/to/checkv_db")
    threads: THREADS
    shell:
        """
        mkdir -p {output.checkv_dir}
        
        # Combine viral predictions
        cat {input.vs2_fasta} {input.vibrant_fasta} {input.dvf_fasta} > {output.checkv_dir}/all_viral_predictions.fasta
        
        # Remove duplicates
        python -c "
from Bio import SeqIO

# Read all sequences
seqs = {{}}
for record in SeqIO.parse('{output.checkv_dir}/all_viral_predictions.fasta', 'fasta'):
    if record.id not in seqs:
        seqs[record.id] = record

# Write unique sequences
with open('{output.combined_viral}', 'w') as outf:
    SeqIO.write(seqs.values(), outf, 'fasta')
"
        
        # Check if there are any viral contigs
        if [ -s {output.combined_viral} ]; then
            # Run CheckV
            checkv end_to_end {output.combined_viral} {output.checkv_dir} \
                -t {threads} -d {params.checkv_db}
                
            # Extract high-confidence viral contigs
            python -c "
import pandas as pd
from Bio import SeqIO

# Read CheckV quality results
quality_file = '{output.checkv_dir}/quality_summary.tsv'
try:
    df = pd.read_csv(quality_file, sep='\t')
    high_confidence = df[(df['checkv_quality'] == 'High-quality') | 
                        (df['checkv_quality'] == 'Complete')]['contig_id'].tolist()
except Exception as e:
    print(f'Error reading CheckV results: {{e}}')
    high_confidence = []

# Extract high-confidence sequences
with open('{output.combined_viral}', 'r') as input_handle:
    seqs = SeqIO.to_dict(SeqIO.parse(input_handle, 'fasta'))
    
high_conf_seqs = [seqs[contig] for contig in high_confidence if contig in seqs]

with open('{output.high_conf_viral}', 'w') as output_handle:
    SeqIO.write(high_conf_seqs, output_handle, 'fasta')
"
        else
            # Create empty files if no viral contigs
            touch {output.high_conf_viral}
        fi
        """

# vConTACT2 for viral classification and networks
rule vcontact2:
    input:
        checkv_dir = rules.checkv.output.checkv_dir,
        combined_viral = rules.checkv.output.combined_viral,
        high_conf_viral = rules.checkv.output.high_conf_viral,
        proteins = rules.gene_prediction.output.proteins
    output:
        vcontact2_dir = directory(OUTPUT_DIR + "/vcontact2/{sample}_vcontact2"),
        overview = OUTPUT_DIR + "/vcontact2/{sample}_vcontact2/genome_by_genome_overview.csv"
    params:
        vcontact2_db = config.get("vcontact2_db", "/path/to/vcontact2_db")
    threads: THREADS
    shell:
        """
        mkdir -p {output.vcontact2_dir}
        
        # Check if there are viral contigs
        if [ -s {input.combined_viral} ]; then
            # Extract viral proteins
            python -c "
from Bio import SeqIO

# Get viral contig IDs
viral_ids = set()
for record in SeqIO.parse('{input.combined_viral}', 'fasta'):
    viral_ids.add(record.id)

# Extract viral proteins
viral_proteins = []
for record in SeqIO.parse('{input.proteins}', 'fasta'):
    contig_id = record.id.split('_')[0]
    if contig_id in viral_ids:
        viral_proteins.append(record)

# Write viral proteins
with open('viral_proteins.faa', 'w') as outf:
    SeqIO.write(viral_proteins, outf, 'fasta')
"
            
            # Run vConTACT2
            vcontact2 --raw-proteins viral_proteins.faa --rel-mode 'Diamond' \
                --proteins-fp viral_proteins.faa --db {params.vcontact2_db} \
                --output-dir {output.vcontact2_dir} --threads {threads}
        else
            # Create empty directory structure
            mkdir -p {output.vcontact2_dir}
            touch {output.overview}
        fi
        """

#############################################################################
# Part 4: Embedding Models and Feature Extraction
#############################################################################

# Generate protein embeddings for further analysis
rule protein_embeddings:
    input:
        proteins = rules.gene_prediction.output.proteins
    output:
        embeddings = OUTPUT_DIR + "/embeddings/{sample}_protein_embeddings.npy",
        ids = OUTPUT_DIR + "/embeddings/{sample}_protein_ids.txt"
    params:
        model = config.get("protein_embedding_model", "esm2")
    threads: THREADS
    shell:
        """
        mkdir -p {OUTPUT_DIR}/embeddings
        
        # Run protein embedding generation
        python run_protein_embedding.py \
            --input {input.proteins} \
            --model {params.model} \
            --output {output.embeddings} \
            --ids {output.ids} \
            --threads {threads}
        """

# Generate DNA sequence embeddings
rule dna_embeddings:
    input:
        contigs = rules.megahit_assembly.output.contigs
    output:
        embeddings = OUTPUT_DIR + "/embeddings/{sample}_dna_embeddings.npy",
        ids = OUTPUT_DIR + "/embeddings/{sample}_contig_ids.txt"
    params:
        model = config.get("dna_embedding_model", "dnabert"),
        kmer_size = 6
    threads: THREADS
    shell:
        """
        mkdir -p {OUTPUT_DIR}/embeddings
        
        # Extract sequences in the right format for embedding model
        python -c "
from Bio import SeqIO
import numpy as np

seqs = []
ids = []
for record in SeqIO.parse('{input.contigs}', 'fasta'):
    seqs.append(str(record.seq))
    ids.append(record.id)

# Save IDs for later reference
with open('{output.ids}', 'w') as f:
    f.write('\\n'.join(ids))

# Save sequences for embedding processing
with open('sequences.txt', 'w') as f:
    f.write('\\n'.join(seqs))
"
        
        # Run DNA embedding generation
        python run_dna_embedding.py \
            --input sequences.txt \
            --model {params.model} \
            --output {output.embeddings} \
            --threads {threads} \
            --kmer_size {params.kmer_size}
        """

#############################################################################
# Part 5: Integration and Analysis
#############################################################################

# Anomaly detection using embeddings
rule anomaly_detection:
    input:
        protein_embeddings = rules.protein_embeddings.output.embeddings,
        protein_ids = rules.protein_embeddings.output.ids,
        dna_embeddings = rules.dna_embeddings.output.embeddings,
        contig_ids = rules.dna_embeddings.output.ids,
        cat_dir = rules.cat_classification.output.cat_dir,
        cat_taxonomy = rules.cat_classification.output.taxonomy,
        mmseqs2_dir = rules.mmseqs2_classification.output.mmseqs2_dir,
        mmseqs2_taxonomy = rules.mmseqs2_classification.output.taxonomy,
        checkv_dir = rules.checkv.output.checkv_dir,
        combined_viral = rules.checkv.output.combined_viral,
        high_conf_viral = rules.checkv.output.high_conf_viral
    output:
        anomalies_protein = OUTPUT_DIR + "/anomaly_detection/{sample}_anomalies_protein.tsv",
        anomalies_dna = OUTPUT_DIR + "/anomaly_detection/{sample}_anomalies_dna.tsv",
        anomaly_models = directory(OUTPUT_DIR + "/anomaly_detection/{sample}_anomaly_models"),
        visualization = directory(OUTPUT_DIR + "/anomaly_detection/{sample}_visualization")
    params:
        methods = config.get("anomaly_detection_methods", "isolation_forest,vae,dbscan"),
        contamination = 0.05,
        dim_reduction = "umap,tsne",
        perplexity = 30
    threads: THREADS
    shell:
        """
        mkdir -p {output.anomaly_models}
        mkdir -p {output.visualization}
        
        # Prepare annotation files for contextual anomaly detection
        python -c "
import pandas as pd
import numpy as np

# Load CAT taxonomy
try:
    cat_df = pd.read_csv('{input.cat_taxonomy}', sep='\t')
    cat_df = cat_df.rename(columns={{'# contig': 'contig_id'}})
    cat_dict = dict(zip(cat_df['contig_id'], cat_df['classification']))
except Exception as e:
    print(f'Error loading CAT taxonomy: {{e}}')
    cat_dict = {{}}

# Load MMseqs2 taxonomy
try:
    mmseqs_df = pd.read_csv('{input.mmseqs2_taxonomy}', sep='\t', header=None)
    mmseqs_df.columns = ['query', 'taxonomy', 'rank', 'taxid', 'score']
    mmseqs_dict = dict(zip(mmseqs_df['query'], mmseqs_df['taxonomy']))
except Exception as e:
    print(f'Error loading MMseqs2 taxonomy: {{e}}')
    mmseqs_dict = {{}}

# Load viral contigs
viral_contigs = set()
try:
    with open('{input.combined_viral}') as f:
        for line in f:
            if line.startswith('>'):
                contig = line.strip().lstrip('>')
                viral_contigs.add(contig)
except Exception as e:
    print(f'Error loading viral contigs: {{e}}')

# Load protein IDs and contig IDs
protein_ids = []
with open('{input.protein_ids}') as f:
    protein_ids = [line.strip() for line in f]

contig_ids = []
with open('{input.contig_ids}') as f:
    contig_ids = [line.strip() for line in f]

# Create annotation files
protein_anno = {{}}
for pid in protein_ids:
    contig = pid.split('_')[0]  # Extract contig ID from protein ID
    
    if contig in viral_contigs:
        protein_anno[pid] = 'viral'
    elif contig in cat_dict:
        protein_anno[pid] = cat_dict[contig]
    elif pid in mmseqs_dict:
        protein_anno[pid] = mmseqs_dict[pid]
    else:
        protein_anno[pid] = 'unknown'

contig_anno = {{}}
for cid in contig_ids:
    if cid in viral_contigs:
        contig_anno[cid] = 'viral'
    elif cid in cat_dict:
        contig_anno[cid] = cat_dict[cid]
    else:
        contig_anno[cid] = 'unknown'

# Write annotation files
with open('protein_annotations.tsv', 'w') as f:
    f.write('protein_id\\ttaxonomy\\n')
    for pid, anno in protein_anno.items():
        f.write(f'{{pid}}\\t{{anno}}\\n')

with open('contig_annotations.tsv', 'w') as f:
    f.write('contig_id\\ttaxonomy\\n')
    for cid, anno in contig_anno.items():
        f.write(f'{{cid}}\\t{{anno}}\\n')
"
        
        # Run anomaly detection using both protein and DNA embeddings
        python run_anomaly_detection.py \
            --protein_embeddings {input.protein_embeddings} \
            --protein_ids {input.protein_ids} \
            --protein_annotations protein_annotations.tsv \
            --dna_embeddings {input.dna_embeddings} \
            --contig_ids {input.contig_ids} \
            --contig_annotations contig_annotations.tsv \
            --methods "{params.methods}" \
            --output_protein {output.anomalies_protein} \
            --output_dna {output.anomalies_dna} \
            --save_models {output.anomaly_models} \
            --visualize {output.visualization} \
            --threads {threads} \
            --contamination {params.contamination} \
            --dimensions_reduction {params.dim_reduction} \
            --perplexity {params.perplexity}
        """

# Generate comprehensive report integrating all results
rule comprehensive_report:
    input:
        anomalies_protein = rules.anomaly_detection.output.anomalies_protein,
        anomalies_dna = rules.anomaly_detection.output.anomalies_dna,
        anomaly_models = rules.anomaly_detection.output.anomaly_models,
        visualization = rules.anomaly_detection.output.visualization,
        cat_dir = rules.cat_classification.output.cat_dir,
        cat_taxonomy = rules.cat_classification.output.taxonomy,
        gtdbtk_dir = rules.gtdbtk.output.gtdbtk_dir,
        ani_results = rules.ani.output.ani_dir,
        checkv_dir = rules.checkv.output.checkv_dir,
        combined_viral = rules.checkv.output.combined_viral,
        high_conf_viral = rules.checkv.output.high_conf_viral,
        vcontact2_dir = rules.vcontact2.output.vcontact2_dir
    output:
        report = OUTPUT_DIR + "/reports/{sample}_comprehensive_report.html",
        novel_pathogens = OUTPUT_DIR + "/reports/{sample}_novel_pathogens.tsv"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/reports
        
        # Generate comprehensive report
        python generate_comprehensive_report.py \
            --sample_id {wildcards.sample} \
            --anomalies_protein {input.anomalies_protein} \
            --anomalies_dna {input.anomalies_dna} \
            --visualization_dir {input.visualization} \
            --cat_taxonomy {input.cat_taxonomy} \
            --gtdbtk_dir {input.gtdbtk_dir} \
            --ani_dir {input.ani_results} \
            --checkv_dir {input.checkv_dir} \
            --vcontact2_dir {input.vcontact2_dir} \
            --high_confidence_viral {input.high_conf_viral} \
            --output_report {output.report} \
            --output_novel_pathogens {output.novel_pathogens}
        """

