# Snakemake workflow for metagenomics virus detection pipeline
# Filename: Snakefile

# Configuration
configfile: "config.yaml"

# Define wildcards for samples
SAMPLES = config["samples"]
HOST_GENOME = config["host_genome"]

# Define final output files that should be generated
rule all:
    input:
        # Preprocessing
        expand("results/preprocessing/{sample}_host_removed_1.fastq.gz", sample=SAMPLES),
        expand("results/preprocessing/{sample}_host_removed_2.fastq.gz", sample=SAMPLES),
        # Quality control
        expand("results/qc/{sample}_fastp_report.html", sample=SAMPLES),
        # Known virus detection with Kraken2
        expand("results/kraken2/{sample}_kraken2_report.txt", sample=SAMPLES),
        # Viral assembly
        expand("results/assembly/{sample}/contigs.fasta", sample=SAMPLES),
        # DIAMOND alignment against viral proteins
        expand("results/diamond/{sample}_viral_hits.txt", sample=SAMPLES),
        # VirSorter2 viral contig identification
        expand("results/virsorter2/{sample}/final-viral-score.tsv", sample=SAMPLES),
        # CheckV viral quality control
        expand("results/checkv/{sample}/quality_summary.tsv", sample=SAMPLES),
        # Viral annotation with Pharokka and VGAS
        expand("results/pharokka/{sample}/pharokka_output.gbk", sample=SAMPLES),
        expand("results/vgas/{sample}/annotation_results.txt", sample=SAMPLES),
        # Viral taxonomy with vConTACT2
        expand("results/vcontact2/{sample}/genome_by_genome_overview.csv", sample=SAMPLES),
        # Final report
        "results/final_report.html"

# Quality control with fastp
rule fastp:
    input:
        r1 = "raw_data/{sample}_R1.fastq.gz",
        r2 = "raw_data/{sample}_R2.fastq.gz"
    output:
        r1 = temp("results/preprocessing/{sample}_trimmed_1.fastq.gz"),
        r2 = temp("results/preprocessing/{sample}_trimmed_2.fastq.gz"),
        html = "results/qc/{sample}_fastp_report.html",
        json = "results/qc/{sample}_fastp_report.json"
    threads: 8
    resources:
        mem_mb = 8000
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
        --out1 {output.r1} --out2 {output.r2} \
        --html {output.html} --json {output.json} \
        --thread {threads} \
        --detect_adapter_for_pe \
        --cut_front --cut_tail --cut_mean_quality 20 \
        --length_required 50 \
        --correction \
        > {log} 2>&1
        """

# Host DNA removal with Bowtie2
rule remove_host:
    input:
        r1 = "results/preprocessing/{sample}_trimmed_1.fastq.gz",
        r2 = "results/preprocessing/{sample}_trimmed_2.fastq.gz",
        host_index = HOST_GENOME
    output:
        r1 = "results/preprocessing/{sample}_host_removed_1.fastq.gz",
        r2 = "results/preprocessing/{sample}_host_removed_2.fastq.gz",
        bam = temp("results/preprocessing/{sample}.host_aligned.bam")
    threads: 16
    resources:
        mem_mb = 16000
    log:
        "logs/remove_host/{sample}.log"
    shell:
        """
        # Align to host genome
        bowtie2 -p {threads} -x {input.host_index} \
        -1 {input.r1} -2 {input.r2} \
        --un-conc-gz results/preprocessing/{wildcards.sample}_host_removed_%.fastq.gz \
        --very-sensitive \
        | samtools view -bS - > {output.bam} \
        2> {log}
        """

# Known virus detection with Kraken2
rule kraken2:
    input:
        r1 = "results/preprocessing/{sample}_host_removed_1.fastq.gz",
        r2 = "results/preprocessing/{sample}_host_removed_2.fastq.gz",
        db = config["kraken2_db"]
    output:
        report = "results/kraken2/{sample}_kraken2_report.txt",
        output = "results/kraken2/{sample}_kraken2_output.txt"
    threads: 16
    resources:
        mem_mb = 64000
    log:
        "logs/kraken2/{sample}.log"
    shell:
        """
        kraken2 --db {input.db} \
        --threads {threads} \
        --output {output.output} \
        --report {output.report} \
        --paired {input.r1} {input.r2} \
        --use-names \
        2> {log}
        """

# Metagenomic assembly with metaSPAdes
rule metaspades:
    input:
        r1 = "results/preprocessing/{sample}_host_removed_1.fastq.gz",
        r2 = "results/preprocessing/{sample}_host_removed_2.fastq.gz"
    output:
        contigs = "results/assembly/{sample}/contigs.fasta",
        scaffolds = "results/assembly/{sample}/scaffolds.fasta"
    threads: 32
    resources:
        mem_mb = 128000,
        time = "24:00:00"
    log:
        "logs/metaspades/{sample}.log"
    shell:
        """
        # Create temporary unpacked files
        mkdir -p temp/{wildcards.sample}
        gzip -dc {input.r1} > temp/{wildcards.sample}/R1.fastq
        gzip -dc {input.r2} > temp/{wildcards.sample}/R2.fastq
        
        # Run metaSPAdes
        metaspades.py -1 temp/{wildcards.sample}/R1.fastq \
        -2 temp/{wildcards.sample}/R2.fastq \
        -o results/assembly/{wildcards.sample} \
        --meta \
        -t {threads} \
        -m {resources.mem_mb} \
        > {log} 2>&1
        
        # Clean up temporary files
        rm -rf temp/{wildcards.sample}
        """

# DIAMOND alignment against viral protein database
rule diamond:
    input:
        contigs = "results/assembly/{sample}/contigs.fasta",
        db = config["diamond_viral_db"]
    output:
        hits = "results/diamond/{sample}_viral_hits.txt",
        aln = "results/diamond/{sample}_alignments.daa"
    threads: 16
    resources:
        mem_mb = 32000
    log:
        "logs/diamond/{sample}.log"
    shell:
        """
        diamond blastx \
        --query {input.contigs} \
        --db {input.db} \
        --threads {threads} \
        --out {output.hits} \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle \
        --sensitive \
        --daa {output.aln} \
        --max-target-seqs 5 \
        --evalue 1e-5 \
        > {log} 2>&1
        """

# VirSorter2 for viral contig identification
rule virsorter2:
    input:
        contigs = "results/assembly/{sample}/contigs.fasta",
    output:
        score = "results/virsorter2/{sample}/final-viral-score.tsv",
        boundary = "results/virsorter2/{sample}/final-viral-boundary.tsv"
    threads: 16
    resources:
        mem_mb = 32000
    log:
        "logs/virsorter2/{sample}.log"
    shell:
        """
        virsorter run \
        --seqfile {input.contigs} \
        --min-length 1000 \
        --min-score 0.5 \
        -j {threads} \
        --include-groups "dsDNAphage,ssDNA" \
        --working-dir results/virsorter2/{wildcards.sample} \
        --keep-db-preparation \
        > {log} 2>&1
        """

# Extract viral contigs based on VirSorter2 results
rule extract_viruses:
    input:
        contigs = "results/assembly/{sample}/contigs.fasta",
        vs_results = "results/virsorter2/{sample}/final-viral-score.tsv"
    output:
        viral_contigs = "results/virsorter2/{sample}/viral_contigs.fasta"
    log:
        "logs/extract_viruses/{sample}.log"
    script:
        "scripts/extract_viral_contigs.py"

# CheckV for viral quality control
rule checkv:
    input:
        viral_contigs = "results/virsorter2/{sample}/viral_contigs.fasta",
        db = config["checkv_db"]
    output:
        quality = "results/checkv/{sample}/quality_summary.tsv",
        completeness = "results/checkv/{sample}/completeness.tsv",
        contamination = "results/checkv/{sample}/contamination.tsv"
    threads: 16
    resources:
        mem_mb = 32000
    log:
        "logs/checkv/{sample}.log"
    shell:
        """
        checkv end_to_end \
        {input.viral_contigs} \
        results/checkv/{wildcards.sample} \
        -t {threads} \
        -d {input.db} \
        > {log} 2>&1
        """

# Extract high-quality viral contigs
rule extract_hq_viruses:
    input:
        viral_contigs = "results/virsorter2/{sample}/viral_contigs.fasta",
        checkv_quality = "results/checkv/{sample}/quality_summary.tsv"
    output:
        hq_viruses = "results/checkv/{sample}/high_quality_viruses.fasta",
        mq_viruses = "results/checkv/{sample}/medium_quality_viruses.fasta"
    log:
        "logs/extract_hq_viruses/{sample}.log"
    script:
        "scripts/extract_quality_viruses.py"

# Pharokka for phage annotation
rule pharokka:
    input:
        hq_viruses = "results/checkv/{sample}/high_quality_viruses.fasta",
        db = config["pharokka_db"]
    output:
        gbk = "results/pharokka/{sample}/pharokka_output.gbk",
        gff = "results/pharokka/{sample}/pharokka_output.gff"
    threads: 16
    resources:
        mem_mb = 32000
    log:
        "logs/pharokka/{sample}.log"
    shell:
        """
        pharokka.py \
        --infile {input.hq_viruses} \
        --outdir results/pharokka/{wildcards.sample} \
        --threads {threads} \
        --database_dir {input.db} \
        --force \
        > {log} 2>&1
        """

# VGAS for general viral annotation
rule vgas:
    input:
        viruses = "results/checkv/{sample}/medium_quality_viruses.fasta"
    output:
        results = "results/vgas/{sample}/annotation_results.txt"
    threads: 16
    resources:
        mem_mb = 32000
    log:
        "logs/vgas/{sample}.log"
    shell:
        """
        VGAS \
        --infile {input.viruses} \
        --outdir results/vgas/{wildcards.sample} \
        --threads {threads} \
        --min_length 1000 \
        > {log} 2>&1
        """

# vConTACT2 for viral taxonomy and clustering
rule vcontact2:
    input:
        hq_viruses = "results/checkv/{sample}/high_quality_viruses.fasta",
        mq_viruses = "results/checkv/{sample}/medium_quality_viruses.fasta",
        gene_2_genome = "results/vcontact2/{sample}/gene_to_genome.csv",
        protein_file = "results/vcontact2/{sample}/viral_proteins.faa"
    output:
        overview = "results/vcontact2/{sample}/genome_by_genome_overview.csv",
        clusters = "results/vcontact2/{sample}/vcontact2_clusters.csv"
    threads: 16
    resources:
        mem_mb = 64000
    log:
        "logs/vcontact2/{sample}.log"
    shell:
        """
        vcontact2 \
        --raw-proteins {input.protein_file} \
        --proteins-fp {input.gene_2_genome} \
        --db 'ProkaryoticViralRefSeq94-Merged' \
        --output-dir results/vcontact2/{wildcards.sample} \
        --threads {threads} \
        --rel-mode 'Diamond' \
        --clusters-out {output.clusters} \
        --pcs-mode MCL \
        --vcs-mode ClusterONE \
        > {log} 2>&1
        """

# Generate proteins and gene-to-genome mapping for vConTACT2
rule prepare_vcontact2:
    input:
        hq_viruses = "results/checkv/{sample}/high_quality_viruses.fasta",
        mq_viruses = "results/checkv/{sample}/medium_quality_viruses.fasta"
    output:
        gene_2_genome = "results/vcontact2/{sample}/gene_to_genome.csv",
        protein_file = "results/vcontact2/{sample}/viral_proteins.faa"
    threads: 8
    resources:
        mem_mb = 16000
    log:
        "logs/prepare_vcontact2/{sample}.log"
    shell:
        """
        # Combine high and medium quality viruses
        cat {input.hq_viruses} {input.mq_viruses} > results/vcontact2/{wildcards.sample}/all_viruses.fasta
        
        # Run prodigal for gene prediction
        prodigal -i results/vcontact2/{wildcards.sample}/all_viruses.fasta \
        -a {output.protein_file} \
        -d results/vcontact2/{wildcards.sample}/viral_genes.fna \
        -f gff -p meta \
        > results/vcontact2/{wildcards.sample}/prodigal.gff 2> {log}
        
        # Generate gene-to-genome mapping file
        python scripts/create_gene2genome.py \
        --proteins {output.protein_file} \
        --output {output.gene_2_genome}
        """

# Generate final report
rule final_report:
    input:
        kraken_reports = expand("results/kraken2/{sample}_kraken2_report.txt", sample=SAMPLES),
        diamond_hits = expand("results/diamond/{sample}_viral_hits.txt", sample=SAMPLES),
        vs2_results = expand("results/virsorter2/{sample}/final-viral-score.tsv", sample=SAMPLES),
        checkv_results = expand("results/checkv/{sample}/quality_summary.tsv", sample=SAMPLES),
        vcontact2_results = expand("results/vcontact2/{sample}/genome_by_genome_overview.csv", sample=SAMPLES)
    output:
        report = "results/final_report.html"
    script:
        "scripts/generate_report.py"


