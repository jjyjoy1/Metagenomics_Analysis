# Host-Microbiome Interactions Analysis Pipeline
# Identifies host-microbiome interactions, host gene content, and microbiome effects

# Configuration
configfile: "host_microbiome_config.yaml"

# Import sample information from manifest file
import pandas as pd
import os

# Read manifest file
manifest = pd.read_csv(config["manifest_file"])
SAMPLES = manifest["sample_id"].tolist()
SAMPLE_TYPES = manifest["sample_type"].tolist() if "sample_type" in manifest.columns else ["unknown"] * len(SAMPLES)

# Define the final output files
rule all:
    input:
        # Host read separation
        expand("results/host_microbiome/host_separation/{sample}.host.fastq.gz", sample=SAMPLES),
        expand("results/host_microbiome/host_separation/{sample}.microbiome.fastq.gz", sample=SAMPLES),
        
        # Host variant calling
        expand("results/host_microbiome/host_variants/{sample}.variants.vcf.gz", sample=SAMPLES),
        
        # Host gene expression (for RNA-seq data)
        "results/host_microbiome/host_expression/counts_matrix.tsv" if config["include_rnaseq"] else [],
        "results/host_microbiome/host_expression/deseq2_results.tsv" if config["include_rnaseq"] else [],
        
        # Microbial gene content
        expand("results/host_microbiome/microbial_genes/{sample}.genes.fna", sample=SAMPLES),
        "results/host_microbiome/microbial_genes/gene_catalog.fna",
        "results/host_microbiome/microbial_genes/gene_abundance_matrix.tsv",
        
        # Microbiome-host interactions
        "results/host_microbiome/interactions/putative_interactions.tsv",
        
        # Biomarker identification
        "results/host_microbiome/biomarkers/biomarker_species.tsv",
        "results/host_microbiome/biomarkers/biomarker_genes.tsv",
        "results/host_microbiome/biomarkers/biomarker_pathways.tsv",
        
        # Visualizations
        "results/host_microbiome/visualization/host_microbiome_network.png",
        "results/host_microbiome/visualization/biomarker_heatmap.png",
        
        # Final report
        "results/host_microbiome/reports/host_microbiome_report.html"

# Separate host reads from microbiome reads
rule separate_host_reads:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        host_r1="results/host_microbiome/host_separation/{sample}.host_R1.fastq.gz",
        host_r2="results/host_microbiome/host_separation/{sample}.host_R2.fastq.gz",
        microbiome_r1="results/host_microbiome/host_separation/{sample}.microbiome_R1.fastq.gz",
        microbiome_r2="results/host_microbiome/host_separation/{sample}.microbiome_R2.fastq.gz",
        host_single="results/host_microbiome/host_separation/{sample}.host.fastq.gz",
        microbiome_single="results/host_microbiome/host_separation/{sample}.microbiome.fastq.gz"
    params:
        host_genome=config["host_genome"],
        prefix="results/host_microbiome/host_separation/{sample}"
    threads: 16
    shell:
        """
        mkdir -p $(dirname {output.host_r1})
        
        # Align to host reference genome using Bowtie2
        bowtie2 -p {threads} -x {params.host_genome} \
            -1 {input.r1} -2 {input.r2} \
            --un-conc-gz {params.prefix}.microbiome_%.fastq.gz \
            --al-conc-gz {params.prefix}.host_%.fastq.gz \
            -S /dev/null
        
        # Combine reads for single-end processing
        cat {output.host_r1} {output.host_r2} > {output.host_single}
        cat {output.microbiome_r1} {output.microbiome_r2} > {output.microbiome_single}
        """

# Map host reads to reference genome for variant calling
rule map_host_reads:
    input:
        r1="results/host_microbiome/host_separation/{sample}.host_R1.fastq.gz",
        r2="results/host_microbiome/host_separation/{sample}.host_R2.fastq.gz"
    output:
        bam="results/host_microbiome/host_mapping/{sample}.sorted.bam",
        bai="results/host_microbiome/host_mapping/{sample}.sorted.bam.bai"
    params:
        host_genome=config["host_genome"]
    threads: 16
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        # Map to reference with BWA-MEM
        bwa mem -t {threads} {params.host_genome} {input.r1} {input.r2} | \
        samtools sort -@ {threads} -o {output.bam}
        
        # Index BAM file
        samtools index {output.bam}
        """

# Call variants in host genome
rule call_host_variants:
    input:
        bam="results/host_microbiome/host_mapping/{sample}.sorted.bam",
        bai="results/host_microbiome/host_mapping/{sample}.sorted.bam.bai"
    output:
        vcf="results/host_microbiome/host_variants/{sample}.variants.vcf.gz",
        tbi="results/host_microbiome/host_variants/{sample}.variants.vcf.gz.tbi"
    params:
        host_genome=config["host_genome"],
        host_genome_fai=config["host_genome"] + ".fai"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        # Call variants with FreeBayes
        freebayes -f {params.host_genome} \
                 --ploidy 2 \
                 {input.bam} | \
        bgzip > {output.vcf}
        
        # Index VCF
        tabix -p vcf {output.vcf}
        """

# RNA-seq processing (optional)
if config["include_rnaseq"]:
    rule quantify_host_expression:
        input:
            r1="results/host_microbiome/host_separation/{sample}.host_R1.fastq.gz",
            r2="results/host_microbiome/host_separation/{sample}.host_R2.fastq.gz"
        output:
            quant="results/host_microbiome/host_expression/{sample}/quant.sf"
        params:
            host_transcriptome=config["host_transcriptome"],
            outdir="results/host_microbiome/host_expression/{sample}"
        threads: 16
        shell:
            """
            mkdir -p {params.outdir}
            
            # Quantify with Salmon
            salmon quant -i {params.host_transcriptome} \
                         -l A \
                         -1 {input.r1} \
                         -2 {input.r2} \
                         -p {threads} \
                         --validateMappings \
                         -o {params.outdir}
            """
    
    rule merge_expression_counts:
        input:
            quants=expand("results/host_microbiome/host_expression/{sample}/quant.sf", sample=SAMPLES)
        output:
            matrix="results/host_microbiome/host_expression/counts_matrix.tsv"
        params:
            samples=SAMPLES,
            tx2gene=config["tx2gene_map"]
        script:
            "scripts/merge_expression.R"
    
    rule differential_expression:
        input:
            counts="results/host_microbiome/host_expression/counts_matrix.tsv"
        output:
            results="results/host_microbiome/host_expression/deseq2_results.tsv",
            plots="results/host_microbiome/host_expression/deseq2_plots.pdf"
        params:
            sample_info=config["sample_metadata"],
            design=config["deseq2_design"]
        script:
            "scripts/differential_expression.R"

# Predict genes from microbiome
rule predict_microbial_genes:
    input:
        contigs="results/assembly/{sample}/contigs.fasta"
    output:
        genes="results/host_microbiome/microbial_genes/{sample}.genes.fna",
        proteins="results/host_microbiome/microbial_genes/{sample}.proteins.faa"
    shell:
        """
        mkdir -p $(dirname {output.genes})
        
        # Predict genes with Prodigal
        prodigal -i {input.contigs} \
                -d {output.genes} \
                -a {output.proteins} \
                -p meta
        """

# Create non-redundant gene catalog
rule create_gene_catalog:
    input:
        genes=expand("results/host_microbiome/microbial_genes/{sample}.genes.fna", sample=SAMPLES)
    output:
        catalog="results/host_microbiome/microbial_genes/gene_catalog.fna",
        clusters="results/host_microbiome/microbial_genes/gene_clusters.tsv"
    params:
        identity=config["gene_catalog_identity"],
        coverage=config["gene_catalog_coverage"]
    threads: 16
    shell:
        """
        # Concatenate all genes
        cat {input.genes} > results/host_microbiome/microbial_genes/all_genes.fna
        
        # Cluster with CD-HIT
        cd-hit-est -i results/host_microbiome/microbial_genes/all_genes.fna \
                  -o {output.catalog} \
                  -c {params.identity} \
                  -aS {params.coverage} \
                  -T {threads} \
                  -M 0 \
                  -d 0 \
                  -g 1
        
        # Rename output clusters file
        mv {output.catalog}.clstr {output.clusters}
        """

# Quantify gene abundance
rule quantify_gene_abundance:
    input:
        catalog="results/host_microbiome/microbial_genes/gene_catalog.fna",
        r1="results/host_microbiome/host_separation/{sample}.microbiome_R1.fastq.gz",
        r2="results/host_microbiome/host_separation/{sample}.microbiome_R2.fastq.gz"
    output:
        bam="results/host_microbiome/microbial_genes/abundance/{sample}.bam",
        counts="results/host_microbiome/microbial_genes/abundance/{sample}.counts"
    threads: 16
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        # Build index for the gene catalog
        bwa index -a bwtsw {input.catalog}
        
        # Map reads to gene catalog
        bwa mem -t {threads} {input.catalog} {input.r1} {input.r2} | \
        samtools view -F 4 -b | \
        samtools sort -o {output.bam}
        
        # Count reads per gene
        samtools idxstats {output.bam} | cut -f 1,3 > {output.counts}
        """

# Merge gene abundance across samples
rule merge_gene_abundance:
    input:
        counts=expand("results/host_microbiome/microbial_genes/abundance/{sample}.counts", sample=SAMPLES)
    output:
        matrix="results/host_microbiome/microbial_genes/gene_abundance_matrix.tsv"
    params:
        samples=SAMPLES
    script:
        "scripts/merge_gene_abundance.py"

# Predict host-microbiome interactions
rule predict_interactions:
    input:
        gene_abundance="results/host_microbiome/microbial_genes/gene_abundance_matrix.tsv",
        host_variants=expand("results/host_microbiome/host_variants/{sample}.variants.vcf.gz", sample=SAMPLES),
        host_expression="results/host_microbiome/host_expression/counts_matrix.tsv" if config["include_rnaseq"] else []
    output:
        interactions="results/host_microbiome/interactions/putative_interactions.tsv"
    params:
        interaction_db=config["interaction_db"],
        sample_metadata=config["sample_metadata"]
    script:
        "scripts/predict_interactions.py"

# Identify biomarkers
rule identify_biomarkers:
    input:
        metaphlan=expand("results/metaphlan/{sample}.metaphlan.profile.txt", sample=SAMPLES),
        gene_abundance="results/host_microbiome/microbial_genes/gene_abundance_matrix.tsv",
        pathways="results/humann3/merged/pathabundance.tsv" if config["include_functional"] else []
    output:
        biomarker_species="results/host_microbiome/biomarkers/biomarker_species.tsv",
        biomarker_genes="results/host_microbiome/biomarkers/biomarker_genes.tsv",
        biomarker_pathways="results/host_microbiome/biomarkers/biomarker_pathways.tsv" if config["include_functional"] else []
    params:
        sample_metadata=config["sample_metadata"],
        group_column=config["group_column"]
    script:
        "scripts/identify_biomarkers.py"

# Create visualizations
rule visualize_host_microbiome:
    input:
        interactions="results/host_microbiome/interactions/putative_interactions.tsv",
        biomarker_species="results/host_microbiome/biomarkers/biomarker_species.tsv",
        biomarker_genes="results/host_microbiome/biomarkers/biomarker_genes.tsv"
    output:
        network="results/host_microbiome/visualization/host_microbiome_network.png",
        heatmap="results/host_microbiome/visualization/biomarker_heatmap.png",
        pca="results/host_microbiome/visualization/host_microbiome_pca.png"
    params:
        sample_metadata=config["sample_metadata"]
    script:
        "scripts/visualize_host_microbiome.py"

# Generate HTML report
rule host_microbiome_report:
    input:
        network="results/host_microbiome/visualization/host_microbiome_network.png",
        heatmap="results/host_microbiome/visualization/biomarker_heatmap.png",
        pca="results/host_microbiome/visualization/host_microbiome_pca.png",
        biomarker_species="results/host_microbiome/biomarkers/biomarker_species.tsv",
        interactions="results/host_microbiome/interactions/putative_interactions.tsv"
    output:
        report="results/host_microbiome/reports/host_microbiome_report.html"
    params:
        sample_metadata=config["sample_metadata"]
    script:
        "scripts/host_microbiome_report.py"
