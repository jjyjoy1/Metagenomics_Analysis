# Mobile Genetic Elements Analysis Pipeline
# Identifies plasmids, phages, transposons, and other mobile genetic elements

# Configuration
configfile: "mge_config.yaml"

# Import sample information from manifest file
import pandas as pd
import os

# Read manifest file
manifest = pd.read_csv(config["manifest_file"])
SAMPLES = manifest["sample_id"].tolist()

# Define the final output files
rule all:
    input:
        # Plasmid detection
        expand("results/mge/plasmids/plasmidfinder/{sample}.plasmids.tsv", sample=SAMPLES),
        expand("results/mge/plasmids/platon/{sample}/platon_summary.tsv", sample=SAMPLES),
        
        # Phage detection
        expand("results/mge/phages/vibrant/{sample}/VIBRANT_phages_{sample}/VIBRANT_results_{sample}/VIBRANT_phages_{sample}.tsv", sample=SAMPLES),
        expand("results/mge/phages/phigaro/{sample}/phigaro.tsv", sample=SAMPLES),
        
        # Transposon detection
        expand("results/mge/transposons/{sample}.transposons.tsv", sample=SAMPLES),
        
        # Insertion sequences
        expand("results/mge/insertion_sequences/{sample}.is.tsv", sample=SAMPLES),
        
        # Integrated MGE detection
        expand("results/mge/integrated/{sample}.integrated_mge.tsv", sample=SAMPLES),
        
        # Combined MGE report
        "results/mge/reports/combined_mge_report.tsv",
        
        # Visualizations
        "results/mge/visualization/mge_heatmap.png",
        "results/mge/visualization/mge_network.png",
        
        # Final MGE report
        "results/mge/reports/mge_analysis_report.html"

# Use assemblies from the main pipeline
rule get_assemblies:
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/mge/assemblies/{sample}.contigs.fasta"
    shell:
        """
        mkdir -p results/mge/assemblies
        cp {input} {output}
        """

# Plasmid detection using PlasmidFinder
rule plasmidfinder:
    input:
        contigs="results/mge/assemblies/{sample}.contigs.fasta"
    output:
        tsv="results/mge/plasmids/plasmidfinder/{sample}.plasmids.tsv",
        json="results/mge/plasmids/plasmidfinder/{sample}.plasmids.json"
    params:
        outdir="results/mge/plasmids/plasmidfinder/{sample}",
        db=config["plasmidfinder_db"]
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run PlasmidFinder
        plasmidfinder.py \
            -i {input.contigs} \
            -o {params.outdir} \
            -p {params.db} \
            -mp blast \
            -t {threads}
        
        # Rename output files
        mv {params.outdir}/results_tab.tsv {output.tsv}
        mv {params.outdir}/results.json {output.json}
        """

# Plasmid detection and characterization using Platon
rule platon:
    input:
        contigs="results/mge/assemblies/{sample}.contigs.fasta"
    output:
        summary="results/mge/plasmids/platon/{sample}/platon_summary.tsv",
        plasmids="results/mge/plasmids/platon/{sample}/platon_plasmids.fasta"
    params:
        outdir="results/mge/plasmids/platon/{sample}",
        db=config["platon_db"]
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run Platon
        platon --db {params.db} \
               --threads {threads} \
               --output {params.outdir} \
               --verbose \
               {input.contigs}
        """

# Phage detection using VIBRANT
rule vibrant:
    input:
        contigs="results/mge/assemblies/{sample}.contigs.fasta"
    output:
        phages="results/mge/phages/vibrant/{sample}/VIBRANT_phages_{sample}/VIBRANT_results_{sample}/VIBRANT_phages_{sample}.tsv",
        genome_quality="results/mge/phages/vibrant/{sample}/VIBRANT_phages_{sample}/VIBRANT_results_{sample}/VIBRANT_genome_quality_{sample}.tsv"
    params:
        outdir="results/mge/phages/vibrant/{sample}",
        vibrant_db=config["vibrant_db"]
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run VIBRANT
        VIBRANT_run.py -i {input.contigs} \
                      -t {threads} \
                      -folder {params.outdir} \
                      -d {params.vibrant_db} \
                      -virome
        """

# Phage detection using Phigaro
rule phigaro:
    input:
        contigs="results/mge/assemblies/{sample}.contigs.fasta"
    output:
        tsv="results/mge/phages/phigaro/{sample}/phigaro.tsv",
        html="results/mge/phages/phigaro/{sample}/phigaro.html"
    params:
        outdir="results/mge/phages/phigaro/{sample}",
        config=config["phigaro_config"]
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run Phigaro
        phigaro -f {input.contigs} \
                -o {params.outdir} \
                -e tsv html \
                -t {threads} \
                -c {params.config}
        """

# Transposon detection
rule transposon_detection:
    input:
        contigs="results/mge/assemblies/{sample}.contigs.fasta"
    output:
        tsv="results/mge/transposons/{sample}.transposons.tsv"
    params:
        transposon_db=config["transposon_db"]
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        
        # Extract ORFs using Prodigal
        prodigal -i {input.contigs} -a results/mge/transposons/{wildcards.sample}.proteins.faa -p meta
        
        # Search against transposon database using DIAMOND
        diamond blastp --query results/mge/transposons/{wildcards.sample}.proteins.faa \
                      --db {params.transposon_db} \
                      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
                      --threads {threads} \
                      --evalue 1e-10 \
                      --id 80 \
                      --query-cover 80 \
                      --out {output.tsv}
        """

# Insertion sequence detection
rule is_detection:
    input:
        contigs="results/mge/assemblies/{sample}.contigs.fasta"
    output:
        tsv="results/mge/insertion_sequences/{sample}.is.tsv"
    params:
        is_db=config["insertion_sequence_db"]
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        
        # Search against IS database using BLAST
        blastn -query {input.contigs} \
               -db {params.is_db} \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
               -num_threads {threads} \
               -evalue 1e-10 \
               -perc_identity 80 \
               -out {output.tsv}
        """

# Integrated MGE detection (genomic islands, integrons, etc.)
rule integrated_mge:
    input:
        contigs="results/mge/assemblies/{sample}.contigs.fasta"
    output:
        gi="results/mge/integrated/{sample}.genomic_islands.tsv",
        integrons="results/mge/integrated/{sample}.integrons.tsv",
        combined="results/mge/integrated/{sample}.integrated_mge.tsv"
    params:
        outdir="results/mge/integrated",
        integron_finder_db=config["integron_finder_db"]
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run IslandPath-DIMOB for genomic island prediction
        islandpath-dimob {input.contigs} > {output.gi}
        
        # Run IntegronFinder for integron detection
        integron_finder --cpu {threads} \
                       --dbdir {params.integron_finder_db} \
                       --outdir {params.outdir}/{wildcards.sample}_integrons \
                       --mute \
                       {input.contigs}
        
        # Copy integron results to output
        cp {params.outdir}/{wildcards.sample}_integrons/Results_Integron_Finder_{wildcards.sample}.integrons {output.integrons}
        
        # Combine results
        python scripts/combine_integrated_mge.py \
            --gi {output.gi} \
            --integrons {output.integrons} \
            --output {output.combined} \
            --sample {wildcards.sample}
        """

# Combine all MGE results
rule combine_mge_results:
    input:
        plasmids=expand("results/mge/plasmids/plasmidfinder/{sample}.plasmids.tsv", sample=SAMPLES),
        platon=expand("results/mge/plasmids/platon/{sample}/platon_summary.tsv", sample=SAMPLES),
        phages_vibrant=expand("results/mge/phages/vibrant/{sample}/VIBRANT_phages_{sample}/VIBRANT_results_{sample}/VIBRANT_phages_{sample}.tsv", sample=SAMPLES),
        phages_phigaro=expand("results/mge/phages/phigaro/{sample}/phigaro.tsv", sample=SAMPLES),
        transposons=expand("results/mge/transposons/{sample}.transposons.tsv", sample=SAMPLES),
        is_elements=expand("results/mge/insertion_sequences/{sample}.is.tsv", sample=SAMPLES),
        integrated=expand("results/mge/integrated/{sample}.integrated_mge.tsv", sample=SAMPLES)
    output:
        combined="results/mge/reports/combined_mge_report.tsv"
    script:
        "scripts/combine_mge_results.py"

# Create MGE visualizations
rule visualize_mge:
    input:
        combined="results/mge/reports/combined_mge_report.tsv"
    output:
        heatmap="results/mge/visualization/mge_heatmap.png",
        network="results/mge/visualization/mge_network.png",
        barplot="results/mge/visualization/mge_barplot.png",
        circular="results/mge/visualization/mge_circular.png"
    script:
        "scripts/visualize_mge.py"

# Generate HTML report for MGE analysis
rule mge_report:
    input:
        combined="results/mge/reports/combined_mge_report.tsv",
        heatmap="results/mge/visualization/mge_heatmap.png",
        network="results/mge/visualization/mge_network.png",
        barplot="results/mge/visualization/mge_barplot.png",
        circular="results/mge/visualization/mge_circular.png"
    output:
        report="results/mge/reports/mge_analysis_report.html"
    script:
        "scripts/mge_report.py"
