#!/usr/bin/env nextflow

/*
================================================================================
    MetaVirome Pipeline - A Nextflow workflow for virus detection in metagenomes
================================================================================
*/

// Define parameters with defaults
params {
    // Input options
    reads          = null
    host_genome    = null
    single_end     = false
    
    // Output directory
    outdir         = "./results"
    
    // Resources
    max_memory     = '128.GB'
    max_cpus       = 16
    max_time       = '240.h'
    
    // Databases
    kraken2_db     = null  // Path to Kraken2 viral database
    diamond_db     = null  // Path to DIAMOND viral protein database
    
    // Tool parameters
    min_contig_len = 1000  // Minimum contig length for viral analysis
    
    // Help message
    help = false
}

// Print help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (params.reads == null) {
    error "Input reads not specified with --reads argument"
    exit 1
}

if (params.host_genome == null) {
    error "Host genome not specified with --host_genome argument"
    exit 1
}

if (params.kraken2_db == null) {
    error "Kraken2 database not specified with --kraken2_db argument"
    exit 1
}

if (params.diamond_db == null) {
    error "DIAMOND database not specified with --diamond_db argument"
    exit 1
}

// Create a channel for input reads
Channel
    .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
    .set { read_pairs_ch }

// Define host genome
host_genome_file = file(params.host_genome)

/*
 * PREPROCESSING - Remove host DNA
 */
process REMOVE_HOST {
    tag "$sample_id"
    publishDir "${params.outdir}/01_host_removed", mode: 'copy'
    
    input:
    tuple val(sample_id), file(reads) from read_pairs_ch
    file host_genome from host_genome_file
    
    output:
    tuple val(sample_id), file("${sample_id}_host_removed_*.fastq.gz") into host_removed_reads_ch
    
    script:
    def read_args = params.single_end ? "-1 ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def out_args = params.single_end ? "-o ${sample_id}_host_removed_1.fastq.gz" : "-o ${sample_id}_host_removed_1.fastq.gz -o2 ${sample_id}_host_removed_2.fastq.gz"
    """
    bowtie2-build $host_genome host_reference
    bowtie2 -p $task.cpus -x host_reference $read_args --un-conc-gz ${sample_id}_host_removed_%.fastq.gz
    """
}

/*
 * KNOWN VIRUS DETECTION - Kraken2
 */
process KRAKEN2 {
    tag "$sample_id"
    publishDir "${params.outdir}/02_kraken2", mode: 'copy'
    
    input:
    tuple val(sample_id), file(reads) from host_removed_reads_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_kraken2.report") into kraken2_reports_ch
    tuple val(sample_id), file(reads) into diamond_reads_ch
    
    script:
    def read_args = params.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    """
    kraken2 --db $params.kraken2_db --threads $task.cpus --output ${sample_id}_kraken2.output --report ${sample_id}_kraken2.report $read_args
    """
}

/*
 * SENSITIVE VIRAL DETECTION - DIAMOND
 */
process DIAMOND {
    tag "$sample_id"
    publishDir "${params.outdir}/03_diamond", mode: 'copy'
    
    input:
    tuple val(sample_id), file(reads) from diamond_reads_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_diamond.tsv") into diamond_results_ch
    tuple val(sample_id), file(reads) into assembly_reads_ch
    
    script:
    def read_input = params.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    """
    if [ "${params.single_end}" == "true" ]; then
        diamond blastx --db $params.diamond_db --query ${reads[0]} --outfmt 6 --out ${sample_id}_diamond.tsv --threads $task.cpus --sensitive
    else
        cat ${reads[0]} ${reads[1]} > ${sample_id}_combined.fastq.gz
        diamond blastx --db $params.diamond_db --query ${sample_id}_combined.fastq.gz --outfmt 6 --out ${sample_id}_diamond.tsv --threads $task.cpus --sensitive
    fi
    """
}

/*
 * VIRAL ASSEMBLY - metaSPAdes
 */
process METASPADES {
    tag "$sample_id"
    publishDir "${params.outdir}/04_assembly", mode: 'copy'
    
    input:
    tuple val(sample_id), file(reads) from assembly_reads_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_contigs.fasta") into assembly_contigs_ch
    
    script:
    def memory = task.memory.toGiga()
    def read_args = params.single_end ? "-s ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    metaspades.py $read_args --meta -t $task.cpus -m $memory -o assembly
    
    # Rename and filter contigs by minimum length
    awk '/^>/ {print ">$sample_id" ++i; next} {print}' < assembly/contigs.fasta | \
    awk 'BEGIN {RS=">"} length(\$0) >= $params.min_contig_len {print ">"substr(\$0,1)}' > ${sample_id}_contigs.fasta
    """
}

/*
 * VIRAL CONTIG IDENTIFICATION - VirSorter2
 */
process VIRSORTER2 {
    tag "$sample_id"
    publishDir "${params.outdir}/05_virsorter2", mode: 'copy'
    
    input:
    tuple val(sample_id), file(contigs) from assembly_contigs_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_virsorter"), file("${sample_id}_viral_contigs.fasta") into virsorter_contigs_ch
    
    script:
    """
    virsorter run --keep-tmpdir --seqfile $contigs --jobs $task.cpus --working-dir ${sample_id}_virsorter --include-groups "dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae"
    
    # Extract viral contigs
    cat ${sample_id}_virsorter/final-viral-combined.fa > ${sample_id}_viral_contigs.fasta
    """
}

/*
 * VIRAL QUALITY CONTROL - CheckV
 */
process CHECKV {
    tag "$sample_id"
    publishDir "${params.outdir}/06_checkv", mode: 'copy'
    
    input:
    tuple val(sample_id), file(virsorter_dir), file(viral_contigs) from virsorter_contigs_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_checkv"), file("${sample_id}_high_quality.fasta") into checkv_contigs_ch
    tuple val(sample_id), file("${sample_id}_high_quality.fasta") into (pharokka_contigs_ch, vcontact_contigs_ch)
    
    script:
    """
    checkv end_to_end $viral_contigs ${sample_id}_checkv -t $task.cpus
    
    # Get high-quality and complete viral contigs
    cat ${sample_id}_checkv/quality_summary.tsv | \
        awk -F'\t' '\$8=="High-quality" || \$8=="Complete" {print \$1}' > high_quality_ids.txt
    
    seqtk subseq $viral_contigs high_quality_ids.txt > ${sample_id}_high_quality.fasta
    """
}

/*
 * VIRAL ANNOTATION - Pharokka (for phages)
 */
process PHAROKKA {
    tag "$sample_id"
    publishDir "${params.outdir}/07_annotation", mode: 'copy'
    
    input:
    tuple val(sample_id), file(high_quality_contigs) from pharokka_contigs_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_pharokka") into pharokka_results_ch
    
    script:
    """
    pharokka.py -i $high_quality_contigs -o ${sample_id}_pharokka -t $task.cpus -f
    """
}

/*
 * NOVELTY ASSESSMENT - vConTACT2
 */
process VCONTACT2 {
    tag "$sample_id"
    publishDir "${params.outdir}/08_taxonomy", mode: 'copy'
    
    input:
    tuple val(sample_id), file(high_quality_contigs) from vcontact_contigs_ch
    
    output:
    tuple val(sample_id), file("${sample_id}_vcontact") into vcontact_results_ch
    
    script:
    """
    # Extract proteins using Prodigal
    prodigal -i $high_quality_contigs -a ${sample_id}_proteins.faa -p meta
    
    # Prepare input file for vConTACT2
    echo "protein_id,contig_id,keywords" > ${sample_id}_gene2genome.csv
    grep ">" ${sample_id}_proteins.faa | sed 's/>//g' | awk -F ' # ' '{print \$1","substr(\$1,1,index(\$1,"_")-1)",,"}' >> ${sample_id}_gene2genome.csv
    
    # Run vConTACT2
    vcontact2 --raw-proteins ${sample_id}_proteins.faa \
              --rel-mode 'Diamond' \
              --proteins-fp ${sample_id}_gene2genome.csv \
              --db 'ProkaryoticViralRefSeq201-Merged' \
              --output-dir ${sample_id}_vcontact \
              --threads $task.cpus
    """
}

/*
 * Function to display help message
 */
def helpMessage() {
    log.info"""
    ==========================================================
      MetaVirome Pipeline - A Nextflow workflow for virus detection in metagenomes
    ==========================================================
    
    Usage:
    nextflow run metavirome.nf --reads "path/to/reads/*_R{1,2}.fastq.gz" --host_genome "path/to/host.fasta" --kraken2_db "path/to/kraken2_db" --diamond_db "path/to/diamond_db"
    
    Required arguments:
      --reads          Path to input reads (must be enclosed in quotes and include the glob pattern)
      --host_genome    Path to host genome for host DNA removal
      --kraken2_db     Path to Kraken2 viral database
      --diamond_db     Path to DIAMOND viral protein database
    
    Optional arguments:
      --outdir         Output directory (default: ./results)
      --single_end     Specify if reads are single-end (default: false)
      --min_contig_len Minimum contig length for viral analysis (default: 1000)
      
    Resource allocation:
      --max_memory     Maximum memory to use (default: 128.GB)
      --max_cpus       Maximum CPUs to use (default: 16)
      --max_time       Maximum time to run (default: 240.h)
    """
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}


