// Nextflow pipeline: Novel Bacteria/Virus Detection from NGS data

params.reads = "data/*_R{1,2}.fastq.gz"
params.ref_host = "reference/hg38.fa"
params.output = "results"

process QC {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastqc.html")

    script:
    """
    fastqc -o . ${reads}
    """
}

process HostRemoval {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_clean_R1.fastq.gz"), path("${sample_id}_clean_R2.fastq.gz")

    script:
    """
    bowtie2 -x ${params.ref_host} -1 ${reads[0]} -2 ${reads[1]} \
        --un-conc-gz ${sample_id}_clean_R%.fastq.gz -S /dev/null
    """
}

process Assembly {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} -o ${sample_id}_megahit
    cp ${sample_id}_megahit/final.contigs.fa ${sample_id}_contigs.fasta
    """
}

process Binning {
    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_bins")

    script:
    """
    metabat2 -i ${contigs} -o ${sample_id}_bins/bin
    """
}

process GenePrediction {
    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_genes.gff")

    script:
    """
    prodigal -i ${contigs} -o ${sample_id}_genes.gff -a ${sample_id}_proteins.faa -p meta
    """
}

process Annotation {
    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}_annotation.tsv")

    script:
    """
    eggnog-mapper -i ${proteins} -o ${sample_id}_annotation --itype proteins
    """
}

process NoveltyDetection {
    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}_novelty_scores.tsv")

    script:
    """
    python3 run_novelty_detection.py --input ${proteins} --output ${sample_id}_novelty_scores.tsv
    """
}

workflow {
    Channel.fromFilePairs(params.reads, flat: true).set { samples }

    samples |>
        QC |>
        HostRemoval |>
        Assembly |>
        Binning |>
        GenePrediction |>
        Annotation & NoveltyDetection
}  


