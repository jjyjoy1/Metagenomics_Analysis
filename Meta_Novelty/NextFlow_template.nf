// Anomaly detection using embeddings
process AnomalyDetection {
    publishDir "${params.output}/anomaly_detection", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(protein_embeddings), path(protein_ids)
    tuple val(sample_id), path(dna_embeddings), path(contig_ids)
    tuple val(sample_id), path(cat_dir), path(cat_taxonomy)
    tuple val(sample_id), path(mmseqs2_dir), path(mmseqs2_taxonomy)
    tuple val(sample_id), path(checkv_dir), path(combined_viral), path(high_conf_viral)
    
    output:
    tuple val(sample_id), path("${sample_id}_anomalies_protein.tsv"), path("${sample_id}_anomalies_dna.tsv"), path("${sample_id}_anomaly_models"), path("${sample_id}_visualization")
    
    script:
    """
    # Create output directories
    mkdir -p ${sample_id}_anomaly_models
    mkdir -p ${sample_id}_visualization
    
    # Prepare annotation files for contextual anomaly detection
    python3 -c "
    import pandas as pd
    import numpy as np
    
    # Load CAT taxonomy
    try:
        cat_df = pd.read_csv('${cat_taxonomy}', sep='\t')
        cat_df = cat_df.rename(columns={'# contig': 'contig_id'})
        cat_dict = dict(zip(cat_df['contig_id'], cat_df['classification']))
    except:
        cat_dict = {}
    
    # Load MMseqs2 taxonomy
    try:
        mmseqs_df = pd.read_csv('${mmseqs2_taxonomy}', sep='\t', header=None)
        mmseqs_df.columns = ['query', 'taxonomy', 'rank', 'taxid', 'score']
        mmseqs_dict = dict(zip(mmseqs_df['query'], mmseqs_df['taxonomy']))
    except:
        mmseqs_dict = {}
    
    # Load viral contigs
    viral_contigs = set()
    try:
        with open('${combined_viral}') as f:
            for line in f:
                if line.startswith('>'):
                    contig = line.strip().lstrip('>')
                    viral_contigs.add(contig)
    except:
        pass
    
    # Load protein IDs and contig IDs
    protein_ids = []
    with open('${protein_ids}') as f:
        protein_ids = [line.strip() for line in f]
    
    contig_ids = []
    with open('${contig_ids}') as f:
        contig_ids = [line.strip() for line in f]
    
    # Create annotation files
    protein_anno = {}
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
    
    contig_anno = {}
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
            f.write(f'{pid}\\t{anno}\\n')
    
    with open('contig_annotations.tsv', 'w') as f:
        f.write('contig_id\\ttaxonomy\\n')
        for cid, anno in contig_anno.items():
            f.write(f'{cid}\\t{anno}\\n')
    "
    
    # Run anomaly detection using both protein and DNA embeddings
    python3 run_anomaly_detection.py \
        --protein_embeddings ${protein_embeddings} \
        --protein_ids ${protein_ids} \
        --protein_annotations protein_annotations.tsv \
        --dna_embeddings ${dna_embeddings} \
        --contig_ids ${contig_ids} \
        --contig_annotations contig_annotations.tsv \
        --methods "${params.anomaly_detection_methods}" \
        --output_protein ${sample_id}_anomalies_protein.tsv \
        --output_dna ${sample_id}_anomalies_dna.tsv \
        --save_models ${sample_id}_anomaly_models \
        --visualize ${sample_id}_visualization \
        --threads ${params.threads} \
        --contamination 0.05 \
        --dimensions_reduction umap,tsne \
        --perplexity 30
    """
}

// Generate comprehensive report integrating all results
process ComprehensiveReport {
    publishDir "${params.output}/reports", mode: 'copy'
    
    input:
    tuple val(sample_id), path(anomalies_protein), path(anomalies_dna), path(anomaly_models), path(visualization)
    tuple val(sample_id), path(cat_dir), path(cat_taxonomy)
    tuple val(sample_id), path(gtdbtk_dir)
    tuple val(sample_id), path(ani_results)
    tuple val(sample_id), path(checkv_dir), path(combined_viral), path(high_conf_viral)
    tuple val(sample_id), path(vcontact2_dir)
    
    output:
    tuple val(sample_id), path("${sample_id}_comprehensive_report.html"), path("${sample_id}_novel_pathogens.tsv")
    
    script:
    """
    # Generate comprehensive report
    python3 generate_comprehensive_report.py \
        --sample_id ${sample_id} \
        --anomalies_protein ${anomalies_protein} \
        --anomalies_dna ${anomalies_dna} \
        --visualization_dir ${visualization} \
        --cat_taxonomy ${cat_taxonomy} \
        --gtdbtk_dir ${gtdbtk_dir} \
        --ani_dir ${ani_results} \
        --checkv_dir ${checkv_dir} \
        --vcontact2_dir ${vcontact2_dir} \
        --high_confidence_viral ${high_conf_viral} \
        --output_report ${sample_id}_comprehensive_report.html \
        --output_novel_pathogens ${sample_id}_novel_pathogens.tsv
    """
}

// Main workflow
workflow {
    // Create output directory
    createOutputDir()
    
    // Create input channel from paired reads
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, reads -> tuple(sample_id, reads) }
    
    // Run QC on raw reads
    qc_ch = QC(read_pairs_ch)
    
    // Remove host sequences
    clean_reads_ch = HostRemoval(read_pairs_ch)
    
    // Assemble non-host reads
    contigs_ch = Assembly(clean_reads_ch)
    
    // Map reads back to contigs
    mapping_ch = ReadMapping(contigs_ch, clean_reads_ch)
    
    // Bin contigs into potential genomes
    bins_ch = Binning(contigs_ch, mapping_ch)
    
    // Improve binning using assembly graph
    graphbin_ch = GraphBin(contigs_ch, bins_ch)
    
    // Predict genes from contigs
    genes_ch = GenePrediction(contigs_ch)
    
    // Run different pathogen detection tools in parallel
    
    // 1. General taxonomic classification tools
    cat_ch = CAT(contigs_ch, genes_ch)
    bat_ch = BAT(graphbin_ch, genes_ch)
    mmseqs2_ch = MMseqs2(contigs_ch, genes_ch)
    
    // 2. Bacterial specific analysis
    gtdbtk_ch = GTDBTk(graphbin_ch)
    ani_ch = ANI(graphbin_ch, gtdbtk_ch)
    
    // 3. Viral specific tools
    virsorter2_ch = VirSorter2(contigs_ch)
    vibrant_ch = VIBRANT(contigs_ch)
    dvf_ch = DeepVirFinder(contigs_ch)
    
    // 4. Combine viral predictions and assess quality
    checkv_ch = CheckV(virsorter2_ch, vibrant_ch, dvf_ch)
    
    // 5. Viral classification
    vcontact2_ch = vConTACT2(checkv_ch, genes_ch)
    
    // 6. Generate embeddings for anomaly detection
    protein_embeddings_ch = ProteinEmbeddings(genes_ch)
    dna_embeddings_ch = DNAEmbeddings(contigs_ch)
    
    // 7. Anomaly detection
    anomalies_ch = AnomalyDetection(
        protein_embeddings_ch,
        dna_embeddings_ch,
        cat_ch,
        mmseqs2_ch,
        checkv_ch
    )
    
    // 8. Generate comprehensive report
    report_ch = ComprehensiveReport(
        anomalies_ch,
        cat_ch,
        gtdbtk_ch,
        ani_ch,
        checkv_ch,
        vcontact2_ch
    )
} Nextflow pipeline: Novel Bacteria/Virus Detection from NGS data

// Define parameters
params.reads = "data/*_R{1,2}.fastq.gz"
params.ref_host = "reference/hg38.fa"
params.output = "results"
params.protein_embedding_model = "esm2" // Options: esm2, prot_t5, esm1b
params.dna_embedding_model = "dnabert" // Options: dnabert, nucleotide_transformer, k_mer
params.anomaly_detection_methods = "isolation_forest,vae,dbscan" // Comma-separated list of methods to use
params.threads = 8 // Number of threads to use for computation-heavy processes
params.min_contig_size = 1000 // Minimum contig size for downstream analysis
params.checkv_db = "/path/to/checkv_db" // Path to CheckV database
params.gtdb_db = "/path/to/gtdb_db" // Path to GTDB database
params.vcontact2_db = "/path/to/vcontact2_db" // Path to vConTACT2 database

// Create output directory
process createOutputDir {
    exec:
    new File(params.output).mkdirs()
}

// Quality control for raw reads
process QC {
    publishDir "${params.output}/qc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastqc.html"), path("${sample_id}_fastqc.zip")

    script:
    """
    fastqc -o . ${reads[0]} ${reads[1]}
    mv *_fastqc.html ${sample_id}_fastqc.html
    mv *_fastqc.zip ${sample_id}_fastqc.zip
    """
}

// Remove host sequences
process HostRemoval {
    publishDir "${params.output}/cleaned_reads", mode: 'copy'
    
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

// Assemble non-host reads
process Assembly {
    publishDir "${params.output}/assembly", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta"), path("${sample_id}_assembly_graph.gfa")

    script:
    """
    megahit -1 ${r1} -2 ${r2} -o ${sample_id}_megahit --num-cpu-threads ${params.threads}
    cp ${sample_id}_megahit/final.contigs.fa ${sample_id}_contigs.fasta
    
    # Generate assembly graph in GFA format for GraphBin
    # Note: If MEGAHIT doesn't produce GFA, might need to use Bandage to convert
    python3 -c "
    import os
    from Bio import SeqIO
    
    # Create a minimal GFA if MEGAHIT doesn't provide one
    if not os.path.exists('${sample_id}_megahit/assembly_graph.gfa'):
        with open('${sample_id}_assembly_graph.gfa', 'w') as f:
            f.write('H\\tVN:Z:1.0\\n')
            for record in SeqIO.parse('${sample_id}_contigs.fasta', 'fasta'):
                f.write(f'S\\t{record.id}\\t{record.seq}\\tLN:i:{len(record.seq)}\\n')
    else:
        os.system('cp ${sample_id}_megahit/assembly_graph.gfa ${sample_id}_assembly_graph.gfa')
    "
    
    # Filter contigs by minimum size
    python3 -c "
    from Bio import SeqIO
    records = [r for r in SeqIO.parse('${sample_id}_contigs.fasta', 'fasta') if len(r.seq) >= ${params.min_contig_size}]
    SeqIO.write(records, '${sample_id}_contigs.fasta', 'fasta')
    "
    """
}

// Map reads back to contigs for coverage information
process ReadMapping {
    publishDir "${params.output}/mapping", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}_coverage.txt")

    script:
    """
    # Index contigs
    bwa index ${contigs}
    
    # Map reads to contigs
    bwa mem -t ${params.threads} ${contigs} ${r1} ${r2} | samtools view -bS - > ${sample_id}.bam
    
    # Sort and index BAM file
    samtools sort -@ ${params.threads} ${sample_id}.bam -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    
    # Calculate coverage
    samtools depth -a ${sample_id}.bam | awk '{sum[\$1]+=\$3; cnt[\$1]++} END {for (contig in sum) print contig, sum[contig]/cnt[contig]}' > ${sample_id}_coverage.txt
    """
}

// Bin contigs into potential genomes
process Binning {
    publishDir "${params.output}/bins", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)
    tuple val(sample_id), path(bam), path(bai), path(coverage)

    output:
    tuple val(sample_id), path("${sample_id}_bins"), path("${sample_id}_bins.summary.tsv")

    script:
    """
    # Run metabat2 for binning
    mkdir -p ${sample_id}_bins
    metabat2 -i ${contigs} -a ${coverage} -o ${sample_id}_bins/bin \
        -m ${params.min_contig_size} -t ${params.threads} --unbinned
    
    # Create a summary of the bins
    echo -e "Bin_Id\\tContigs\\tTotal_Length\\tN50" > ${sample_id}_bins.summary.tsv
    for bin in ${sample_id}_bins/bin*.fa; do
        bin_id=\$(basename \${bin%.fa})
        contig_count=\$(grep -c ">" \$bin)
        total_length=\$(grep -v ">" \$bin | tr -d '\\n' | wc -c)
        
        # Calculate N50
        python3 -c "
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
        n50=\$(cat n50.txt)
        
        echo -e "\$bin_id\\t\$contig_count\\t\$total_length\\t\$n50" >> ${sample_id}_bins.summary.tsv
    done
    """
}

// Improve binning using assembly graph
process GraphBin {
    publishDir "${params.output}/graphbin", mode: 'copy'
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)
    tuple val(sample_id), path(bins), path(bin_summary)

    output:
    tuple val(sample_id), path("${sample_id}_graphbin"), path("${sample_id}_graphbin_bins.csv")

    script:
    """
    # Prepare the initial binning result in the format required by GraphBin
    python3 -c "
    import os
    import glob
    
    bins_path = '${bins}'
    bin_files = glob.glob(os.path.join(bins_path, 'bin*.fa'))
    
    with open('initial_contig_bins.csv', 'w') as outfile:
        outfile.write('contig_id,bin_id\\n')
        for bin_file in bin_files:
            bin_name = os.path.basename(bin_file).replace('.fa', '')
            with open(bin_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        contig_id = line.strip()[1:].split()[0]
                        outfile.write(f'{contig_id},{bin_name}\\n')
    "
    
    # Run GraphBin
    graphbin --assembler megahit --graph ${assembly_graph} --contigs ${contigs} \
        --binned initial_contig_bins.csv --output ${sample_id}_graphbin
    
    # Reorganize bins based on GraphBin output
    mkdir -p ${sample_id}_graphbin
    python3 -c "
    import os
    import shutil
    from Bio import SeqIO
    
    # Read GraphBin output
    graphbin_result = {}
    with open('${sample_id}_graphbin_bins.csv', 'r') as f:
        next(f)  # Skip header
        for line in f:
            contig_id, bin_id = line.strip().split(',')
            if bin_id != 'unbinned':
                graphbin_result[contig_id] = bin_id
    
    # Create new bins
    os.makedirs('${sample_id}_graphbin', exist_ok=True)
    bin_contents = {}
    
    for record in SeqIO.parse('${contigs}', 'fasta'):
        contig_id = record.id
        if contig_id in graphbin_result:
            bin_id = graphbin_result[contig_id]
            if bin_id not in bin_contents:
                bin_contents[bin_id] = []
            bin_contents[bin_id].append(record)
    
    # Write new bin files
    for bin_id, records in bin_contents.items():
        SeqIO.write(records, os.path.join('${sample_id}_graphbin', f'{bin_id}.fa'), 'fasta')
    "
    """
}

// Predict genes from contigs
process GenePrediction {
    publishDir "${params.output}/genes", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)

    output:
    tuple val(sample_id), path("${sample_id}_genes.gff"), path("${sample_id}_proteins.faa"), path("${sample_id}_genes.fna")

    script:
    """
    prodigal -i ${contigs} -o ${sample_id}_genes.gff -a ${sample_id}_proteins.faa -d ${sample_id}_genes.fna -p meta -q
    """
}

// CAT/BAT taxonomic classification
process CAT {
    publishDir "${params.output}/cat", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)
    tuple val(sample_id), path(genes_gff), path(proteins), path(genes_fna)

    output:
    tuple val(sample_id), path("${sample_id}_cat"), path("${sample_id}_cat.taxonomy.tsv")

    script:
    """
    # Run CAT on contigs
    CAT contigs -c ${contigs} -p ${proteins} -o ${sample_id}_cat -t ${params.threads}
    
    # Generate taxonomy file
    CAT add_names -i ${sample_id}_cat.contig2classification.txt -o ${sample_id}_cat.taxonomy.tsv
    """
}

// BAT for bin classification
process BAT {
    publishDir "${params.output}/bat", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(graphbin_dir), path(graphbin_csv)
    tuple val(sample_id), path(genes_gff), path(proteins), path(genes_fna)

    output:
    tuple val(sample_id), path("${sample_id}_bat"), path("${sample_id}_bat.taxonomy.tsv")

    script:
    """
    # Run BAT on bins
    BAT bins -i ${graphbin_dir} -p ${proteins} -o ${sample_id}_bat -t ${params.threads}
    
    # Generate taxonomy file
    CAT add_names -i ${sample_id}_bat.bin2classification.txt -o ${sample_id}_bat.taxonomy.tsv
    """
}

// GTDB-Tk for bacterial phylogenetic placement
process GTDBTk {
    publishDir "${params.output}/gtdbtk", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(graphbin_dir), path(graphbin_csv)

    output:
    tuple val(sample_id), path("${sample_id}_gtdbtk")

    script:
    """
    export GTDBTK_DATA_PATH=${params.gtdb_db}
    
    # Run GTDB-Tk classify workflow
    gtdbtk classify_wf --genome_dir ${graphbin_dir} \
        --out_dir ${sample_id}_gtdbtk \
        --cpus ${params.threads} \
        --extension fa
    """
}

// Calculate Average Nucleotide Identity (ANI)
process ANI {
    publishDir "${params.output}/ani", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(graphbin_dir), path(graphbin_csv)
    tuple val(sample_id), path(gtdbtk_dir)
    
    output:
    tuple val(sample_id), path("${sample_id}_ani_results")
    
    script:
    """
    mkdir -p ${sample_id}_ani_results
    
    # Extract reference genomes based on GTDB-Tk results
    python3 -c "
    import os
    import pandas as pd
    import glob
    
    # Check if classification file exists
    bac_file = '${gtdbtk_dir}/classify/gtdbtk.bac120.summary.tsv'
    ar_file = '${gtdbtk_dir}/classify/gtdbtk.ar53.summary.tsv'
    
    ref_genomes = {}
    
    if os.path.exists(bac_file):
        bac_df = pd.read_csv(bac_file, sep='\t')
        for idx, row in bac_df.iterrows():
            bin_id = row['user_genome']
            closest_ref = row['closest_reference']
            ref_genomes[bin_id] = closest_ref
    
    if os.path.exists(ar_file):
        ar_df = pd.read_csv(ar_file, sep='\t')
        for idx, row in ar_df.iterrows():
            bin_id = row['user_genome']
            closest_ref = row['closest_reference']
            ref_genomes[bin_id] = closest_ref
    
    with open('reference_genomes.txt', 'w') as f:
        for bin_id, ref in ref_genomes.items():
            f.write(f'{bin_id}\t{ref}\\n')
    "
    
    # Run FastANI for each bin against its closest reference
    while read bin_id ref_genome; do
        # We would need to have access to the GTDB reference genomes
        # For now, we'll create a placeholder
        fastANI -q ${graphbin_dir}/${bin_id}.fa -r /path/to/gtdb_genomes/${ref_genome}.fna \
            -o ${sample_id}_ani_results/${bin_id}_vs_${ref_genome}.ani
    done < reference_genomes.txt
    
    # Create a summary of ANI results
    echo -e "Bin_ID\tReference_Genome\tANI\tQuery_Coverage\tReference_Coverage" > ${sample_id}_ani_results/ani_summary.tsv
    for ani_file in ${sample_id}_ani_results/*.ani; do
        if [ -s \$ani_file ]; then
            bin_id=\$(basename \$ani_file | cut -d_ -f1)
            ref_id=\$(basename \$ani_file | cut -d_ -f3 | cut -d. -f1)
            ani_value=\$(cut -f3 \$ani_file)
            query_cov=\$(cut -f4 \$ani_file)
            ref_cov=\$(cut -f5 \$ani_file)
            echo -e "\$bin_id\t\$ref_id\t\$ani_value\t\$query_cov\t\$ref_cov" >> ${sample_id}_ani_results/ani_summary.tsv
        else
            bin_id=\$(basename \$ani_file | cut -d_ -f1)
            ref_id=\$(basename \$ani_file | cut -d_ -f3 | cut -d. -f1)
            echo -e "\$bin_id\t\$ref_id\tNA\tNA\tNA" >> ${sample_id}_ani_results/ani_summary.tsv
        fi
    done
    """
}

// Viral detection using VirSorter2
process VirSorter2 {
    publishDir "${params.output}/virsorter2", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)

    output:
    tuple val(sample_id), path("${sample_id}_virsorter2"), path("${sample_id}_viral_contigs_vs2.fasta")

    script:
    """
    # Run VirSorter2
    virsorter run --seqfile ${contigs} --working-dir ${sample_id}_virsorter2 \
        --include-groups dsDNAphage,ssDNA,RNA --min-length ${params.min_contig_size} \
        -j ${params.threads} --keep-original-seq
    
    # Extract viral sequences
    python3 -c "
    import os
    import pandas as pd
    from Bio import SeqIO
    
    # Read the VirSorter2 results
    results_file = '${sample_id}_virsorter2/final-viral-score.tsv'
    if os.path.exists(results_file):
        df = pd.read_csv(results_file, sep='\t')
        viral_contigs = df[df['max_score'] >= 0.9]['seqname'].tolist()
    else:
        viral_contigs = []
    
    # Extract the viral sequences
    with open('${contigs}', 'r') as input_handle:
        seqs = SeqIO.to_dict(SeqIO.parse(input_handle, 'fasta'))
        
    viral_seqs = [seqs[contig] for contig in viral_contigs if contig in seqs]
    
    with open('${sample_id}_viral_contigs_vs2.fasta', 'w') as output_handle:
        SeqIO.write(viral_seqs, output_handle, 'fasta')
    "
    """
}

// Viral detection using VIBRANT
process VIBRANT {
    publishDir "${params.output}/vibrant", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)

    output:
    tuple val(sample_id), path("${sample_id}_vibrant"), path("${sample_id}_viral_contigs_vibrant.fasta")

    script:
    """
    # Run VIBRANT
    VIBRANT_run.py -i ${contigs} -folder ${sample_id}_vibrant -t ${params.threads} \
        -l ${params.min_contig_size}
    
    # Extract viral sequences
    if [ -f "${sample_id}_vibrant/VIBRANT_phages_${contigs}/VIBRANT_phages_${contigs}.fasta" ]; then
        cp "${sample_id}_vibrant/VIBRANT_phages_${contigs}/VIBRANT_phages_${contigs}.fasta" "${sample_id}_viral_contigs_vibrant.fasta"
    else
        touch "${sample_id}_viral_contigs_vibrant.fasta"
    fi
    """
}

// Additional viral detection with DeepVirFinder
process DeepVirFinder {
    publishDir "${params.output}/deepvirfinder", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)

    output:
    tuple val(sample_id), path("${sample_id}_dvf"), path("${sample_id}_viral_contigs_dvf.fasta")

    script:
    """
    # Run DeepVirFinder
    mkdir -p ${sample_id}_dvf
    python /path/to/DeepVirFinder/dvf.py -i ${contigs} -o ${sample_id}_dvf -c ${params.threads} -l ${params.min_contig_size}
    
    # Extract viral sequences with score > 0.9 and p-value < 0.05
    python3 -c "
    import pandas as pd
    from Bio import SeqIO
    
    # Read DeepVirFinder results
    result_file = '${sample_id}_dvf/${contigs}.dvfpred.txt'
    try:
        df = pd.read_csv(result_file, sep='\t')
        viral_contigs = df[(df['score'] > 0.9) & (df['pvalue'] < 0.05)]['name'].tolist()
    except:
        viral_contigs = []
    
    # Extract viral sequences
    with open('${contigs}', 'r') as input_handle:
        seqs = SeqIO.to_dict(SeqIO.parse(input_handle, 'fasta'))
        
    viral_seqs = [seqs[contig] for contig in viral_contigs if contig in seqs]
    
    with open('${sample_id}_viral_contigs_dvf.fasta', 'w') as output_handle:
        SeqIO.write(viral_seqs, output_handle, 'fasta')
    "
    """
}

// Combine viral predictions and assess quality with CheckV
process CheckV {
    publishDir "${params.output}/checkv", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vs2_dir), path(vs2_fasta)
    tuple val(sample_id), path(vibrant_dir), path(vibrant_fasta)
    tuple val(sample_id), path(dvf_dir), path(dvf_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_checkv"), path("${sample_id}_combined_viral_contigs.fasta"), path("${sample_id}_high_confidence_viral_contigs.fasta")

    script:
    """
    # Combine viral predictions
    cat ${vs2_fasta} ${vibrant_fasta} ${dvf_fasta} > ${sample_id}_all_viral_predictions.fasta
    
    # Remove duplicates
    python3 -c "
    from Bio import SeqIO
    
    # Read all sequences
    seqs = {}
    for record in SeqIO.parse('${sample_id}_all_viral_predictions.fasta', 'fasta'):
        if record.id not in seqs:
            seqs[record.id] = record
    
    # Write unique sequences
    with open('${sample_id}_combined_viral_contigs.fasta', 'w') as outf:
        SeqIO.write(seqs.values(), outf, 'fasta')
    "
    
    # Run CheckV
    checkv end_to_end ${sample_id}_combined_viral_contigs.fasta ${sample_id}_checkv \
        -t ${params.threads} -d ${params.checkv_db}
    
    # Extract high-confidence viral contigs
    python3 -c "
    import pandas as pd
    from Bio import SeqIO
    
    # Read CheckV quality results
    quality_file = '${sample_id}_checkv/quality_summary.tsv'
    try:
        df = pd.read_csv(quality_file, sep='\t')
        high_confidence = df[(df['checkv_quality'] == 'High-quality') | 
                            (df['checkv_quality'] == 'Complete')]['contig_id'].tolist()
    except:
        high_confidence = []
    
    # Extract high-confidence sequences
    with open('${sample_id}_combined_viral_contigs.fasta', 'r') as input_handle:
        seqs = SeqIO.to_dict(SeqIO.parse(input_handle, 'fasta'))
        
    high_conf_seqs = [seqs[contig] for contig in high_confidence if contig in seqs]
    
    with open('${sample_id}_high_confidence_viral_contigs.fasta', 'w') as output_handle:
        SeqIO.write(high_conf_seqs, output_handle, 'fasta')
    "
    """
}

// vConTACT2 for viral classification and networks
process vConTACT2 {
    publishDir "${params.output}/vcontact2", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(checkv_dir), path(combined_viral), path(high_conf_viral)
    tuple val(sample_id), path(genes_gff), path(proteins), path(genes_fna)

    output:
    tuple val(sample_id), path("${sample_id}_vcontact2")

    script:
    """
    # Extract viral proteins
    python3 -c "
    from Bio import SeqIO
    
    # Get viral contig IDs
    viral_ids = set()
    for record in SeqIO.parse('${combined_viral}', 'fasta'):
        viral_ids.add(record.id)
    
    # Extract viral proteins
    viral_proteins = []
    for record in SeqIO.parse('${proteins}', 'fasta'):
        contig_id = record.id.split('_')[0]
        if contig_id in viral_ids:
            viral_proteins.append(record)
    
    # Write viral proteins
    with open('viral_proteins.faa', 'w') as outf:
        SeqIO.write(viral_proteins, outf, 'fasta')
    "
    
    # Run vConTACT2
    vcontact2 --raw-proteins viral_proteins.faa --rel-mode 'Diamond' \
        --proteins-fp viral_proteins.faa --db ${params.vcontact2_db} \
        --output-dir ${sample_id}_vcontact2 --threads ${params.threads}
    """
}

// MMseqs2 clustering and classification
process MMseqs2 {
    publishDir "${params.output}/mmseqs2", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)
    tuple val(sample_id), path(genes_gff), path(proteins), path(genes_fna)

    output:
    tuple val(sample_id), path("${sample_id}_mmseqs2"), path("${sample_id}_mmseqs2_taxonomy.tsv")

    script:
    """
    # Create MMseqs2 database
    mkdir -p ${sample_id}_mmseqs2
    mmseqs createdb ${proteins} ${sample_id}_mmseqs2/query_db
    
    # Run MMseqs2 taxonomy assignment
    mmseqs taxonomy ${sample_id}_mmseqs2/query_db /path/to/mmseqs2_db/mmseqs2_db \
        ${sample_id}_mmseqs2/taxonomy_results tmp \
        --threads ${params.threads} --tax-lineage true
    
    # Get taxonomy report
    mmseqs createtsv ${sample_id}_mmseqs2/query_db ${sample_id}_mmseqs2/taxonomy_results \
        ${sample_id}_mmseqs2_taxonomy.tsv --threads ${params.threads}
    """
}

// Generate protein embeddings for further analysis
process ProteinEmbeddings {
    publishDir "${params.output}/embeddings", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(genes_gff), path(proteins), path(genes_fna)

    output:
    tuple val(sample_id), path("${sample_id}_protein_embeddings.npy"), path("${sample_id}_protein_ids.txt")

    script:
    """
    python3 run_protein_embedding.py \
        --input ${proteins} \
        --model ${params.protein_embedding_model} \
        --output ${sample_id}_protein_embeddings.npy \
        --ids ${sample_id}_protein_ids.txt \
        --threads ${params.threads}
    """
}

// Generate DNA sequence embeddings
process DNAEmbeddings {
    publishDir "${params.output}/embeddings", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(contigs), path(assembly_graph)

    output:
    tuple val(sample_id), path("${sample_id}_dna_embeddings.npy"), path("${sample_id}_contig_ids.txt")

    script:
    """
    # Extract sequences in the right format for embedding model
    python3 -c "
    from Bio import SeqIO
    import numpy as np
    
    seqs = []
    ids = []
    for record in SeqIO.parse('${contigs}', 'fasta'):
        seqs.append(str(record.seq))
        ids.append(record.id)
    
    # Save IDs for later reference
    with open('${sample_id}_contig_ids.txt', 'w') as f:
        f.write('\\n'.join(ids))
    
    # Save sequences for embedding processing
    with open('sequences.txt', 'w') as f:
        f.write('\\n'.join(seqs))
    "
    
    # Run DNA embedding generation
    python3 run_dna_embedding.py \
        --input sequences.txt \
        --model ${params.dna_embedding_model} \
        --output ${sample_id}_dna_embeddings.npy \
        --threads ${params.threads} \
        --kmer_size 6
    """
}

//


