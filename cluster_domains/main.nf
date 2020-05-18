#!/usr/bin/env nextflow

/*
nextflow main.nf -resume
 */

// This file is created by `find_cath_superfamily_domains_in_mgnify`
mgy_domain_sequences = Channel
    .fromPath("${params.project_dir}/mgy_clusters_${params.superfamily}.dom.fa")

// This file is made by a script that reads from the CATH SQL DB
cath_domains_tsv = Channel
    .fromPath("${params.cath_dir}/superfamily/cath_v4_2_0_${params.superfamily}.dom.tsv")

process cath_domain_sequences {
    publishDir params.publish_dir, mode: "copy"

    container "julia:1.2.0-buster"

    input:
    file cath_domains_tsv

    output:
    file "cath_v4_2_0_${params.superfamily}.dom.fa" into cath_domain_sequences

    """
    #!/usr/bin/env julia
    using CATHBase

    records = parse_sf_tsv("$cath_domains_tsv")
    outputfasta = "cath_v4_2_0_${params.superfamily}.dom.fa"
    writefasta(outputfasta, records)
    """
}

process concatenate_mgnify_and_cath_domains {
    publishDir params.publish_dir, mode: "copy"

    container "julia:1.2.0-buster"

    input:
    file cath_domain_sequences
    file mgy_domain_sequences

    output:
    file "cath_mgy_${params.superfamily}.dom.fa" into cath_mgy_domain_sequences

    """
    #!/usr/bin/env julia

    using CATHBase

    cath = readfasta("$cath_domain_sequences")
    mgy = filter!(readfasta("$mgy_domain_sequences"), FullLengthSequence())
    writefasta("cath_mgy_${params.superfamily}.dom.fa", [cath; mgy])
    """

    // """
    // cat $cath_domain_sequences $mgy_domain_sequences > cath_mgy_${params.superfamily}.dom.fa
    // """
}

cath_mgy_domain_sequences = cath_mgy_domain_sequences
    // .splitFasta(by: 1000).take(1) // TEMP
    .first()

seq_id_percentage = Channel.from(30, 50, 70, 90)

process linclust {
    publishDir "${params.publish_dir}/domain_sequence_clusters", mode: "copy"

    container "soedinglab/mmseqs2:version-10"
    cpus 8
    memory "8 GB"

    input:
    file cath_mgy_domain_sequences
    val seq_id_percentage

    output:
    file "cath_mgy_${params.superfamily}.dom.${seq_id_percentage}.clust.tsv"

    script:
    seq_id_decimal = seq_id_percentage / 100

    """
    mmseqs createdb $cath_mgy_domain_sequences seq.db
    mmseqs linclust \
        -c $seq_id_decimal --min-seq-id $seq_id_decimal \
        --cov-mode 1 \
        --kmer-per-seq 100 \
        --threads ${task.cpus} \
        seq.db clust.db tmp
    mmseqs createtsv seq.db seq.db clust.db \
        cath_mgy_${params.superfamily}.dom.${seq_id_percentage}.clust.tsv
    rm -rf tmp
    """
}
