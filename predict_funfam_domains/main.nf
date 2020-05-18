#!/usr/bin/env nextflow

/*
nextflow main.nf -w ~/data/mgy_clusters_3.40.50.1820/funfam/work -resume
*/

Channel
    .fromPath("${params.project_dir}/mgy_clusters_${params.superfamily}.seq.crh.fl.fa")
    .into{ sequences1; sequences2 }

hmms = Channel
    .fromPath("${params.cath_dir}/funfam/funfam-hmm3-v4_2_0.lib")
    .first()

process db_size {
    container "debian:buster-slim"

    input:
    file "sequences" from sequences1

    output:
    file db_size

    """
    grep -c ">" $sequences | tr -d '\n' > db_size
    """
}

db_size = db_size.first()

sequences = sequences2
    .splitFasta(by: params.n_records_per_split, file: true)
    // .take(params.n_splits) // TEMP

process hmmsearch {
    cpus 4
    container "biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1"

    input:
    file sequences
    file hmms
    file db_size

    output:
    file domtbl
    file hmmsearch_out

    """
    hmmsearch --cpu $task.cpus --cut_tc -Z \$(cat $db_size) \
        -o hmmsearch_out --domtblout domtbl $hmms $sequences
    """
}

process cath_resolve_hits {
    container "harryscholes/cath-resolve-hits:0.16.2"

    input:
    file hmmsearch_out

    output:
    file resolvedhits

    """
    cath-resolve-hits --input-format hmmsearch_out --output-hmmer-aln \
        --min-dc-hmm-coverage=80 --worst-permissible-bitscore 25 \
        --long-domains-preference 2 \
         $hmmsearch_out > resolvedhits
    """
}

process merge_hmmsearch_output {
    publishDir params.publish_dir, mode: "copy"

    container "debian:buster-slim"

    input:
    file "domtbl" from domtbl.collect()

    output:
    file "${params.superfamily}.funfam.domtbl.gz" into merged_domtbl

    """
    sed '/^#/d' $domtbl | gzip > ${params.superfamily}.funfam.domtbl.gz
    """
}

process merge_cath_resolve_hits_output {
    publishDir params.publish_dir, mode: "copy"

    container "debian:buster-slim"

    input:
    file "resolvedhits" from resolvedhits.collect()

    output:
    file "${params.superfamily}.funfam.crh.gz" into merged_resolvedhits

    """
    sed '/^#/d' $resolvedhits | gzip > ${params.superfamily}.funfam.crh.gz
    """
}
