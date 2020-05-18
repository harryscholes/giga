#!/usr/bin/env nextflow

/*
nextflow main.nf -resume
*/

outfile = "mgy_clusters_${params.superfamily}"

/*
Predict CATH superfamily domains
*/

model_to_family_map = Channel
    .fromPath("${params.cath_dir}/hmms/model_to_family_map.csv")

Channel
    .fromPath("${params.cath_dir}/hmms/main.hmm")
    .first() // convert to value channel
    .into{hmms1; hmms2}


process superfamily_s95_models {
    publishDir params.publish_dir, mode: "copy"

    container "biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1"

    input:
    file model_to_family_map
    file hmms1

    output:
    file "${params.superfamily}.hmm" into s95_models

    """
    grep ${params.superfamily} $model_to_family_map \
        | cut -b 2-11 \
        | hmmfetch -f $hmms1 - > ${params.superfamily}.hmm
    """
}

s95_models = s95_models.first()

// source: ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2019_05/mgy_clusters.fa.gz
Channel.fromPath("${params.mgnify_dir}/mgy_clusters_20190531.fa.gz")
    .into{mgnify_fasta1; mgnify_fasta2}

mgnify_fasta1 = mgnify_fasta1
    .splitFasta(by: params.n_records_per_split_sf, file: true)
    // .take(params.n_splits) // TEMP

process hmmsearch_superfamily_s95_models {
    cpus 4
    container "biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1"

    input:
    file mgnify_fasta1
    file s95_models

    output:
    file domtbl

    """
    hmmsearch --cpu $task.cpus --domE 0.001 --incdomE 0.001 -Z $params.mgnify_db_size \
        -o /dev/null --domtblout domtbl $s95_models $mgnify_fasta1
    """
}

process merge_hmmsearch_superfamily_s95_models_output {
    publishDir params.publish_dir, mode: "copy"

    container "debian:buster-slim"

    input:
    file "domtbl" from domtbl.collect()

    output:
    file "${outfile}.domtbl.gz" into merged_domtbl

    """
    sed '/^#/d' $domtbl | gzip > ${outfile}.domtbl.gz
    """
}

process extract_sequences_matching_superfamily_s95_models {
    publishDir params.publish_dir, mode: "copy"

    container "julia:1.2.0-buster"

    input:
    file merged_domtbl
    file mgnify_fasta2

    output:
    file "${outfile}.seq.fa" into mgnify_sf_seq_fasta1, mgnify_sf_seq_fasta2, mgnify_sf_seq_fasta3

    """
    #!/usr/bin/env julia

    using CATHBase

    function main()
        inputfasta = "$mgnify_fasta2"
        domainmatches = matches(DOMTBLFile("$merged_domtbl"))
        outputfasta = "${outfile}.seq.fa"
        sequences(inputfasta, domainmatches, outputfasta)
    end

    main()
    """
}

process db_size {
    container "debian:buster-slim"

    input:
    file mgnify_sf_seq_fasta1

    output:
    file mgnify_sf_seq_db_size

    """
    grep -c ">" $mgnify_sf_seq_fasta1 | tr -d '\n' > mgnify_sf_seq_db_size
    """
}

mgnify_sf_seq_fasta2 = mgnify_sf_seq_fasta2
    .splitFasta(by: params.n_records_per_split_all, file: true)
    // .take(params.n_splits) // TEMP

mgnify_sf_seq_db_size = mgnify_sf_seq_db_size.first() // convert to value channel

process hmmsearch_all_s95_models {
    cpus 4
    container "biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1"

    input:
    file mgnify_sf_seq_fasta2
    file mgnify_sf_seq_db_size
    file hmms2

    output:
    file domtbl2
    file hmmsearchf

    """
    hmmsearch --cpu $task.cpus --domE 0.001 --incdomE 0.001 -Z \$(cat $mgnify_sf_seq_db_size) \
        -o hmmsearchf --domtblout domtbl2 $hmms2 $mgnify_sf_seq_fasta2
    """
}

process cath_resolve_hits {
    container "harryscholes/cath-resolve-hits:0.16.2"

    input:
    file hmmsearchf

    output:
    file "resolvedhits" into resolvedhits1, resolvedhits2

    """
    cath-resolve-hits --min-dc-hmm-coverage=80 --worst-permissible-bitscore 25 \
        --output-hmmer-aln --input-format hmmsearch_out $hmmsearchf > resolvedhits
    """
}

process assign_cath_superfamilies {
    container "python:2.7.15"

    input:
    file resolvedhits2

    output:
    file assigned_superfamilies

    """
    python $baseDir/../gene3d/assign_cath_superfamilies.py $resolvedhits2
    cat ${resolvedhits2}.csv > assigned_superfamilies
    """
}


process merge_hmmsearch_all_s95_models_output {
    publishDir params.publish_dir, mode: "copy"

    container "debian:buster-slim"

    input:
    file "domtbl2" from domtbl2.collect()

    output:
    file "${outfile}_all.domtbl.gz" into merged_domtbl2

    """
    sed '/^#/d' $domtbl2 | gzip > ${outfile}_all.domtbl.gz
    """
}

process merge_cath_resolve_hits_output {
    publishDir params.publish_dir, mode: "copy"

    container "debian:buster-slim"

    input:
    file "resolvedhits" from resolvedhits1.collect()

    output:
    file "${outfile}_all.crh.gz" into merged_resolvedhits

    """
    sed '/^#/d' $resolvedhits | gzip > ${outfile}_all.crh.gz
    """
}

process merge_assign_cath_superfamilies_output {
    publishDir params.publish_dir, mode: "copy"

    container "debian:buster-slim"

    input:
    file "assigned_superfamilies" from assigned_superfamilies.collect()

    output:
    file "${outfile}_all.mda.gz" into merged_assigned_superfamilies

    """
    sed '/^#/d' $assigned_superfamilies | gzip > ${outfile}_all.mda.gz
    """
}

process extract_domain_sequences {
    publishDir params.publish_dir, mode: "copy"

    container "julia:1.2.0-buster"

    input:
    file merged_assigned_superfamilies
    file mgnify_sf_seq_fasta3

    output:
    file "${outfile}.dom.fa" into mgy_domain_sequences

    """
    #!/usr/bin/env julia

    using CATHBase

    function main()
        inputfasta = "$mgnify_sf_seq_fasta3"
        domainmatches = matches(MDAFile("$merged_assigned_superfamilies"))
        filter!(x->model(x) == "${params.superfamily}", domainmatches)
        outputfasta = "${outfile}.dom.fa"
        domains(inputfasta, domainmatches, outputfasta)
    end

    main()
    """
}
