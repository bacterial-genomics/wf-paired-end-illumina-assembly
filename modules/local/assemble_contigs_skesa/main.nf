process ASSEMBLE_CONTIGS_SKESA {

    publishDir "${params.outdir}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "${meta.id}"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Raw_Assembly_File.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${meta.id}.${task.process}${filename}" }

    label "process_high"
    tag { "${meta.id}" }

    container "gregorysprenger/skesa@sha256:4455882b5d0fd968630325428729395422be7340301c31d15874a295904b7f26"

    input:
    tuple val(meta), path(R1), path(R2), path(single), path(qc_nonoverlap_filecheck)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                       , emit: versions
    path "${meta.id}.Raw_Assembly_File.tsv"                   , emit: qc_raw_assembly_filecheck
    tuple val(meta), path("contigs.fasta"), path("*File*.tsv"), emit: contigs

    shell:
    allow_snps = params.skesa_allow_snps ? "--allow snps" : ""
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_nonoverlap_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    msg "INFO: Assembling contigs using SKESA"

    if [[ ! -f contigs.fasta ]]; then
      skesa \
        --reads !{R1},!{R2} \
        --reads !{single} \
        !{allow_snps} \
        --cores !{task.cpus} \
        --memory !{task.memory} \
        --contigs_out contigs.fasta \
        --steps !{params.skesa_steps} \
        --kmer !{params.skesa_kmer_length} \
        --fraction !{params.skesa_fraction} \
        --max_snp_len !{params.skesa_max_snp_length} \
        --min_contig !{params.skesa_min_contig_length} \
        --vector_percent !{params.skesa_vector_percent}
    fi

    if verify_minimum_file_size "contigs.fasta" 'Raw Assembly File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{meta.id}\tRaw Assembly File\tPASS" > !{meta.id}.Raw_Assembly_File.tsv
    else
      echo -e "!{meta.id}\tRaw Assembly File\tFAIL" > !{meta.id}.Raw_Assembly_File.tsv
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      skesa: $(skesa --version 2>&1 | grep 'SKESA' | cut -d ' ' -f 2)
    END_VERSIONS
    '''
}
