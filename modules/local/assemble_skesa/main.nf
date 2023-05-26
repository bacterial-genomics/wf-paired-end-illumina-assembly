process ASSEMBLE_SKESA {

    publishDir "${params.outdir}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "${prefix}"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Raw_Assembly_File.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    label "process_high"
    tag { "${prefix}" }

    container "gregorysprenger/skesa@sha256:4455882b5d0fd968630325428729395422be7340301c31d15874a295904b7f26"

    input:
    tuple val(prefix), path(R1), path(R2), path(single), path(qc_nonoverlap_filecheck)

    output:
    tuple val(prefix), path("contigs.fasta"), path("*File*.tsv"), emit: contigs
    path "${prefix}.Raw_Assembly_File.tsv", emit: qc_raw_assembly_filecheck
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

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

    msg "INFO: Running SKESA with !{task.cpus} threads and !{task.memory} memory"

    if [[ ! -f contigs.fasta ]]; then
      RAMSIZE_TOT=$(echo !{task.memory} | cut -d ' ' -f 1)
      msg "INFO: RAMSIZE = ${RAMSIZE_TOT}"

      skesa \
        --reads !{R1},!{R2} \
        --reads !{single} \
        --use_paired_ends \
        --contigs_out contigs.fasta \
        --memory !{task.memory} \
        --cores !{task.cpus} \
        --kmer !{params.skesa_kmer_length} \
        --vector_percent !{params.skesa_vector_percent} \
        --steps !{params.skesa_steps} \
        --fraction !{params.skesa_fraction} \
        --max_snp_len !{params.skesa_max_snp_length} \
        --min_contig !{params.skesa_min_contig_length} \
        !{allow_snps}

    fi

    if verify_minimum_file_size "contigs.fasta" 'Raw Assembly File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{prefix}\tRaw Assembly File\tPASS" > !{prefix}.Raw_Assembly_File.tsv
    else
      echo -e "!{prefix}\tRaw Assembly File\tFAIL" > !{prefix}.Raw_Assembly_File.tsv
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      skesa: $(skesa --version 2>&1 | grep 'SKESA' | cut -d ' ' -f 2)
    END_VERSIONS
    '''
}
