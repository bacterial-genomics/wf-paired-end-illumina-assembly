process ALIGN_16S_BLAST {

    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.16S_BLASTn_Output_File.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    tag { "${prefix}" }

    container "gregorysprenger/ncbi-blast-plus@sha256:2d3e226d2eb31e3e0d5a80d7325b3a2ffd873ad1f2bd81215fd0b43727019279"

    input:
    tuple val(prefix), path(extracted_base), path(qc_extracted_filecheck), path(assembly)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                              , emit: versions
    path "${prefix}.16S_BLASTn_Output_File.tsv"                      , emit: qc_blastn_filecheck
    tuple val(prefix), path("${prefix}.blast.tsv"), path("*File.tsv"), emit: blast_tsv

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_extracted_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Classify each 16S sequence record
    if [[ -d "!{params.blast_db}" ]]; then
      database="!{params.blast_db}"
      msg "INFO: Using user specified BLAST database: !{params.blast_db}"
    else
      database="/db"
      msg "INFO: Using pre-loaded 16S rRNA database for BLAST"
    fi

    # Set BLAST database as an environment variable
    export BLASTDB=${database}

    # Confirm the 16S db exists
    for ext in nin nsq nhr; do
      if ! verify_minimum_file_size "${BLASTDB}/16S_ribosomal_RNA.${ext}" '16S BLASTn database' "!{params.min_filesize_blastn_db}"; then
        msg "ERROR: pre-formatted BLASTn database (.${ext}) for 16S rRNA genes is missing" >&2
        exit 1
      fi
    done

    msg "INFO: Running blastn with !{task.cpus} threads"

    blastn \
      -word_size 10 \
      -task blastn \
      -db 16S_ribosomal_RNA \
      -num_threads "!{task.cpus}" \
      -query "!{extracted_base}" \
      -out "!{prefix}.blast.tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ssciname"

    if verify_minimum_file_size "!{prefix}.blast.tsv" '16S BLASTn Output File' "!{params.min_filesize_blastn_output}"; then
      echo -e "!{prefix}\t16S BLASTn Output File\tPASS" > !{prefix}.16S_BLASTn_Output_File.tsv
    else
      echo -e "!{prefix}\t16S BLASTn Output File\tFAIL" > !{prefix}.16S_BLASTn_Output_File.tsv
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      blast: $(blastn -version | head -n 1 | awk '{print $2}')
    END_VERSIONS
    '''
}
