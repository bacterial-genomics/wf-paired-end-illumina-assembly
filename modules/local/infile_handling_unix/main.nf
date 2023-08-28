process INFILE_HANDLING_UNIX {

    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Raw_Initial_FastQ_Files.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${meta.id}.${task.process}${filename}" }

    container "ubuntu:jammy"
    tag { "${meta.id}" }

    input:
    tuple val(meta), path(reads)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                             , emit: versions
    tuple val(meta), path(reads), path("*File*.tsv"), emit: input
    path "${meta.id}.Raw_Initial_FastQ_Files.tsv"   , emit: qc_input_filecheck

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Read 1: !{reads[0]}"
    msg "INFO: Read 2: !{reads[1]}"

    i=1
    for fastq in !{reads}; do
      if verify_minimum_file_size "${fastq}" 'Raw Initial FastQ Files' "!{params.min_filesize_fastq_input}"; then
        echo -e "!{meta.id}\tRaw Initial FastQ (R${i}) File\tPASS" \
        >> !{meta.id}.Raw_Initial_FastQ_Files.tsv
      else
        echo -e "!{meta.id}\tRaw Initial FastQ (R${i}) File\tFAIL" \
        >> !{meta.id}.Raw_Initial_FastQ_Files.tsv
      fi
      ((i++))
    done

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
