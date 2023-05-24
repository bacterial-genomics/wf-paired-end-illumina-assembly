process INFILE_HANDLING_UNIX {

    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Raw_Initial_FastQ_Files.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}" }

    container "ubuntu:jammy"
    tag { "${prefix}" }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(prefix), path(reads), path("*File*.tsv"), emit: input
    path "${prefix}.Raw_Initial_FastQ_Files.tsv", emit: qc_input_filecheck
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions
        
    shell:
    // Split meta.id on first underscore if applicable
    prefix="${meta.id}".split('_')[0];
    '''
    source bash_functions.sh
    
    msg "INFO: R1 = !{reads[0]}"
    msg "INFO: R2 = !{reads[1]}"

    i=1
    for fastq in !{reads}; do
      if verify_minimum_file_size "${fastq}" 'Raw Initial FastQ Files' "!{params.min_filesize_fastq_input}"; then
        echo -e "!{prefix}\tRaw Initial FastQ (R${i}) File\tPASS" \
        >> !{prefix}.Raw_Initial_FastQ_Files.tsv
      else
        echo -e "!{prefix}\tRaw Initial FastQ (R${i}) File\tFAIL" \
        >> !{prefix}.Raw_Initial_FastQ_Files.tsv
      fi
      ((i++))
    done

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}