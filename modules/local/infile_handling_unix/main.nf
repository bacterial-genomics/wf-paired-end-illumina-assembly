process INFILE_HANDLING_UNIX {

    tag { "${meta.id}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.Raw_Initial_FastQ_File.tsv"), emit: qc_filecheck
    tuple val(meta), path(reads)                                  , emit: input
    path(".command.{out,err}")
    path("versions.yml")                                          , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Read 1: !{reads[0]}"
    msg "INFO: Read 2: !{reads[1]}"

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}.Raw_Initial_FastQ_File.tsv"

    i=1
    for fastq in !{reads}; do
      # Check if input FastQ file is corrupted
      if [[ ${fastq} =~ .gz ]] && [[ $(gunzip -t ${fastq} 2>&1) ]]; then
        msg "ERROR: Input file ${fastq} is corrupted and assembly cannot be performed!"
        exit 1
      elif [[ ${fastq} =~ .fastq ]] && [[ ! $(cat ${fastq} > /dev/null 2>&1) ]]; then
        msg "ERROR: Input file ${fastq} is corrupted and assembly cannot be performed!"
        exit 1
      fi

      # Check if input FastQ file meets minimum file size requirement
      if verify_minimum_file_size "${fastq}" 'Raw Initial FastQ Files' "!{params.min_filesize_fastq_input}"; then
        echo -e "!{meta.id}\tRaw Initial FastQ (R${i}) File\tPASS" >> "!{meta.id}.Raw_Initial_FastQ_File.tsv"
      else
        echo -e "!{meta.id}\tRaw Initial FastQ (R${i}) File\tFAIL" >> "!{meta.id}.Raw_Initial_FastQ_File.tsv"
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
