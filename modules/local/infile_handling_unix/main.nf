process INFILE_HANDLING_UNIX {

    // errorStrategy 'terminate'

    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Raw_Initial_FastQ_Files.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}" }

    container "ubuntu:jammy"
    tag { "${base}" }

    input:
        tuple val(basename), path(input)

    output:
        tuple val(base), path(input), path("*File*.tsv"), emit: input
        path "${base}.Raw_Initial_FastQ_Files.tsv", emit: qc_input_filecheck
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions
        
    shell:
        // Split basename on first underscore if applicable
        base=basename.split('_')[0];
        '''
        source bash_functions.sh
        
        msg "INFO: R1 = !{input[0]}"
        msg "INFO: R2 = !{input[1]}"

        i=1
        for fastq in !{input}; do
          if verify_minimum_file_size "${fastq}" 'Raw Initial FastQ Files' "!{params.min_filesize_fastq_input}"; then
            echo -e "!{base}\tRaw Initial FastQ (R${i}) File\tPASS" \
             >> !{base}.Raw_Initial_FastQ_Files.tsv
            touch !{base}.filecheck.txt
          else
            echo -e "!{base}\tRaw Initial FastQ (R${i}) File\tFAIL" \
             >> !{base}.Raw_Initial_FastQ_Files.tsv
          fi
          ((i++))
        done

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
        END_VERSIONS
        '''
}