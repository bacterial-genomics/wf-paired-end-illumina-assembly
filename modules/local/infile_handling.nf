process INFILE_HANDLING {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${basename}.${task.process}${filename}" }

    container "ubuntu:focal"

    input:
        tuple val(basename), path(input)

    output:
        path input, emit: input
        val basename, emit: base
        val size, emit: size
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions
        
    shell:
        if (params.size) {
            size=params.size
        }
        else {
            size=input[0].size();
        }
        '''
        source bash_functions.sh
        
        msg "INFO: R1 = !{input[0]}"
        msg "INFO: R2 = !{input[1]}"

        for fastq in !{input}; do
            verify_file_minimum_size ${fastq} 'fastq' '10M'
        done

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
        END_VERSIONS
        '''
}
