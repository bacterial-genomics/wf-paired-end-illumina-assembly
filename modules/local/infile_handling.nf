process INFILE_HANDLING {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}" }

    container "ubuntu:focal"
    tag { "${base}" }

    input:
        tuple val(basename), path(input)

    output:
        tuple val(base), val(size), path(input), emit: input
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions
        
    shell:
        // Get filesize from input file or size set by params.size
        if (params.size != "null") {
            size=params.size
        }
        else {
            size=input[0].size();
        }

        // Split basename on first underscore if applicable
        base=basename.split('_')[0];
        '''
        source bash_functions.sh
        
        msg "INFO: R1 = !{input[0]}"
        msg "INFO: R2 = !{input[1]}"

        for fastq in !{input}; do
            verify_file_minimum_size "${fastq}" 'fastq' '10M' '100'
        done

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}(!{basename})":
            ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
        END_VERSIONS
        '''
}
