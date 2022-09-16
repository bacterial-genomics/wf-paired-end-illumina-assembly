process INFILE_HANDLING {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    input:
        tuple val(basename), path(input)

    output:
        path input, emit: input
        val basename, emit: base
        val size, emit: size
        path ".command.out"
        path ".command.err"
        
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

    verify_file_minimum_size !{input[0]} 'fastq' '10M'
    verify_file_minimum_size !{input[1]} 'fastq' '10M'

    '''
}