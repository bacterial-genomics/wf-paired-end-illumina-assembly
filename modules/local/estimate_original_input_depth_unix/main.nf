process ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX {

    publishDir   "${params.process_log_dir}",
        mode:    "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs:  { filename -> "${meta.id}.${task.process}${filename}" }

    tag { "${meta.id}" }

    container "ubuntu:jammy"

    input:
    tuple val(meta), val(total_bp), val(genome_size)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                                                                   , emit: versions
    tuple val(meta), path("${meta.id}.initial_depth.txt"), path("${meta.id}.fraction_of_reads_to_use.txt"), emit: fraction_of_reads_to_use

    shell:
    '''
    source bash_functions.sh

    bp=$(cat !{total_bp})
    size=$(cat !{genome_size})

    initial_depth=$(( ${bp} / ${size} ))
    msg "INFO: Initial input depth of coverage estimated to be ${initial_depth}x"

    # Calculate the fraction of reads to subsample
    fraction_of_reads_to_use=$(awk \
      -v OFMT='%.6f' \
      -v initial_depth="${initial_depth}" \
      -v want_depth="!{params.depth}" \
      'BEGIN {i = want_depth / initial_depth ; print i}')

    if ! [[ ${fraction_of_reads_to_use} =~ ^[0-9.]+$ ]]; then
      msg "ERROR: unable to calculate fraction of reads to use: ${fraction_of_reads_to_use}" >&2
      exit 1
    fi
    msg "INFO: Fraction of reads to use: ${fraction_of_reads_to_use}"

    echo -n "${initial_depth}" > !{meta.id}.initial_depth.txt
    echo -n "${fraction_of_reads_to_use}" > !{meta.id}.fraction_of_reads_to_use.txt

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
