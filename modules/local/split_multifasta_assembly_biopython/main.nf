process SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON {

    publishDir   "${params.process_log_dir}",
        mode:    "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs:  { filename -> "${meta.id}.${task.process}${filename}" }

    tag { "${meta.id}" }

    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(assembly)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"            , emit: versions
    tuple val(meta), path("bins/*"), emit: split_multifasta_assembly_dir

    shell:
    no_gaps = params.split_multifasta_remove_gaps ? "--nogaps" : ""
    '''
    source bash_functions.sh

    # Split assembly multi-record FastA into individual FastA files for each contig
    if [[ -s "!{assembly}" ]]; then
      split.multifasta.py \
        --outdir "bins" \
        --infile "!{assembly}" \
        !{no_gaps} \
        --suffix '' \
        --ext !{params.split_multifasta_extension}

    else
      msg "ERROR: ${assembly} absent" >&2
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
