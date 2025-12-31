process SUMMARIZE_SUBSAMPLING {

    label "process_low"
    tag { "${meta.id}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads)            , emit: reads
    path("${meta.id}.Downsample_Status.tsv"), emit: output
    path(".command.{out,err}")
    path("versions.yml")                    , emit: versions

    script:
    """
    echo -e "Sample_name\tDownsampled_[Yes|No]" > "${meta.id}.Downsample_Status.tsv"
    if [[ "${meta.downsampled}" == "true" ]]; then
      echo -e "${meta.id}\tYes" >> "${meta.id}.Downsample_Status.tsv"
    else
      echo -e "${meta.id}\tNo" >> "${meta.id}.Downsample_Status.tsv"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(awk -F ' ' '{print \$2,\$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    """
}
