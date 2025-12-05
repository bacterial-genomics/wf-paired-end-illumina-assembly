process ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX {

    tag { "${meta.id}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), val(input_total_bp_file), val(genome_size_file)

    output:
    tuple val(meta), path("${meta.id}.Estimated_Initial_Depth.tsv"), path("${meta.id}.Read_Fraction_to_Subsample.tsv"), emit: fraction_of_reads_to_use_file
    path(".command.{out,err}")
    path("versions.yml")                                                                                  , emit: versions

    script:
    """
    source bash_functions.sh

    bp=\$(awk 'NR==2 {print (\$2 != "" ? \$2 : 0)}' ${input_total_bp_file})
    size=\$(awk 'NR==2 {print (\$2 != "" ? \$2 : 0)}' ${genome_size_file})

    initial_depth=\$(( \${bp} / \${size} ))
    msg "INFO: Estimated coverage depth of ${meta.id}: \${initial_depth}x"

    # Calculate the fraction of reads to subsample
    read_fraction_to_use=\$(awk \
      -v OFMT='%.6f' \
      -v initial_depth="\${initial_depth}" \
      -v want_depth="${params.depth}" \
      'BEGIN {i = want_depth / initial_depth ; print i}')

    if ! [[ \${read_fraction_to_use} =~ ^[0-9.]+\$ ]]; then
      msg "ERROR: unable to calculate fraction of reads to use: \${read_fraction_to_use}" >&2
      exit 1
    fi
    msg "INFO: Fraction of reads to use: \${read_fraction_to_use}"

    # Form the summary TSV outputs
    echo -e "Sample_name\tEstimated_Original_Depth_[x]" > "${meta.id}.Estimated_Initial_Depth.tsv"
    echo -e "${meta.id}\t\${initial_depth}" >> "${meta.id}.Estimated_Initial_Depth.tsv"

    echo -e "Sample_name\tRead_Fraction_to_Subsample" > "${meta.id}.Read_Fraction_to_Subsample.tsv"
    echo -e "${meta.id}\t\${read_fraction_to_use}" >> "${meta.id}.Read_Fraction_to_Subsample.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(awk -F ' ' '{print \$2,\$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    """
}
