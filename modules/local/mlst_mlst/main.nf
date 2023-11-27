process MLST_MLST {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/mlst@sha256:69c8c8027474b8f361ef4a579df171702f3ed52f45e3fb388a41ccbf4542706f"

    input:
    tuple val(meta), path(assembly)

    output:
    path(".command.{out,err}")
    path("versions.yml")                                 , emit: versions
    path("${meta.id}-${meta.assembler}.Summary.MLST.tab"), emit: summary_mlst

    shell:
    '''
    source bash_functions.sh

    # MLST for each assembly
    msg "INFO: Performing MLST"

    if [[ -s !{assembly} ]]; then
      mlst \
        "!{assembly}" \
        --threads !{task.cpus} \
        >> "!{meta.id}-!{meta.assembler}.Summary.MLST.tab"

      sed -i \
        '1i Filename\tPubMLST scheme name\tSequence type\tAllele IDs' \
        "!{meta.id}-!{meta.assembler}.Summary.MLST.tab"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        mlst: $(mlst --version | awk '{print $2}')
    END_VERSIONS
    '''
}
