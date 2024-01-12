process MLST_MLST {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/mlst@sha256:69c8c8027474b8f361ef4a579df171702f3ed52f45e3fb388a41ccbf4542706f"

    input:
    tuple val(meta), path(assembly)

    output:
    path("${meta.id}-${meta.assembler}.Summary.MLST.tab"), emit: summary_mlst
    path(".command.{out,err}")
    path("versions.yml")                                 , emit: versions

    shell:
    scheme = params.mlst_scheme ? params.mlst_scheme : ''
    exclude = params.mlst_ignore_scheme ? params.mlst_ignore_scheme : ''
    '''
    source bash_functions.sh

    # MLST for each assembly
    msg "INFO: Performing MLST"

    # Check if input scheme is in mlst's database
    mlst_scheme=!{scheme}
    if [[ !{scheme} != '' ]] && \
      [[ ! $(mlst --list 2>&1 | tail -n 1 | grep -w $scheme) ]]; then
      msg "WARN: Specified MLST scheme is not valid. Defaulting to auto detecting the scheme."
      mlst_scheme=''
    fi

    # Check if scheme to ignore is in mlst's database
    exclude_list=()
    for e in $(echo !{exclude} | tr ',' ' '); do
      if [[ ${e} != '' ]] && \
        [[ $(mlst --list 2>&1 | tail -n 1 | grep -w $e) ]]; then
        exclude_list+=( ${e} )
    done

    # Reformat exclude list
    exclude_list=$(echo ${exclude_list[@]} | tr ' ' ',')

    if [[ -s !{assembly} ]]; then
      mlst \
        "!{assembly}" \
        --threads !{task.cpus} \
        --scheme ${mlst_scheme} \
        --exclude ${exclude_list} \
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
