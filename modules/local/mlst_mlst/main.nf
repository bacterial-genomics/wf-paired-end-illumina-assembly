process MLST_MLST {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/mlst@sha256:69c8c8027474b8f361ef4a579df171702f3ed52f45e3fb388a41ccbf4542706f"

    input:
    tuple val(meta), path(assembly)

    output:
    path("${meta.id}-${meta.assembler}.MLST.tsv"), emit: summary
    path(".command.{out,err}")
    path("versions.yml")                         , emit: versions

    shell:
    scheme       = params.mlst_scheme        ? "${params.mlst_scheme.toLowerCase()}"       : ''
    exclude      = params.mlst_ignore_scheme ? "${params.mlst_ignore_scheme.toLowerCase()}": ''
    min_score    = params.mlst_min_score     ? "--minscore ${params.mlst_min_score}"       : "--minscore '50'"
    min_identity = params.mlst_min_identity  ? "--minid ${params.mlst_min_identity}"       : "--minid '95'"
    min_coverage = params.mlst_min_coverage  ? "--mincov ${params.mlst_min_coverage}"      : "--mincov '10'"
    '''
    source bash_functions.sh

    # MLST for each assembly
    msg "INFO: Performing MLST"

    # Check if input scheme is in mlst's database
    mlst_scheme="!{scheme}"
    if [[ "!{scheme}" != '' ]] && \
      [[ ! $(mlst --list 2>&1 | tail -n 1 | grep -w "${mlst_scheme}") ]]; then
      msg "WARN: Specified MLST scheme is not valid. Defaulting to auto detecting the scheme."
      mlst_scheme=''
    fi

    # Check if scheme to ignore is in mlst's database
    exclude_list=()
    for e in $(echo "!{exclude}" | tr ',' ' '); do
      if [[ "${e}" != '' ]] && \
        [[ $(mlst --list 2>&1 | tail -n 1 | grep -w "${e}") ]]; then
        exclude_list+=( "${e}" )
      fi
    done

    # Reformat exclude list
    if [[ -z ${exclude_list[@]} ]]; then
      exclude_list=''
    else
      exclude_list=$(echo ${exclude_list[@]} | tr ' ' ',')
    fi

    if [[ -s !{assembly} ]]; then
      mlst \
        "!{assembly}" \
        !{min_score} \
        !{min_identity} \
        !{min_coverage} \
        --threads !{task.cpus} \
        --scheme "${mlst_scheme}" \
        --exclude "${exclude_list}" \
        >> "!{meta.id}-!{meta.assembler}.MLST.tsv"

      # Replace the assembly file (with "-<assembler pkg>" in first column with nextflow's
      #   sample_name in the report, to ensure Sample_name consistent in report.
      # NOTE: Sample name alone not clear enough which assembler was used at this stage to
      #       avoid confusion when > 1 assembler is used in same output directory path.
      awk -v id="!{meta.id}-!{meta.assembler}" \
        'BEGIN{FS=OFS="\t"} {$1=id; print}' \
        "${meta.id}-${meta.assembler}.MLST.tsv" \
        > tmp \
      && mv tmp "${meta.id}-${meta.assembler}.MLST.tsv"

      # Add header line to data output
      # NOTE: use awk instead of sed to avoid error on () characters
      awk \
        'BEGIN{print "Sample_name-Assembler\tPubMLST_scheme_name\tSequence_type_(ST-#)\tAllele_numbers"}
         {print}
        ' \
        "${meta.id}-${meta.assembler}.MLST.tsv" \
        > tmp \
      && \
      mv tmp "${meta.id}-${meta.assembler}.MLST.tsv"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        mlst: $(mlst --version | awk '{print $2}')
    END_VERSIONS
    '''
}
