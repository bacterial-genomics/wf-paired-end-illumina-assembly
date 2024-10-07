process MLST_MLST {

    tag { "${meta.id}-${meta.assembler}" }
    container "staphb/mlst@sha256:17e78a25fc5171706b22c8c3d4b1ca2352593b56fef8f28401dd5da3e2e7abe8"  // staphb/mlst:2.23.0-2024-09-01

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

    msg "INFO: Looking for MLST schemes to exclude ..."

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

    msg "INFO: Excluding MLST schemes: ${exclude_list}"

    if [[ -s !{assembly} ]]; then
      msg "INFO: Performing MLST ..."

      mlst \
        "!{assembly}" \
        !{min_score} \
        !{min_identity} \
        !{min_coverage} \
        --novel "!{meta.id}-!{meta.assembler}.MLST.novel.fasta" \
        --threads !{task.cpus} \
        --scheme "${mlst_scheme}" \
        --exclude "${exclude_list}" \
        > "!{meta.id}-!{meta.assembler}.MLST.tsv"

      msg "INFO: Completed MLST genotyping"

      # Print header line and add in Sample_name identifier to data row
      awk -F $'\t' -v id="!{meta.id}" \
        'BEGIN{
          OFS=FS
          print "Sample_name" OFS "PubMLST_scheme_name" OFS "Sequence_type_(ST-#)" OFS "Allele_numbers"
        }
        {$1=id; print}' \
        "!{meta.id}-!{meta.assembler}.MLST.tsv" \
        > tmp \
        && \
        mv tmp "!{meta.id}-!{meta.assembler}.MLST.tsv"

      msg "INFO: Appended header to MLST summary output file"

    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        mlst: $(mlst --version | awk '{print $2}')
    END_VERSIONS
    '''
}
