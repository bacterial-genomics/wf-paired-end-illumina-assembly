process FILTER_CONTIGS_BIOPYTHON {

    tag { "${meta.id}" }
    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(R1), path(R2), path(single), path(qc_nonoverlap_filecheck), path(contigs), path(qc_assembly_filecheck)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                , emit: versions
    tuple val(meta), path("${meta.id}.uncorrected.fna"), emit: uncorrected_contigs

    shell:
    gcskew = params.filter_contigs_gcskew ? "" : "-g"
    keep_low_complexity = params.filter_contigs_keep_low_complexity ? "" : "-m"
    no_sort = params.filter_contigs_no_sort ? "--no-sort" : ""

    if (params.filter_contigs_discard_file == "None") {
      discard_file = ""
    } else {
      discard_file = "-d ${params.filter_contigs_discard_file}"
    }
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_nonoverlap_filecheck} !{qc_assembly_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Remove junk contigs
    filter.contigs.py \
      -i !{contigs} \
      -b "!{meta.id}" \
      -o !{meta.id}.uncorrected.fna \
      !{no_sort} \
      !{gcskew} \
      !{discard_file} \
      !{keep_low_complexity} \
      -c !{params.filter_contigs_coverage} \
      --deflines !{params.filter_contigs_deflines}

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
