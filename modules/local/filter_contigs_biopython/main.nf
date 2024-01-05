process FILTER_CONTIGS_BIOPYTHON {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.uncorrected.fna"), emit: uncorrected_contigs
    path(".command.{out,err}")
    path("versions.yml")                                                 , emit: versions

    shell:
    gcskew = params.filter_contigs_gcskew ? "" : "-g"
    keep_low_complexity = params.filter_contigs_keep_low_complexity ? "" : "-m"
    no_sort = params.filter_contigs_no_sort ? "--no-sort" : ""

    if (params.filter_contigs_discard_file) {
      discard_file = "-d ${params.filter_contigs_discard_file}"
    } else {
      discard_file = ""
    }
    '''
    source bash_functions.sh

    # Remove junk contigs
    filter.contigs.py \
      -i !{contigs} \
      -b "!{meta.id}-!{meta.assembler}" \
      -o "!{meta.id}-!{meta.assembler}.uncorrected.fna" \
      -l !{params.filter_contigs_length} \
      -c !{params.filter_contigs_coverage} \
      --deflines !{params.filter_contigs_deflines} \
      !{no_sort} \
      !{gcskew} \
      !{discard_file} \
      !{keep_low_complexity}

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
