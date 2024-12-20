process FILTER_CONTIGS_BIOPYTHON {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.uncorrected.fna"), emit: uncorrected_contigs
    path("${meta.id}-${meta.assembler}.discarded-contigs.fa.gz")         , emit: discarded_contigs
    path("${meta.id}-${meta.assembler}.filter-contigs-stats.txt")        , emit: filter_stats
    path(".command.{out,err}")
    path("versions.yml")                                                 , emit: versions

    shell:
    gcskew = params.filter_contigs_gcskew ? "" : "--gcskew"
    keep_low_complexity = params.filter_contigs_keep_low_complexity ? "" : "--complex"
    no_sort = params.filter_contigs_no_sort ? "--no-sort" : ""

    '''
    source bash_functions.sh

    msg "INFO: Filtering contigs from !{contigs} ..."

    # Remove junk contigs
    filter.contigs.py \
      --infile !{contigs} \
      --baseheader "!{meta.id}-!{meta.assembler}" \
      --outfile "!{meta.id}-!{meta.assembler}.uncorrected.fna" \
      --len !{params.filter_contigs_length} \
      --cov !{params.filter_contigs_coverage} \
      --deflines !{params.filter_contigs_deflines} \
      --discarded "!{meta.id}-!{meta.assembler}.discarded-contigs.fa" \
      !{no_sort} \
      !{gcskew} \
      !{keep_low_complexity} \
      2> "!{meta.id}-!{meta.assembler}.filter-contigs-stats.txt"

    msg "INFO: Completed contig filtering for !{meta.id}"

    if [ -s "!{meta.id}-!{meta.assembler}.discarded-contigs.fa" ]; then
      gzip -9f "!{meta.id}-!{meta.assembler}.discarded-contigs.fa"
      msg "INFO: discarded contigs saved as !{meta.id}-!{meta.assembler}.discarded-contigs.fa.gz"
    else
      msg "INFO: no contigs were discarded, therefore not storing empty !{meta.id}-!{meta.assembler}.discarded-contigs.fa.gz file"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
