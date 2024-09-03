process READ_CLASSIFY_KRAKEN_ONE {

    label "process_high"
    tag { "${meta.id}" }
    container "gregorysprenger/kraken@sha256:b5ab4b75fb197b16e81d8cc3878e08479bc7d105ac0b2e948e6f6a9985cfc93e"

    input:
    tuple val(meta), path(cleaned_fastq_files)
    path database

    output:
    path("*.kraken_output.tsv.gz")
    path("*.kraken_summary.tsv")  , emit: summary
    path(".command.{out,err}")
    path("versions.yml")          , emit: versions

    shell:
    '''
    source bash_functions.sh
    source summarize_kraken.sh

    # Investigate taxonomic identity of cleaned reads
    if [ ! -s !{meta.id}.kraken_summary.tsv ]; then
      msg "INFO: Performing Kraken1 classifications"
      kraken \
        --fastq-input \
        --db !{database} \
        --gzip-compressed \
        --threads !{task.cpus} \
        !{cleaned_fastq_files[0]} !{cleaned_fastq_files[1]} !{cleaned_fastq_files[2]} \
        > "!{meta.id}_kraken.output"

      msg "INFO: Creating Kraken Report"
      kraken-report \
        --db !{database} \
        !{meta.id}_kraken.output \
        > kraken.tsv 2>&1 | tr '^M' '\\n' 1>&2

      echo -ne "!{meta.id}\t" > "!{meta.id}.kraken_summary.tsv"
      summarize_kraken 'kraken.tsv' | sed 's/%//g' >> "!{meta.id}.kraken_summary.tsv"

      # Add header to output
      SUMMARY_HEADER=(
        "Sample_name"
        "Reads_(%)" "Reads_(#)" "Unclassified"
        "Reads_(%)" "Reads_(#)" "Genus" "Reads_(%)" "Reads_(#)" "Genus" "Reads_(%)" "Reads_(#)" "Genus"
        "Reads_(%)" "Reads_(#)" "Species" "Reads_(%)" "Reads_(#)" "Species" "Reads_(%)" "Reads_(#)" "Species"
      )
      SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')
      sed -i "1i ${SUMMARY_HEADER}" "!{meta.id}.kraken_summary.tsv"

      mv kraken.tsv "!{meta.id}.kraken_output.tsv"
      gzip "!{meta.id}.kraken_output.tsv"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        kraken: $(kraken --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
