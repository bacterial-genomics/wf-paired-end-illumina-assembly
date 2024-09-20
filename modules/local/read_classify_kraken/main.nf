process READ_CLASSIFY_KRAKEN_ONE {

    label "process_high"
    tag { "${meta.id}" }
    container "staphb/kraken@sha256:6f426bbe8ba0b49b6285d773392a94aa79f424ddc50bfb7a00bb52552ea77267"

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
      msg "INFO: Performing Kraken1 classifications of !{cleaned_fastq_files[0]} !{cleaned_fastq_files[1]} !{cleaned_fastq_files[2]} ..."

      kraken \
        --fastq-input \
        --db !{database} \
        --gzip-compressed \
        --threads !{task.cpus} \
        !{cleaned_fastq_files[0]} !{cleaned_fastq_files[1]} !{cleaned_fastq_files[2]} \
        > "!{meta.id}_kraken.output"

      msg "INFO: Created Kraken1 classifications of !{cleaned_fastq_files[0]} !{cleaned_fastq_files[1]} !{cleaned_fastq_files[2]}"

      msg "INFO: Making Kraken1 report from !{meta.id}_kraken.output ..."

      msg "INFO: Creating Kraken1 Report"
      kraken-report \
        --db !{database} \
        !{meta.id}_kraken.output \
        > kraken.tsv

      msg "INFO: Created Kraken1 report kraken.tsv"

      msg "INFO: Summarizing Kraken1 report kraken.tsv ..."

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

      msg "INFO: Created Kraken1 report !{meta.id}.kraken_summary.tsv"

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
