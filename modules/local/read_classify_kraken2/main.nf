process READ_CLASSIFY_KRAKEN_TWO {

    label "process_high"
    tag { "${meta.id}" }
    container "gregorysprenger/kraken2@sha256:213e70b0f465464b2e52f9f128dcb0cc6761705f6e99b7ce48a5a27a6851083a"

    input:
    tuple val(meta), path(cleaned_fastq_files)
    path database

    output:
    path(".command.{out,err}")
    path("${meta.id}.kraken2_output.tab.gz")
    path("${meta.id}.kraken2_summary.tsv")
    path("versions.yml")                    , emit: versions

    shell:
    '''
    source bash_functions.sh
    source summarize_kraken.sh

    # Investigate taxonomic identity of cleaned reads
    if [ ! -s !{meta.id}.kraken2_summary.tsv ]; then
      msg "INFO: Performing Kraken2 classifications"
      kraken2 \
        --use-names \
        --gzip-compressed \
        --db !{database} \
        --output /dev/null \
        --report kraken2.tab \
        --threads !{task.cpus} \
        !{cleaned_fastq_files[0]} !{cleaned_fastq_files[1]} !{cleaned_fastq_files[2]}

      msg "INFO: Summarizing Kraken2"
      summarize_kraken 'kraken2.tab' > "!{meta.id}.kraken2_summary.tsv"

      # Add header to kraken summary
      sed -i \
        '1i % Reads\t# Reads\tUnclassified\t% Reads\t# Reads\tGenus\t% Reads\t# Reads\tGenus\t% Reads\t# Reads\tSpecies\t% Reads\t# Reads\tSpecies\t% Reads\t# Reads' \
        !{meta.id}.kraken2_summary.tsv

      mv kraken2.tab !{meta.id}.kraken2_output.tab
      gzip !{meta.id}.kraken2_output.tab
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        kraken2: $(kraken2 --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
