process READ_CLASSIFY_KRAKEN_ONE {

    label "process_high"
    label "process_high_memory"
    tag { "${meta.id}" }
    container "gregorysprenger/kraken@sha256:b5ab4b75fb197b16e81d8cc3878e08479bc7d105ac0b2e948e6f6a9985cfc93e"

    input:
    tuple val(meta), path(cleaned_fastq_files)
    path database

    output:
    path(".command.{out,err}")
    path("${meta.id}.kraken_output.tab.gz")
    path("${meta.id}.kraken_summary.tsv")
    path("versions.yml")                   , emit: versions

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
        > kraken.tab 2>&1 | tr '^M' '\n' 1>&2

      msg "INFO: Summarizing Kraken1"
      summarize_kraken 'kraken.tab' > "!{meta.id}.kraken_summary.tsv"

      # Add header to kraken summary
      sed -i \
        '1i % Reads\t# Reads\tUnclassified\t% Reads\t# Reads\tGenus\t% Reads\t# Reads\tGenus\t% Reads\t# Reads\tSpecies\t% Reads\t# Reads\tSpecies\t% Reads\t# Reads' \
        !{meta.id}.kraken_summary.tsv

      mv kraken.tab !{meta.id}.kraken_output.tab
      gzip !{meta.id}.kraken_output.tab
    fi

    # Get process version information
    cat <<-"    END_VERSIONS" | sed -r 's/^ {4}//' > versions.yml
    "!{task.process}":
        kraken: $(kraken --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
