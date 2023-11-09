process READ_CLASSIFY_KRAKEN_ONE {

    label "process_high"
    label "process_high_memory"
    tag { "${meta.id}" }
    container "gregorysprenger/kraken@sha256:650ce8ce4a5e313dfafa1726168bb4f7942e543075743766afe1f21ae19abf9c"

    input:
    tuple val(meta), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck)
    path database

    output:
    path ".command.out"
    path ".command.err"
    path "${meta.id}.kraken_output.tab.gz"
    path "${meta.id}.Summary.tsv"
    path "versions.yml"                  , emit: versions

    shell:
    '''
    source bash_functions.sh
    source summarize_kraken.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_nonoverlap_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Investigate taxonomic identity of cleaned reads
    if [ ! -s !{meta.id}.Summary.tsv ]; then
      msg "INFO: Performing Kraken1 classifications"
      kraken \
        --fastq-input \
        --db !{database} \
        --gzip-compressed \
        --threads !{task.cpus} \
        !{paired_R1_gz} !{paired_R2_gz} !{single_gz} \
        > !{meta.id}_kraken.output

      msg "INFO: Creating Kraken Report"
      kraken-report \
        --db !{database} \
        !{meta.id}_kraken.output \
        > kraken.tab 2>&1 | tr '^M' '\n' 1>&2

      msg "INFO: Summarizing Kraken1"
      summarize_kraken 'kraken.tab' > !{meta.id}.Summary.tsv

      # Add header to kraken summary
      sed -i \
        '1i % Reads\t# Reads\tUnclassified\t% Reads\t# Reads\tGenus\t% Reads\t# Reads\tGenus\t% Reads\t# Reads\tSpecies\t% Reads\t# Reads\tSpecies\t% Reads\t# Reads' \
        !{meta.id}.Summary.tsv

      mv kraken.tab !{meta.id}.kraken_output.tab
      gzip !{meta.id}.kraken_output.tab
    fi

    # Get process version information
    cat <<-END_VERSIONS | sed -r 's/^ {4}//' | sed "s/\bEND_VERSIONS\b//" > versions.yml
    "!{task.process}":
        kraken: $(kraken --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
