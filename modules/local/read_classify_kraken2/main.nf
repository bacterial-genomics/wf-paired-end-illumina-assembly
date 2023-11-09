process READ_CLASSIFY_KRAKEN_TWO {

    label "process_high"
    tag { "${meta.id}" }
    container "gregorysprenger/kraken2@sha256:e4282b158c9899382a23b53b8e07a791b83d8fd89e502b640a8cd0a411f6ca72"

    input:
    tuple val(meta), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck)
    path database

    output:
    path ".command.out"
    path ".command.err"
    path "${meta.id}_kraken2.tab.gz"
    path "${meta.id}.taxonomy2-reads.tab"
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
    if [ ! -s !{meta.id}.taxonomy2-reads.tab ]; then
      msg "INFO: Performing Kraken2 classifications"
      kraken2 \
        --use-names \
        --gzip-compressed \
        --db !{database} \
        --output /dev/null \
        --report kraken2.tab \
        --threads !{task.cpus} \
        !{paired_R1_gz} !{paired_R2_gz} !{single_gz}

      msg "INFO: Summarizing Kraken2"
      summarize_kraken 'kraken2.tab' > !{meta.id}.taxonomy2-reads.tab

      mv kraken2.tab !{meta.id}_kraken2.tab
      gzip !{meta.id}_kraken2.tab
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        kraken2: $(kraken2 --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
