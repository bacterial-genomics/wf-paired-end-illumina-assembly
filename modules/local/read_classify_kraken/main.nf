process READ_CLASSIFY_KRAKEN_ONE {

    publishDir "${params.outdir}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tab*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    label "process_high"
    label "process_high_memory"
    tag { "${prefix}" }

    container "gregorysprenger/kraken@sha256:650ce8ce4a5e313dfafa1726168bb4f7942e543075743766afe1f21ae19abf9c"

    input:
    tuple val(prefix), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck)

    output:
    path "${prefix}.taxonomy1-reads.tab"
    path "${prefix}_kraken1.tab.gz"
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

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

    # If user doesn't provide a non-default db path, use the path in the Docker container,
    #  which contains a smaller minikraken database
    if [[ -d "!{params.kraken1_db}" ]]; then
      database="!{params.kraken1_db}"
      msg "INFO: Using user specified Kraken 1 database: !{params.kraken1_db}"
    else
      database="/kraken-database/"
      msg "INFO: Using pre-loaded MiniKraken database for Kraken 1"
    fi

    # Confirm the db exists
    for ext in idx kdb; do
      if ! verify_minimum_file_size "${database}/database.${ext}" 'kraken database' '10c'; then
        msg "ERROR: pre-formatted kraken database (.${ext}) for read classification is missing" >&2
        exit 1
      fi
    done

    # Investigate taxonomic identity of cleaned reads
    if [ ! -s !{prefix}.taxonomy1-reads.tab ]; then
      msg "INFO: Running Kraken1 with !{task.cpus} threads"
      kraken \
        --db ${database} \
        --threads !{task.cpus} \
        --fastq-input \
        --gzip-compressed \
        !{paired_R1_gz} !{paired_R2_gz} !{single_gz} \
        > !{prefix}_kraken.output

      msg "INFO: Running kraken-report"
      kraken-report \
        --db ${database} \
        !{prefix}_kraken.output \
        > kraken.tab 2>&1 | tr '^M' '\n' 1>&2

      msg "INFO: Summarizing Kraken1"
      summarize_kraken 'kraken.tab' > !{prefix}.taxonomy1-reads.tab

      mv kraken.tab !{prefix}_kraken1.tab
      gzip !{prefix}_kraken1.tab
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      kraken: $(kraken --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}

process READ_CLASSIFY_KRAKEN_TWO {

    publishDir "${params.outdir}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tab*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    label "process_high"
    tag { "${prefix}" }

    container "gregorysprenger/kraken2@sha256:94517863a0aa1fa820275df349777859d5c194523335978ae154e8ea8190c71d"

    input:
    tuple val(prefix), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck)

    output:
    path "${prefix}.taxonomy2-reads.tab"
    path "${prefix}_kraken2.tab.gz"
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

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

    # If user doesn't provide a non-default db path, use the path in the Docker container,
    #  which contains a smaller minikraken database
    if [[ -d "!{params.kraken2_db}" ]]; then
      database="!{params.kraken2_db}"
      msg "INFO: Using user specified Kraken 2 database: !{params.kraken2_db}"
    else
      database="/kraken2-database"
      msg "INFO: Using pre-loaded MiniKraken2 database for Kraken 2"
    fi

    # Confirm the db exists
    for pref in hash opts taxo; do
      if ! verify_minimum_file_size "${database}/${pref}".k2d 'kraken2 database' '10c'; then
        msg "ERROR: pre-formatted kraken2 database (${pref}.k2d) for read classification is missing" >&2
        exit 1
      fi
    done

    # Investigate taxonomic identity of cleaned reads
    if [ ! -s !{prefix}.taxonomy2-reads.tab ]; then
      msg "INFO: Running Kraken2 with !{task.cpus} threads"
      kraken2 \
        --db "${database}" \
        --threads !{task.cpus} \
        --gzip-compressed \
        --output /dev/null \
        --use-names \
        --report kraken2.tab \
        !{paired_R1_gz} !{paired_R2_gz} !{single_gz}

      msg "INFO: Summarizing Kraken2"
      summarize_kraken 'kraken2.tab' > !{prefix}.taxonomy2-reads.tab

      mv kraken2.tab !{prefix}_kraken2.tab
      gzip !{prefix}_kraken2.tab
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      kraken2: $(kraken2 --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
