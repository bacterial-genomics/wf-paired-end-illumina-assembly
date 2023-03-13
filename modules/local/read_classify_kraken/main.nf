process READ_CLASSIFY_KRAKEN_ONE {

    // errorStrategy 'terminate'

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tab*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    label "process_high"
    tag { "${base}" }

    container "staphb/kraken@sha256:d372099288c3a7c0cc90ea7e516c643e7096c90a551b45d531bd26b4e7f46255"

    input:
        tuple val(base), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_overlap_filecheck)
        path kraken1_db

    output:
        path "${base}.taxonomy1-reads.tab"
        path "${base}_kraken1.tab.gz"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        # Exit if previous process fails qc filechecks
        if [ $(grep "FAIL" !{base}*File*.tsv) ]; then
          exit 1
        fi

        source bash_functions.sh
        source summarize_kraken.sh

        # If user doesn't provide a non-default db path, use the path in the Docker container,
        #  which contains a smaller minikraken database
        if [[ -d "!{kraken1_db}" ]]; then
            database="!{kraken1_db}"
            msg "INFO: Using user specified Kraken 1 database: !{params.kraken1_db}"
        else
            database="/kraken-database/minikraken_20171013_4GB"
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
        if [ ! -s !{base}.taxonomy1-reads.tab ]; then
          msg "INFO: Running Kraken1 with !{task.cpus} threads"
          kraken \
          --db ${database} \
          --threads !{task.cpus} \
          --fastq-input \
          --gzip-compressed \
          !{paired_R1_gz} !{paired_R2_gz} !{single_gz} \
          > !{base}_kraken.output

          msg "INFO: Running kraken-report"
          kraken-report \
          --db ${database} \
          !{base}_kraken.output \
          > kraken.tab 2>&1 | tr '^M' '\n' 1>&2

          msg "INFO: Summarizing Kraken1"
          summarize_kraken 'kraken.tab' > !{base}.taxonomy1-reads.tab

          mv kraken.tab !{base}_kraken1.tab
          gzip !{base}_kraken1.tab
        fi

        # Get process version
        echo -e '"!{task.process} (!{base})":' > versions.yml
        echo -e "    kraken: $(kraken --version | head -n 1 | awk 'NF>1{print $NF}')" >> versions.yml
        '''
}

process READ_CLASSIFY_KRAKEN_TWO {

    // errorStrategy 'terminate'

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tab*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    label "process_high"
    tag { "${base}" }

    container "staphb/kraken2@sha256:5b107d0141d6042a6b0ac6a5852990dc541fbff556a85eb0c321a7771200ba56"

    input:
        tuple val(base), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_overlap_filecheck)
        path kraken2_db

    output:
        path "${base}.taxonomy2-reads.tab"
        path "${base}_kraken2.tab.gz"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        # Exit if previous process fails qc filechecks
        if [ $(grep "FAIL" !{base}*File*.tsv) ]; then
          exit 1
        fi

        source bash_functions.sh
        source summarize_kraken.sh

        # If user doesn't provide a non-default db path, use the path in the Docker container,
        #  which contains a smaller minikraken database
        if [[ -d "!{kraken2_db}" ]]; then
          database="!{kraken2_db}"
          msg "INFO: Using user specified Kraken 2 database: !{params.kraken2_db}"
        else
          database="/kraken2-db"
          msg "INFO: Using pre-loaded MiniKraken2 database for Kraken 2"
        fi

        # Confirm the db exists
        for pref in hash opts taxo; do
          if ! verify_minimum_file_size "${database}/${pref}".k2d 'kraken2 database' '10c'; then
            msg "ERROR: pre-formatted kraken2 database (${ext}.k2d) for read classification is missing" >&2
            exit 1
          fi
        done

        # Investigate taxonomic identity of cleaned reads
        if [ ! -s !{base}.taxonomy2-reads.tab ]; then
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
          summarize_kraken 'kraken2.tab' > !{base}.taxonomy2-reads.tab

          mv kraken2.tab !{base}_kraken2.tab
          gzip !{base}_kraken2.tab
        fi

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            kraken2: $(kraken2 --version | head -n 1 | awk 'NF>1{print $NF}')
        END_VERSIONS
        '''
}