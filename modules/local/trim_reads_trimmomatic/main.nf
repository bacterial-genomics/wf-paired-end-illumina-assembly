process TRIM_READS_TRIMMOMATIC {

    publishDir "${params.outdir}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.trimmo.tsv"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{Adapters_FastA_File,Adapter-removed_FastQ_Files}.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    label "process_high"
    tag { "${prefix}" }

    container "snads/trimmomatic@sha256:afbb19fdf540e6bd508b657e8dafffb6b411b5b0bf0e302347889220a0b571f1"

    input:
    tuple val(prefix), path(noPhiX_R1), path(noPhiX_R2), path(qc_phix_filecheck)

    output:
    path ".command.out"
    path ".command.err"
    path "${prefix}.single.fq"
    path "${prefix}.trimmo.tsv"
    path "versions.yml"                                                                                  , emit: versions
    path "${prefix}.Adapters_FastA_File.tsv"                                                             , emit: qc_adapters_filecheck
    path "${prefix}.Adapter-removed_FastQ_Files.tsv"                                                     , emit: qc_removed_adapters_filecheck
    tuple val(prefix), path("${prefix}_R1.paired.fq"), path("${prefix}_R2.paired.fq"), path("*File*.tsv"), emit: trimmo

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_phix_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Get Adapters, check if it exists, and verify file size
    if ! check_if_file_exists_allow_seconds !{params.adapter_reference} '60'; then
      exit 1
    fi
    if verify_minimum_file_size !{params.adapter_reference} 'Adapters FastA' "!{params.min_filesize_adapters}"; then
      echo -e "!{prefix}\tAdapters FastA File\tPASS" > !{prefix}.Adapters_FastA_File.tsv
    else
      echo -e "!{prefix}\tAdapters FastA File\tFAIL" > !{prefix}.Adapters_FastA_File.tsv
    fi

    # Adapter clip and quality trim
    msg "INFO: Running trimmomatic with !{task.cpus} threads"

    trimmomatic PE \
      -phred33 \
      -threads !{task.cpus} \
      !{noPhiX_R1} !{noPhiX_R2} \
      !{prefix}_R1.paired.fq !{prefix}_R1.unpaired.fq \
      !{prefix}_R2.paired.fq !{prefix}_R2.unpaired.fq \
      MINLEN:50 \
      LEADING:10 \
      TRAILING:10 \
      SLIDINGWINDOW:6:30 \
      ILLUMINACLIP:!{params.adapter_reference}:2:20:10:8:TRUE

    TRIMMO_DISCARD=$(grep '^Input Read Pairs: ' .command.err \
    | grep ' Dropped: ' | awk '{print $20}')

    msg "INFO: ${TRIMMO_DISCARD} reads are poor quality and were discarded" >&2

    CNT_BROKEN_R1=$(awk '{lines++} END{print lines/4}' !{prefix}_R1.unpaired.fq)
    CNT_BROKEN_R2=$(awk '{lines++} END{print lines/4}' !{prefix}_R2.unpaired.fq)

    if [[ -z "${TRIMMO_DISCARD}" || -z "${CNT_BROKEN_R1}" || -z "${CNT_BROKEN_R2}" ]]; then
      msg 'ERROR: unable to parse discarded read counts from trimmomatic log' >&2
      exit 1
    fi

    CNT_BROKEN=$((${CNT_BROKEN_R1} + ${CNT_BROKEN_R2}))

    msg "INFO: $CNT_BROKEN_R1 forward reads lacked a high quality R2 sister read" >&2
    msg "INFO: $CNT_BROKEN_R2 reverse reads lacked a high quality R1 sister read" >&2
    msg "INFO: $CNT_BROKEN total broken read pairs were saved as singletons" >&2

    echo -e "!{prefix}\t${TRIMMO_DISCARD} reads Discarded\t${CNT_BROKEN} reads Singletons" \
    > !{prefix}.trimmo.tsv

    cat !{prefix}_R1.unpaired.fq !{prefix}_R2.unpaired.fq > !{prefix}.single.fq

    rm -f !{prefix}_R1.unpaired.fq !{prefix}_R2.unpaired.fq

    for suff in R1.paired.fq R2.paired.fq; do
      if verify_minimum_file_size "!{prefix}_${suff}" 'Adapter-removed FastQ Files' "!{params.min_filesize_fastq_adapters_removed}"; then
        echo -e "!{prefix}\tAdapter-removed ($suff) FastQ File\tPASS" \
          >> !{prefix}.Adapter-removed_FastQ_Files.tsv
      else
        echo -e "!{prefix}\tAdapter-removed ($suff) FastQ File\tFAIL" \
          >> !{prefix}.Adapter-removed_FastQ_Files.tsv
      fi
    done

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      trimmomatic: $(trimmomatic -version)
    END_VERSIONS
    '''
}
