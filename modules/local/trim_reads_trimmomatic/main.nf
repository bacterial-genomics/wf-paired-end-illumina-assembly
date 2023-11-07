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
        saveAs: { filename -> "${meta.id}.${task.process}${filename}" }

    label "process_high"
    tag { "${meta.id}" }

    container "snads/trimmomatic@sha256:afbb19fdf540e6bd508b657e8dafffb6b411b5b0bf0e302347889220a0b571f1"

    input:
    tuple val(meta), path(noPhiX_R1), path(noPhiX_R2), path(qc_phix_filecheck)
    path adapter_reference_file

    output:
    path ".command.out"
    path ".command.err"
    path "${meta.id}.single.fq"
    path "${meta.id}.trimmo.tsv"
    path "versions.yml"                                                                                  , emit: versions
    path "${meta.id}.Adapters_FastA_File.tsv"                                                            , emit: qc_adapters_filecheck
    path "${meta.id}.Adapter-removed_FastQ_Files.tsv"                                                    , emit: qc_removed_adapters_filecheck
    tuple val(meta), path("${meta.id}_R1.paired.fq"), path("${meta.id}_R2.paired.fq"), path("*File*.tsv"), emit: trimmo

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

    # Verify adapter reference file size
    if verify_minimum_file_size !{adapter_reference_file} 'Adapters FastA' "!{params.min_filesize_adapters}"; then
      echo -e "!{meta.id}\tAdapters FastA File\tPASS" > !{meta.id}.Adapters_FastA_File.tsv
    else
      echo -e "!{meta.id}\tAdapters FastA File\tFAIL" > !{meta.id}.Adapters_FastA_File.tsv
    fi

    # Adapter clip and quality trim
    msg "INFO: Performing read trimming with Trimmomatic"

    trimmomatic PE \
      -phred33 \
      -threads !{task.cpus} \
      !{noPhiX_R1} !{noPhiX_R2} \
      !{meta.id}_R1.paired.fq !{meta.id}_R1.unpaired.fq \
      !{meta.id}_R2.paired.fq !{meta.id}_R2.unpaired.fq \
      MINLEN:50 \
      LEADING:10 \
      TRAILING:10 \
      SLIDINGWINDOW:6:30 \
      ILLUMINACLIP:!{adapter_reference_file}:2:20:10:8:TRUE

    TRIMMO_DISCARD=$(grep '^Input Read Pairs: ' .command.err \
    | grep ' Dropped: ' | awk '{print $20}')

    msg "INFO: ${TRIMMO_DISCARD} reads are poor quality and were discarded" >&2

    CNT_BROKEN_R1=$(awk '{lines++} END{print lines/4}' !{meta.id}_R1.unpaired.fq)
    CNT_BROKEN_R2=$(awk '{lines++} END{print lines/4}' !{meta.id}_R2.unpaired.fq)

    if [[ -z "${TRIMMO_DISCARD}" || -z "${CNT_BROKEN_R1}" || -z "${CNT_BROKEN_R2}" ]]; then
      msg 'ERROR: unable to parse discarded read counts from trimmomatic log' >&2
      exit 1
    fi

    CNT_BROKEN=$((${CNT_BROKEN_R1} + ${CNT_BROKEN_R2}))

    msg "INFO: $CNT_BROKEN_R1 forward reads lacked a high quality R2 sister read" >&2
    msg "INFO: $CNT_BROKEN_R2 reverse reads lacked a high quality R1 sister read" >&2
    msg "INFO: $CNT_BROKEN total broken read pairs were saved as singletons" >&2

    echo -e "!{meta.id}\t${TRIMMO_DISCARD} reads Discarded\t${CNT_BROKEN} reads Singletons" \
    > !{meta.id}.trimmo.tsv

    cat !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq > !{meta.id}.single.fq

    rm -f !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq

    for suff in R1.paired.fq R2.paired.fq; do
      if verify_minimum_file_size "!{meta.id}_${suff}" 'Adapter-removed FastQ Files' "!{params.min_filesize_fastq_adapters_removed}"; then
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tPASS" \
          >> !{meta.id}.Adapter-removed_FastQ_Files.tsv
      else
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tFAIL" \
          >> !{meta.id}.Adapter-removed_FastQ_Files.tsv
      fi
    done

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        trimmomatic: $(trimmomatic -version)
    END_VERSIONS
    '''
}
