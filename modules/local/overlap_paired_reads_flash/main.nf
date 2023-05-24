process OVERLAP_PAIRED_READS_FLASH {

    publishDir "${params.outdir}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{overlap.tsv,clean-reads.tsv,gz}"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Non-overlapping_FastQ_Files.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    label "process_low"
    tag { "${prefix}" }

    container "snads/flash@sha256:363b2f44d040c669191efbc3d3ba99caf5efd3fdef370af8f00f3328932143a6"

    input:
    tuple val(prefix), path(reads), path(qc_input_filecheck), path(paired_R1), path(paired_R2), path(qc_adapter_filecheck)

    output:
    tuple val(prefix), path("${prefix}_R1.paired.fq.gz"), path("${prefix}_R2.paired.fq.gz"), path("${prefix}.single.fq.gz"), path("*File*.tsv"), emit: gzip_reads
    path "${prefix}.Non-overlapping_FastQ_Files.tsv", emit: qc_nonoverlap_filecheck
    path "${prefix}.overlap.tsv"
    path "${prefix}.clean-reads.tsv"
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_input_filecheck} !{qc_adapter_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Determine read length based on the first 100 reads
    echo -e "$(zcat "!{reads[0]}" | head -n 400 > read_R1_len.txt)"
    READ_LEN=$(awk 'NR%4==2 {if(length > x) {x=length; y=$0}} END{print length(y)}' read_R1_len.txt)

    OVERLAP_LEN=$(echo | awk -v n=${READ_LEN} '{print int(n*0.8)}')
    msg "INFO: ${READ_LEN} bp read length detected from raw input" >&2

    # Merge overlapping sister reads into singleton reads
    if [ ${OVERLAP_LEN} -gt 0 ]; then
      msg "INFO: ${OVERLAP_LEN} bp overlap will be required for sister reads to be merged" >&2

      msg "INFO: Running flash with !{task.cpus} threads"
      flash \
        -m ${OVERLAP_LEN} \
        -M ${READ_LEN} \
        -o flash \
        -t !{task.cpus} \
        !{paired_R1} !{paired_R2}

      for suff in notCombined_1.fastq notCombined_2.fastq; do
        if verify_minimum_file_size "flash.${suff}" 'Non-overlapping FastQ Files' "!{params.min_filesize_non_overlapping_fastq}"; then
          echo -e "!{prefix}\tNon-overlapping FastQ File (${suff})\tPASS" \
            >> !{prefix}.Non-overlapping_FastQ_Files.tsv
        else
          echo -e "!{prefix}\tNon-overlapping FastQ File (${suff})\tFAIL" \
            >> !{prefix}.Non-overlapping_FastQ_Files.tsv
        fi
      done

      rm !{paired_R1} !{paired_R2}
      mv flash.notCombined_1.fastq !{prefix}_R1.paired.fq
      mv flash.notCombined_2.fastq !{prefix}_R2.paired.fq

      if [ -f flash.extendedFrags.fastq ] && \
      [ -s flash.extendedFrags.fastq ]; then
        CNT_READS_OVERLAPPED=$(awk '{lines++} END{print lines/4}' \
        flash.extendedFrags.fastq)

        cat flash.extendedFrags.fastq >> !{prefix}.single.fq
        rm flash.extendedFrags.fastq
      fi

      msg "INFO: ${CNT_READS_OVERLAPPED:-0} pairs overlapped into singleton reads" >&2
      echo -e "!{prefix}\t${CNT_READS_OVERLAPPED:-0} reads Overlapped" \
        > !{prefix}.overlap.tsv
    fi

    # Summarize final read set and compress
    count_R1=$(echo $(cat !{prefix}_R1.paired.fq | wc -l))
    CNT_CLEANED_PAIRS=$(echo $((${count_R1}/4)))
    msg "INFO: CNT_CLEANED_PAIRS ${CNT_CLEANED_PAIRS}"

    count_single=$(echo $(cat !{prefix}.single.fq | wc -l))
    CNT_CLEANED_SINGLETON=$(echo $((${count_single}/4)))
    msg "INFO: CNT_CLEANED_SINGLETON ${CNT_CLEANED_SINGLETON}"

    echo -e "!{prefix}\t${CNT_CLEANED_PAIRS} cleaned pairs\t${CNT_CLEANED_SINGLETON} cleaned singletons" \
      > !{prefix}.clean-reads.tsv

    gzip !{prefix}.single.fq \
    !{prefix}_R1.paired.fq \
    !{prefix}_R2.paired.fq

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      flash: $(flash --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}