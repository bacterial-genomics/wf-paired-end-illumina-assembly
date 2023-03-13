process OVERLAP_PAIRED_READS_FLASH {

    // errorStrategy 'terminate'

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{overlap.tsv,clean-reads.tsv,gz}"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Non-overlapping_FastQ_Files.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    label "process_low"
    tag { "${base}" }

    container "snads/flash@sha256:363b2f44d040c669191efbc3d3ba99caf5efd3fdef370af8f00f3328932143a6"

    input:
        tuple val(base), path(input), path(qc_input_filecheck), path(paired_R1), path(paired_R2), path(qc_adapter_filechecks)

    output:
        tuple val(base), path("${base}_R1.paired.fq.gz"), path("${base}_R2.paired.fq.gz"), path("${base}.single.fq.gz"), path("*Non*.tsv"), emit: gzip_reads
        path "${base}.Non-overlapping_FastQ_Files.tsv", emit: qc_filecheck
        path "${base}.overlap.tsv"
        path "${base}.clean-reads.tsv"
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

        # Determine read length based on the first 100 reads
        echo -e "$(zcat "!{input[0]}" | head -n 400 > read_R1_len.txt)"
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
              echo -e "!{base}\tNon-overlapping FastQ File (${suff})\tPASS" \
                >> !{base}.Non-overlapping_FastQ_Files.tsv
            else
              echo -e "!{base}\tNon-overlapping FastQ File (${suff})\tFAIL" \
                >> !{base}.Non-overlapping_FastQ_Files.tsv
            fi
          done

          rm !{paired_R1} !{paired_R2}
          mv flash.notCombined_1.fastq !{base}_R1.paired.fq
          mv flash.notCombined_2.fastq !{base}_R2.paired.fq

          if [ -f flash.extendedFrags.fastq ] && \
          [ -s flash.extendedFrags.fastq ]; then
            CNT_READS_OVERLAPPED=$(awk '{lines++} END{print lines/4}' \
            flash.extendedFrags.fastq)

            cat flash.extendedFrags.fastq >> !{base}.single.fq
            rm flash.extendedFrags.fastq
          fi

          msg "INFO: ${CNT_READS_OVERLAPPED:-0} pairs overlapped into singleton reads" >&2
          echo -e "!{base}\t${CNT_READS_OVERLAPPED:-0} reads Overlapped" \
          > !{base}.overlap.tsv
        fi

        # Summarize final read set and compress
        count_R1=$(echo $(cat !{base}_R1.paired.fq | wc -l))
        CNT_CLEANED_PAIRS=$(echo $((${count_R1}/4)))
        msg "INFO: CNT_CLEANED_PAIRS ${CNT_CLEANED_PAIRS}"

        count_single=$(echo $(cat !{base}.single.fq | wc -l))
        CNT_CLEANED_SINGLETON=$(echo $((${count_single}/4)))
        msg "INFO: CNT_CLEANED_SINGLETON ${CNT_CLEANED_SINGLETON}"

        echo -e "!{base}\t${CNT_CLEANED_PAIRS} cleaned pairs\t${CNT_CLEANED_SINGLETON} cleaned singletons" \
        > !{base}.clean-reads.tsv

        gzip !{base}.single.fq \
        !{base}_R1.paired.fq \
        !{base}_R2.paired.fq

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            flash: $(flash --version | head -n 1 | awk 'NF>1{print $NF}')
        END_VERSIONS
        '''
}