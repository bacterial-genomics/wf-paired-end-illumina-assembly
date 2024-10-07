process OVERLAP_PAIRED_READS_FLASH {

    label "process_high"  //extra RAM unnecessary but this is a bottleneck for CPU speed in the workflow
    tag { "${meta.id}" }
    container "staphb/flash@sha256:44889120b49d3f8eefdde8f6040b096d5ee122ceb71d936b596757e4fc16a2c0"

    input:
    tuple val(meta), path(fastq_pairs)

    output:
    tuple val(meta), path("${meta.id}.Non-overlapping_FastQ_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}*{paired,single}.fq.gz")         , emit: cleaned_fastq_files
    path("${meta.id}.Overlap.tsv")                                    , emit: summary
    path("${meta.id}.Clean_Reads_FastQ.SHA512-checksums.tsv")         , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                              , emit: versions

    shell:
    '''
    source bash_functions.sh

    echo !{fastq_pairs}

    # Determine read length based on the first 100 reads
    echo "$(cat !{meta.id}_R1.paired.fq | head -n 400)" > read_R1_len.txt
    READ_LEN=$(awk 'NR%4==2 {if(length > x) {x=length; y=$0}} END{print length(y)}' read_R1_len.txt)

    # Require 80% overlap length relative to input read length
    OVERLAP_LEN=$(echo | awk -v n=${READ_LEN} '{print int(n*0.8)}')
    msg "INFO: ${READ_LEN} bp read length detected from raw input"

    # Merge overlapping sister reads into singleton reads
    if [ ${OVERLAP_LEN} -gt 0 ]; then
      msg "INFO: ${OVERLAP_LEN} bp overlap will be required for sister reads to be merged"

      msg "INFO: Merging paired end reads using FLASH"
      flash \
        -o flash \
        -M ${READ_LEN} \
        -t !{task.cpus} \
        -m ${OVERLAP_LEN} \
        "!{meta.id}_R1.paired.fq" "!{meta.id}_R2.paired.fq"

      # Perform filesize checks and QA report
      echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.Non-overlapping_FastQ_File.tsv"
      for suff in notCombined_1.fastq notCombined_2.fastq; do
        if verify_minimum_file_size "flash.${suff}" 'Non-overlapping FastQ Files' "!{params.min_filesize_non_overlapping_fastq}"; then
          echo -e "!{meta.id}\tNon-overlapping (${suff}) FastQ File\tPASS" \
            >> "!{meta.id}.Non-overlapping_FastQ_File.tsv"
        else
          echo -e "!{meta.id}\tNon-overlapping (${suff}) FastQ File\tFAIL" \
            >> "!{meta.id}.Non-overlapping_FastQ_File.tsv"
        fi
      done

      mv -f flash.notCombined_1.fastq !{meta.id}_R1.paired.fq
      mv -f flash.notCombined_2.fastq !{meta.id}_R2.paired.fq

      CNT_READS_OVERLAPPED=0

      if [ -f flash.extendedFrags.fastq ] && \
      [ -s flash.extendedFrags.fastq ]; then
        CNT_READS_OVERLAPPED=$(awk '{lines++} END{print lines/4}' \
        flash.extendedFrags.fastq)

        cat flash.extendedFrags.fastq >> "!{meta.id}_single.fq"
      else
        # Hack to ensure there's a legit singleton read to pass along to the next steps
        echo "$(tail -n 4 !{meta.id}_R2.paired.fq)" >> "!{meta.id}_single.fq"
      fi

      msg "INFO: ${CNT_READS_OVERLAPPED:-0} pairs overlapped into singleton reads"
    fi

    # Summarize final read set counts
    count_R1=$(wc -l !{meta.id}_R1.paired.fq | awk '{print $1}')
    CNT_CLEANED_PAIRS=$(echo $((${count_R1}/4)))
    msg "INFO: Number of reads cleaned: ${CNT_CLEANED_PAIRS}"

    count_single=$(wc -l "!{meta.id}_single.fq" | awk '{print $1}')
    CNT_CLEANED_SINGLETON=$(echo $((${count_single}/4)))
    msg "INFO: Number of singletons cleaned: ${CNT_CLEANED_SINGLETON}"

    # Report I/O sequence stats
    echo -e "Sample_name\tCleaned_reads_(#_paired)\tCleaned_reads_(#_singletons)\tOverlapped_reads_(#)" \
      > "!{meta.id}.Overlap.tsv"
    echo -e "!{meta.id}\t${CNT_CLEANED_PAIRS}\t${CNT_CLEANED_SINGLETON}\t${CNT_READS_OVERLAPPED:-0}" \
      >> "!{meta.id}.Overlap.tsv"

    # Compress the output FastQ files for outdir storage
    gzip -9f "!{meta.id}_single.fq" \
      "!{meta.id}_R1.paired.fq" \
      "!{meta.id}_R2.paired.fq"

    ### Calculate SHA-512 Checksums of each Input FastQ file ###
    SUMMARY_HEADER=(
      "Sample_name"
      "Checksum_(SHA-512)"
      "File"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.Clean_Reads_FastQ.SHA512-checksums.tsv"

    # Calculate checksums
    for f in "!{meta.id}_R1.paired.fq.gz" "!{meta.id}_R2.paired.fq.gz" "!{meta.id}_single.fq.gz"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.Clean_Reads_FastQ.SHA512-checksums.tsv"
      zcat "${f}" | awk 'NR%2==0' | paste - - | sort -k1,1 | sha512sum | awk '{print $1 "\t" "'"${f}"'"}'
    done >> "!{meta.id}.Clean_Reads_FastQ.SHA512-checksums.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
        flash: $(flash --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
