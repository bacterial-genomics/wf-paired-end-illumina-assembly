process TRIM_READS_TRIMMOMATIC {

    label "process_high"
    tag { "${meta.id}" }
    container "snads/trimmomatic@sha256:afbb19fdf540e6bd508b657e8dafffb6b411b5b0bf0e302347889220a0b571f1"

    input:
    tuple val(meta), path(noPhiX)
    path adapter_reference_file

    output:
    tuple val(meta), path("${meta.id}.Adapter*_File.tsv") , emit: qc_filecheck
    tuple val(meta), path("${meta.id}*{paired,single}.fq"), emit: fastq_adapters_removed
    path("${meta.id}.Trimmomatic.tsv")
    path(".command.{out,err}")
    path("versions.yml")                                  , emit: versions

    shell:
    keep_both_reads           = params.trimmomatic_keep_both_reads                                                      ? 'TRUE'                                       : 'FALSE'
    phred                     = (params.trimmomatic_phred == 64)                                                        ? '-phred64'                                   : '-phred33'
    min_length                = (params.trimmomatic_min_length >= 1)                                                    ? params.trimmomatic_min_length                : 50
    window_size               = (params.trimmomatic_window_size >= 1)                                                   ? params.trimmomatic_window_size               : 6
    seed_mismatches           = (params.trimmomatic_seed_mismatches >= 1)                                               ? params.trimmomatic_seed_mismatches           : 2
    leading_quality           = (params.trimmomatic_leading_quality >= 1 && params.trimmomatic_leading_quality <= 50)   ? params.trimmomatic_leading_quality           : 10
    req_quality               = (params.trimmomatic_required_quality >= 0 && params.trimmomatic_required_quality <= 50) ? params.trimmomatic_required_quality          : 30
    trailing_quality          = (params.trimmomatic_trailing_quality >= 1 && params.trimmomatic_trailing_quality <= 50) ? params.trimmomatic_trailing_quality          : 10
    min_adapter_length        = (params.trimmomatic_min_adapter_length >= 1)                                            ? params.trimmomatic_min_adapter_length        : 8
    simple_clip_threshold     = (params.trimmomatic_simple_clip_threshold >= 1)                                         ? params.trimmomatic_simple_clip_threshold     : 10
    palindrome_clip_threshold = (params.trimmomatic_palindrome_clip_threshold >= 1)                                     ? params.trimmomatic_palindrome_clip_threshold : 20

    illumina_clip_params      = "${seed_mismatches}:${palindrome_clip_threshold}:${simple_clip_threshold}:${min_adapter_length}:${keep_both_reads}"
    '''
    source bash_functions.sh

    # Verify adapter reference file size
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}.Adapters_FastA_File.tsv"
    if verify_minimum_file_size !{adapter_reference_file} 'Adapters FastA' "!{params.min_filesize_adapters}"; then
      echo -e "!{meta.id}\tAdapters FastA File\tPASS" >> "!{meta.id}.Adapters_FastA_File.tsv"
    else
      echo -e "!{meta.id}\tAdapters FastA File\tFAIL" >> "!{meta.id}.Adapters_FastA_File.tsv"
    fi

    # Adapter clip and quality trim
    msg "INFO: Performing read trimming with Trimmomatic"

    trimmomatic PE \
      !{phred} \
      -threads !{task.cpus} \
      !{noPhiX[0]} !{noPhiX[1]} \
      !{meta.id}_R1.paired.fq !{meta.id}_R1.unpaired.fq \
      !{meta.id}_R2.paired.fq !{meta.id}_R2.unpaired.fq \
      MINLEN:!{min_length} \
      LEADING:!{leading_quality} \
      TRAILING:!{trailing_quality} \
      SLIDINGWINDOW:!{window_size}:!{req_quality} \
      ILLUMINACLIP:!{adapter_reference_file}:!{illumina_clip_params}

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

    echo -e "!{meta.id}\t${TRIMMO_DISCARD}\t${CNT_BROKEN}" \
    > "!{meta.id}.Trimmomatic.tsv"

    sed -i '1i Sample name\t# discarded reads\t# singleton reads' !{meta.id}.Trimmomatic.tsv

    cat !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq > "!{meta.id}_single.fq"

    rm -f !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}.Adapter-removed_FastQ_File.tsv"
    for suff in R1.paired.fq R2.paired.fq; do
      if verify_minimum_file_size "!{meta.id}_${suff}" 'Adapter-removed FastQ Files' "!{params.min_filesize_fastq_adapters_removed}"; then
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tPASS" \
          >> "!{meta.id}.Adapter-removed_FastQ_File.tsv"
      else
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tFAIL" \
          >> "!{meta.id}.Adapter-removed_FastQ_File.tsv"
      fi
    done

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        trimmomatic: $(trimmomatic -version)
    END_VERSIONS
    '''
}
