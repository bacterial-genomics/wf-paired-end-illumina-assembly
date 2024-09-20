process TRIM_READS_TRIMMOMATIC {

    label "process_high"
    tag { "${meta.id}" }
    container "staphb/trimmomatic@sha256:57b673e66313e355a447e4fa1a78fd3ba1ae3ddd8c8f91358efe99140acb5ddb"

    input:
    tuple val(meta), path(reads)
    path(adapter_reference_file)

    output:
    tuple val(meta), path("${meta.id}.Adapter*_Fast*_File.tsv"), emit: qc_filecheck  // regex grabs 2 QC Files here
    tuple val(meta), path("${meta.id}*{paired,single}.fq")     , emit: fastq_adapters_removed
    path("${meta.id}.Trimmomatic.tsv")                         , emit: summary
    path("${meta.id}.Trim_FastQ.SHA512-checksums.tsv")         , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                       , emit: versions

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
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.Adapters_FastA_File.tsv"
    if verify_minimum_file_size "!{adapter_reference_file}" 'Adapters FastA' "!{params.min_filesize_adapters}"; then
      echo -e "!{meta.id}\tAdapters FastA File\tPASS" >> "!{meta.id}.Adapters_FastA_File.tsv"
    else
      echo -e "!{meta.id}\tAdapters FastA File\tFAIL" >> "!{meta.id}.Adapters_FastA_File.tsv"
    fi

    # Adapter clip and quality trim
    msg "INFO: Performing read trimming on !{meta.id} with Trimmomatic ..."

    # NOTE: *order* matters on trimming here with Trimmomatic!!!
    trimmomatic PE \
      "!{phred}" \
      -threads "!{task.cpus}" \
      "!{reads[0]}" "!{reads[1]}" \
      "!{meta.id}_R1.paired.fq" "!{meta.id}_R1.unpaired.fq" \
      "!{meta.id}_R2.paired.fq" "!{meta.id}_R2.unpaired.fq" \
      ILLUMINACLIP:"!{adapter_reference_file}":"!{illumina_clip_params}" \
      SLIDINGWINDOW:"!{window_size}":"!{req_quality}" \
      LEADING:"!{leading_quality}" \
      TRAILING:"!{trailing_quality}" \
      MINLEN:"!{min_length}"

    msg "INFO: Completed read trimming on !{meta.id} with Trimmomatic"

    cat !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq > "!{meta.id}_single.fq"

    rm -f !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq

    # Parse input, discard, and output counts
    NUM_INPUT_READS=$(grep '^Input Read Pairs: ' .command.err \
      | awk '{print $4}')

    ### NUM_INPUT_BASES=  **missing**; skip slow calc; SeqKit does this on previous process output

    NUM_REMOVED_READS=$(grep '^Input Read Pairs: ' .command.err \
      | grep ' Dropped: ' | awk '{print $20}')

    PERCENT_REMOVED_READS=$(grep '^Input Read Pairs: ' .command.err \
      | grep ' Dropped: ' | awk '{print $21}' | tr -d '()%')

    NUM_OUTPUT_PAIRED_READS=$(wc -l "!{meta.id}_R1.paired.fq" | awk '{print $1/2}')

    NUM_OUTPUT_SINGLE_READS=$(wc -l "!{meta.id}_single.fq" | awk '{print $1/4}')

    NUM_OUTPUT_READS=$((${NUM_OUTPUT_PAIRED_READS} + ${NUM_OUTPUT_SINGLE_READS}))

    ### NUM_REMOVED_BASES=  **missing**; skip slow calc
    ### PERCENT_REMOVED_BASES=  **missing**; skip slow calc

    PERCENT_OUTPUT_READS=$(echo "${NUM_REMOVED_READS}" "${NUM_INPUT_READS}" \
      | awk '{proportion=$1/$2} END{printf("%.6f", 100-(proportion*100))}')

    ### NUM_OUTPUT_BASES=  **missing**; skip slow calc; SeqKit does this in next process input
    ### PERCENT_OUTPUT_BASES=  **missing**; skip slow calc; SeqKit does this in next process input

    msg "INFO: ${NUM_REMOVED_READS} reads (${PERCENT_REMOVED_READS}% of input) were discarded"

    # Form and create a summary file of input, discarded, and output
    SUMMARY_HEADER=(
      "Sample_name"
      "Input_reads_(#)"
      "Removed_reads_(#)"
      "Removed_reads_(%)"
      "Output_reads_(#)"
      "Output_reads_(%)"
    )
      # Skipped these slow calcs (fastp provides these but not trimmomatic)
      # "Input_basepairs_(#)"
      # "Removed_basepairs_(#)"
      # "Removed_basepairs_(%)"
      # "Output_basepairs_(#)"
      # "Output_basepairs_(%)"

    SUMMARY_OUTPUT=(
      "!{meta.id}"
      "${NUM_INPUT_READS}"
      "${NUM_REMOVED_READS}"
      "${PERCENT_REMOVED_READS}"
      "${NUM_OUTPUT_READS}"
      "${PERCENT_OUTPUT_READS}"
    )
      # Skipped these slow calcs (fastp provides these but not trimmomatic)
      # "${NUM_INPUT_BASES}"
      # "${NUM_REMOVED_BASES}"
      # "${PERCENT_REMOVED_BASES}"
      # "${NUM_OUTPUT_BASES}"
      # "${PERCENT_OUTPUT_BASES}"

    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//1')
    SUMMARY_OUTPUT=$(printf "%s\t" "${SUMMARY_OUTPUT[@]}" | sed 's/\t$//1')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.Trimmomatic.tsv"
    echo "${SUMMARY_OUTPUT}" >> "!{meta.id}.Trimmomatic.tsv"

    # Test/verify paired FastQ outfiles sizes are reasonable to continue
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.Adapter_and_QC_Trimmed_FastQ_File.tsv"
    for suff in R1.paired.fq R2.paired.fq; do
      if verify_minimum_file_size "!{meta.id}_${suff}" 'Adapter-removed FastQ Files' "!{params.min_filesize_fastq_adapters_removed}"; then
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tPASS" \
          >> "!{meta.id}.Adapter_and_QC_Trimmed_FastQ_File.tsv"
      else
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tFAIL" \
          >> "!{meta.id}.Adapter_and_QC_Trimmed_FastQ_File.tsv"
      fi
    done

    ### Calculate SHA-512 Checksums of each FastQ file ###
    msg "INFO: Calculating checksums for !{meta.id}_R1.paired.fq and !{meta.id}_R2.paired.fq !{meta.id}_single.fq ..."

    SUMMARY_HEADER=(
      "Sample_name"
      "Checksum_(SHA-512)"
      "File"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.Trim_FastQ.SHA512-checksums.tsv"

    # Calculate checksums
    for f in "!{meta.id}_R1.paired.fq" "!{meta.id}_R2.paired.fq" "!{meta.id}_single.fq"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.Trim_FastQ.SHA512-checksums.tsv"
      awk 'NR%2==0'  "${f}" | paste - - | sort -k1,1 | sha512sum | awk '{print $1 "\t" "'"${f}"'"}'
    done >> "!{meta.id}.Trim_FastQ.SHA512-checksums.tsv"

    msg "INFO: Calculated checksums for !{meta.id}_R1.paired.fq and !{meta.id}_R2.paired.fq !{meta.id}_single.fq"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
        trimmomatic: $(trimmomatic -version)
    END_VERSIONS
    '''
}
