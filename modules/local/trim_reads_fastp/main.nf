process TRIM_READS_FASTP {

    label "process_high"
    tag { "${meta.id}" }
    container "staphb/fastp@sha256:98bb2bb94bbce4104f7fbbdba72c33c827f5add7cf08cc59fd365c6d82ee4014"

    input:
    tuple val(meta), path(noPhiX)
    path(adapter_reference_file)

    output:
    tuple val(meta), path("${meta.id}.Adapter*_Fast*_File.tsv"), emit: qc_filecheck  // regex grabs 2 QC Files here
    tuple val(meta), path("${meta.id}*{paired,single}.fq")     , emit: fastq_adapters_removed
    path("${meta.id}.Fastp.tsv")                               , emit: summary
    path("${meta.id}.fastp.*")
    path("${meta.id}.Trim_FastQ.SHA512-checksums.tsv")         , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                       , emit: versions

    shell:
    adapter_fasta = adapter_reference_file ? "--adapter_fasta ${adapter_reference_file}" : ""
    '''
    source bash_functions.sh

    # Verify adapter reference file size if provided
    if [[ -s "!{adapter_reference_file}" ]]; then
      if verify_minimum_file_size !{adapter_reference_file} 'Adapters FastA' "!{params.min_filesize_adapters}"; then
        echo -e "!{meta.id}\tAdapters FastA File\tPASS" > !{meta.id}.Adapters_FastA_File.tsv
      else
        echo -e "!{meta.id}\tAdapters FastA File\tFAIL" > !{meta.id}.Adapters_FastA_File.tsv
        exit 1
      fi
    fi

    # Adapter clip and quality trim
    msg "INFO: Performing read trimming on !{meta.id} with Fastp ..."

    # Run fastp
    fastp \
      --in1 !{noPhiX[0]} \
      --in2 !{noPhiX[1]} \
      --out1 !{meta.id}_R1.paired.fq \
      --out2 !{meta.id}_R2.paired.fq \
      --unpaired1 !{meta.id}_R1.unpaired.fq \
      --unpaired2 !{meta.id}_R2.unpaired.fq \
      --length_required !{params.fastp_minimum_read_sequence_length} \
      --cut_right \
      --cut_right_window_size !{params.fastp_window_size} \
      --cut_right_mean_quality !{params.fastp_window_mean_quality} \
      --detect_adapter_for_pe \
      !{adapter_fasta} \
      --json !{meta.id}.fastp.json \
      --html !{meta.id}.fastp.html \
      --thread !{task.cpus}

    msg "INFO: Completed read trimming on !{meta.id} with Fastp"

    # Parse input, discard, and output counts (I really wish `jq` was in this container!)
    NUM_INPUT_READS=$(
                      grep -A 10 '"before_filtering"' "!{meta.id}.fastp.json" \
                        | grep "total_reads" \
                        | sed 's/.*://1;s/,//1'
                      )

    NUM_INPUT_BASES=$(
                      grep -A 10 '"before_filtering"' "!{meta.id}.fastp.json" \
                        | grep "total_bases" \
                        | sed 's/.*://1;s/,//1'
                      )

    NUM_OUTPUT_READS=$(
                        grep -A 10 '"after_filtering"' "!{meta.id}.fastp.json" \
                          | grep "total_reads" \
                          | sed 's/.*://1;s/,//1'
                      )

    NUM_OUTPUT_BASES=$(
                        grep -A 10 '"after_filtering"' "!{meta.id}.fastp.json" \
                          | grep "total_bases" \
                          | sed 's/.*://1;s/,//1'
                      )

    NUM_REMOVED_READS=$((${NUM_INPUT_READS} - ${NUM_OUTPUT_READS}))

    PERCENT_REMOVED_READS=$(echo "${NUM_OUTPUT_READS}" "${NUM_INPUT_READS}" \
      | awk '{proportion=$1/$2} END{printf("%.6f", 100-(proportion*100))}')

    NUM_REMOVED_BASES=$((${NUM_INPUT_BASES} - ${NUM_OUTPUT_BASES}))

    PERCENT_REMOVED_BASES=$(echo "${NUM_OUTPUT_BASES}" "${NUM_INPUT_BASES}" \
      | awk '{proportion=$1/$2} END{printf("%.6f", 100-(proportion*100))}')

    PERCENT_OUTPUT_READS=$(echo "${NUM_REMOVED_READS}" "${NUM_INPUT_READS}" \
      | awk '{proportion=$1/$2} END{printf("%.6f", 100-(proportion*100))}')

    PERCENT_OUTPUT_BASES=$(echo "${NUM_REMOVED_BASES}" "${NUM_INPUT_BASES}" \
      | awk '{proportion=$1/$2} END{printf("%.6f", 100-(proportion*100))}')

    # Form and create a summary file of input, discarded, and output
    SUMMARY_HEADER=(
      "Sample_name"
      "Input_reads_(#)"
      "Input_basepairs_(#)"
      "Removed_reads_(#)"
      "Removed_reads_(%)"
      "Removed_basepairs_(#)"
      "Removed_basepairs_(%)"
      "Output_reads_(#)"
      "Output_reads_(%)"
      "Output_basepairs_(#)"
      "Output_basepairs_(%)"
    )

    SUMMARY_OUTPUT=(
      "!{meta.id}"
      "${NUM_INPUT_READS}"
      "${NUM_INPUT_BASES}"
      "${NUM_REMOVED_READS}"
      "${PERCENT_REMOVED_READS}"
      "${NUM_REMOVED_BASES}"
      "${PERCENT_REMOVED_BASES}"
      "${NUM_OUTPUT_READS}"
      "${PERCENT_OUTPUT_READS}"
      "${NUM_OUTPUT_BASES}"
      "${PERCENT_OUTPUT_BASES}"
    )

    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//1')
    SUMMARY_OUTPUT=$(printf "%s\t" "${SUMMARY_OUTPUT[@]}" | sed 's/\t$//1')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.Fastp.tsv"
    echo "${SUMMARY_OUTPUT}" >> "!{meta.id}.Fastp.tsv"

    # Test/verify paired FastQ outfiles sizes are reasonable to continue
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > !{meta.id}.Adapter_and_QC_Trimmed_FastQ_File.tsv
    for suff in R1.paired.fq R2.paired.fq; do
      if verify_minimum_file_size "!{meta.id}_${suff}" 'Adapter-removed FastQ Files' "!{params.min_filesize_fastq_adapters_removed}"; then
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tPASS" \
          >> !{meta.id}.Adapter_and_QC_Trimmed_FastQ_File.tsv
      else
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tFAIL" \
          >> !{meta.id}.Adapter_and_QC_Trimmed_FastQ_File.tsv
      fi
    done

    # Merge broken sister FastQ reads into one singleton FastQ file
    cat !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq >> !{meta.id}_single.fq
    rm -f !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq

    ### Calculate SHA-512: Checksums of each FastQ file ###
    msg "INFO: Calculating checksums for !{meta.id}_R1.paired.fq and !{meta.id}_R2.paired.fq ..."

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

    msg "INFO: Calculated checksums for !{meta.id}_R1.paired.fq and !{meta.id}_R2.paired.fq"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
        fastp: $(fastp --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
