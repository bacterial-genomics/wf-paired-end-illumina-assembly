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
    path("${meta.id}.Trim_FastQ.SHA256-checksums.tsv")         , emit: checksums
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

    # Calculate number of reads discarded
    READ_COUNT_INPUT=$(
                      grep -A3 '"before_filtering"' !{meta.id}.fastp.json \
                        | grep "total_reads" \
                        | sed 's/.*://g;s/,//g'
                      )

    READ_COUNT_OUTPUT=$(
                        grep -A3 '"after_filtering"' !{meta.id}.fastp.json \
                          | grep "total_reads" \
                          | sed 's/.*://g;s/,//g'
                      )

    READ_COUNT_DISCARDED=$((${READ_COUNT_INPUT} - ${READ_COUNT_OUTPUT}))

    msg "INFO: ${READ_COUNT_DISCARDED} reads are poor quality and were discarded"

    # Count up the total number of broken sister reads
    COUNT_BROKEN_R1=$(awk '{lines++} END{print lines/4}' !{meta.id}_R1.unpaired.fq)
    COUNT_BROKEN_R2=$(awk '{lines++} END{print lines/4}' !{meta.id}_R2.unpaired.fq)
    COUNT_BROKEN_TOTAL=$((${COUNT_BROKEN_R1} + ${COUNT_BROKEN_R2}))

    # Log report the total counts of singletons
    msg "INFO: ${COUNT_BROKEN_R1} forward reads lacked a high quality R2 sister read"
    msg "INFO: ${COUNT_BROKEN_R2} reverse reads lacked a high quality R1 sister read"
    msg "INFO: ${COUNT_BROKEN_TOTAL} total broken read pairs were saved as singletons"

    # Create report file of reads removed and broken
    echo -e "!{meta.id}\t${READ_COUNT_DISCARDED}\t${COUNT_BROKEN_TOTAL}" > !{meta.id}.Fastp.tsv
    sed -i '1i Sample_name\tDiscarded_reads_(#)\tSingleton_reads_(#)' !{meta.id}.Fastp.tsv

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

    ### Calculate SHA-256 Checksums of each FastQ file ###
    msg "INFO: Calculating checksums for !{meta.id}_R1.paired.fq and !{meta.id}_R2.paired.fq ..."

    SUMMARY_HEADER=(
      "Sample_name"
      "Checksum_(SHA-256)"
      "File"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.Trim_FastQ.SHA256-checksums.tsv"

    # Calculate checksums
    for f in "!{meta.id}_R1.paired.fq" "!{meta.id}_R2.paired.fq" "!{meta.id}_single.fq"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.Trim_FastQ.SHA256-checksums.tsv"
      awk 'NR%2==0'  "${f}" | paste - - | sort -k1,1 | sha256sum | awk '{print $1 "\t" "'"${f}"'"}'
    done >> "!{meta.id}.Trim_FastQ.SHA256-checksums.tsv"

    msg "INFO: Calculated checksums for !{meta.id}_R1.paired.fq and !{meta.id}_R2.paired.fq"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha256sum: $(sha256sum --version | grep "^sha256sum" | sed 's/sha256sum //1')
        fastp: $(fastp --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
