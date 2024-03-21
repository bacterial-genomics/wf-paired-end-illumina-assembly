process TRIM_READS_FASTP {

    label "process_high"
    tag { "${meta.id}" }
    container "staphb/fastp@sha256:98bb2bb94bbce4104f7fbbdba72c33c827f5add7cf08cc59fd365c6d82ee4014"

    input:
    tuple val(meta), path(noPhiX)
    path(adapter_reference_file)

    output:
    tuple val(meta), path("${meta.id}.Adapter*_File.tsv") , emit: qc_filecheck
    tuple val(meta), path("${meta.id}*{paired,single}.fq"), emit: fastq_adapters_removed
    path("${meta.id}.fastp.tsv")                          , emit: summary
    path("${meta.id}.fastp.*")
    path(".command.{out,err}")
    path("versions.yml")                                  , emit: versions

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
    msg "INFO: Performing read trimming with fastp"

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

    msg "INFO: ${READ_COUNT_DISCARDED} reads are poor quality and were discarded" >&2

    # Count up the total number of broken sister reads
    COUNT_BROKEN_R1=$(awk '{lines++} END{print lines/4}' !{meta.id}_R1.unpaired.fq)
    COUNT_BROKEN_R2=$(awk '{lines++} END{print lines/4}' !{meta.id}_R2.unpaired.fq)
    COUNT_BROKEN_TOTAL=$((${COUNT_BROKEN_R1} + ${COUNT_BROKEN_R2}))

    # Log report the total counts of singletons
    msg "INFO: $COUNT_BROKEN_R1 forward reads lacked a high quality R2 sister read" >&2
    msg "INFO: $COUNT_BROKEN_R2 reverse reads lacked a high quality R1 sister read" >&2
    msg "INFO: $COUNT_BROKEN_TOTAL total broken read pairs were saved as singletons" >&2

    # Create report file of reads removed and broken
    echo -e "!{meta.id}\t${READ_COUNT_DISCARDED}\t${COUNT_BROKEN_TOTAL}" > !{meta.id}.fastp.tsv
    sed -i '1i Sample name\t# discarded reads\t# singleton reads' !{meta.id}.fastp.tsv

    # Test/verify paired FastQ outfiles sizes are reasonable to continue
    for suff in R1.paired.fq R2.paired.fq; do
      if verify_minimum_file_size "!{meta.id}_${suff}" 'Adapter-removed FastQ Files' "!{params.min_filesize_fastq_adapters_removed}"; then
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tPASS" \
          >> !{meta.id}.Adapter-removed_FastQ_File.tsv
      else
        echo -e "!{meta.id}\tAdapter-removed ($suff) FastQ File\tFAIL" \
          >> !{meta.id}.Adapter-removed_FastQ_File.tsv
      fi
    done

    # Merge broken sister FastQ reads into one singleton FastQ file
    cat !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq >> !{meta.id}_single.fq
    rm -f !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        fastp: $(fastp --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
