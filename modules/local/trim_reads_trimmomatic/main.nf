process TRIM_READS_TRIMMOMATIC {

    label "process_high"
    tag { "${meta.id}" }
    container "snads/trimmomatic@sha256:afbb19fdf540e6bd508b657e8dafffb6b411b5b0bf0e302347889220a0b571f1"

    input:
    tuple val(meta), path(noPhiX)
    path adapter_reference_file

    output:
    path(".command.{out,err}")
    path "${meta.id}_single.fq"
    path "${meta.id}.trimmomatic.tsv"
    path "versions.yml"                                  , emit: versions
    tuple val(meta), path("${meta.id}.Adapter*_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}_R{1,2}.paired.fq") , emit: fastq_adapters_removed

    shell:
    '''
    source bash_functions.sh

    # Verify adapter reference file size
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > !{meta.id}.Adapters_FastA_File.tsv
    if verify_minimum_file_size !{adapter_reference_file} 'Adapters FastA' "!{params.min_filesize_adapters}"; then
      echo -e "!{meta.id}\tAdapters FastA File\tPASS" >> "!{meta.id}.Adapters_FastA_File.tsv"
    else
      echo -e "!{meta.id}\tAdapters FastA File\tFAIL" >> "!{meta.id}.Adapters_FastA_File.tsv"
    fi

    # Adapter clip and quality trim
    msg "INFO: Performing read trimming with Trimmomatic"

    trimmomatic PE \
      -phred33 \
      -threads !{task.cpus} \
      !{noPhiX[0]} !{noPhiX[1]} \
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

    echo -e "!{meta.id}\t${TRIMMO_DISCARD}\t${CNT_BROKEN}" \
    > !{meta.id}.trimmomatic.tsv

    sed -i '1i Sample name\t# discarded reads\t# singleton reads' !{meta.id}.trimmomatic.tsv

    cat !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq > !{meta.id}_single.fq

    rm -f !{meta.id}_R1.unpaired.fq !{meta.id}_R2.unpaired.fq

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > !{meta.id}.Adapter-removed_FastQ_File.tsv
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
