process TRIM_READS_TRIMMOMATIC {

    // errorStrategy 'terminate'

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.trimmo.tsv"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{Adapters_FastA,Adapter-removed_FastQ_Files}.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    label "process_high"
    tag { "${base}" }
    
    container "snads/trimmomatic@sha256:afbb19fdf540e6bd508b657e8dafffb6b411b5b0bf0e302347889220a0b571f1"

    input:
        tuple val(base), path(noPhiX_R1), path(noPhiX_R2), path(phix_qc_filechecks)

    output:
        tuple val(base), path("${base}_R1.paired.fq"), path("${base}_R2.paired.fq"), path("*Adapter*.tsv"), emit: trimmo
        path "${base}.Adapters_FastA.tsv", emit: qc_adapters_filecheck
        path "${base}.Adapter-removed_FastQ_Files.tsv", emit: qc_removed_adapters_filecheck
        path "${base}.trimmo.tsv"
        path "${base}.single.fq"
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
        
        # Get Adapters, check if it exists, and verify file size
        ADAPTERS="${DIR}/adapters_Nextera_NEB_TruSeq_NuGEN_ThruPLEX.fas"
        if ! check_if_file_exists_allow_seconds ${ADAPTERS} '60'; then
          exit 1
        fi
        if verify_minimum_file_size ${ADAPTERS} 'Adapters FastA' "!{params.min_filesize_adapters}"; then
          echo -e "!{base}\tAdapters FastA File\tPASS" > !{base}.Adapters_FastA.tsv
        else
          echo -e "!{base}\tAdapters FastA File\tFAIL" > !{base}.Adapters_FastA.tsv
        fi

        # Adapter clip and quality trim
        msg "INFO: Running trimmomatic with !{task.cpus} threads"

        trimmomatic PE \
        -phred33 \
        -threads !{task.cpus} \
        !{noPhiX_R1} !{noPhiX_R2} \
        !{base}_R1.paired.fq !{base}_R1.unpaired.fq \
        !{base}_R2.paired.fq !{base}_R2.unpaired.fq \
        ILLUMINACLIP:${ADAPTERS}:2:20:10:8:TRUE \
        SLIDINGWINDOW:6:30 \
        LEADING:10 \
        TRAILING:10 \
        MINLEN:50

        TRIMMO_DISCARD=$(grep '^Input Read Pairs: ' .command.err \
        | grep ' Dropped: ' | awk '{print $20}')

        msg "INFO: ${TRIMMO_DISCARD} reads are poor quality and were discarded" >&2

        CNT_BROKEN_R1=$(awk '{lines++} END{print lines/4}' !{base}_R1.unpaired.fq)
        CNT_BROKEN_R2=$(awk '{lines++} END{print lines/4}' !{base}_R2.unpaired.fq)

        if [[ -z "${TRIMMO_DISCARD}" || -z "${CNT_BROKEN_R1}" || -z "${CNT_BROKEN_R2}" ]]; then
          msg 'ERROR: unable to parse discarded read counts from trimmomatic log' >&2
          exit 1
        fi

        CNT_BROKEN=$((${CNT_BROKEN_R1} + ${CNT_BROKEN_R2}))

        msg "INFO: $CNT_BROKEN_R1 forward reads lacked a high quality R2 sister read" >&2
        msg "INFO: $CNT_BROKEN_R2 reverse reads lacked a high quality R1 sister read" >&2
        msg "INFO: $CNT_BROKEN total broken read pairs were saved as singletons" >&2
        
        echo -e "!{base}\t${TRIMMO_DISCARD} reads Discarded\t${CNT_BROKEN} reads Singletons" \
        > !{base}.trimmo.tsv

        cat !{base}_R1.unpaired.fq !{base}_R2.unpaired.fq > !{base}.single.fq

        rm -f !{base}_R1.unpaired.fq !{base}_R2.unpaired.fq

        for suff in R1.paired.fq R2.paired.fq; do
          if verify_minimum_file_size "!{base}_${suff}" 'Adapter-removed FastQ Files' "!{params.min_filesize_fastq_adapters_removed}"; then
            echo -e "!{base}\tAdapter-removed ($suff) FastQ File\tPASS" \
              >> !{base}.Adapter-removed_FastQ_Files.tsv
          else
            echo -e "!{base}\tAdapter-removed ($suff) FastQ File\tFAIL" \
              >> !{base}.Adapter-removed_FastQ_Files.tsv
          fi
        done

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            trimmomatic: $(trimmomatic -version)
        END_VERSIONS
        '''
}
