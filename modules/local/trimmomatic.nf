process TRIMMOMATIC {

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    label "process_high"
    
    container "snads/trimmomatic@sha256:afbb19fdf540e6bd508b657e8dafffb6b411b5b0bf0e302347889220a0b571f1"

    input:
        path ADAPTERS
        path noPhiX_R1
        path noPhiX_R2
        path outpath
        val base
        val size

    output:
        path "*R1.paired.fq", emit: R1_paired
        path "*R2.paired.fq", emit: R2_paired
        path "*trimmo.tsv"
        path "*single.fq"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''

        source bash_functions.sh

        # Adapter clip and quality trim
        verify_file_minimum_size !{ADAPTERS} 'adapters' '10k'

        msg "INFO: Running trimmomatic with !{task.cpus} threads"

        trimmomatic PE -phred33 -threads !{task.cpus}\
        !{noPhiX_R1} !{noPhiX_R2}\
        !{base}_R1.paired.fq !{base}_R1.unpaired.fq\
        !{base}_R2.paired.fq !{base}_R2.unpaired.fq\
        ILLUMINACLIP:!{ADAPTERS}:2:20:10:8:TRUE\
        SLIDINGWINDOW:6:30 LEADING:10 TRAILING:10 MINLEN:50

        TRIMMO_DISCARD=$(grep '^Input Read Pairs: ' .command.err \
        | grep ' Dropped: ' | awk '{print $20}')

        msg "INFO: $TRIMMO_DISCARD reads are poor quality and were discarded" >&2

        CNT_BROKEN_R1=$(awk '{lines++} END{print lines/4}' \
        !{base}_R1.unpaired.fq)
        CNT_BROKEN_R2=$(awk '{lines++} END{print lines/4}' \
        !{base}_R2.unpaired.fq)

        if [[ -z "${TRIMMO_DISCARD}" || -z "${CNT_BROKEN_R1}" || -z "${CNT_BROKEN_R2}" ]]; then
            msg 'ERROR: unable to parse discarded read counts from trimmomatic log' >&2
            exit 1
        fi

        CNT_BROKEN=$((${CNT_BROKEN_R1} + ${CNT_BROKEN_R2}))

        msg "INFO: $CNT_BROKEN_R1 forward reads lacked a high quality R2 sister read" >&2
        msg "INFO: $CNT_BROKEN_R2 reverse reads lacked a high quality R1 sister read" >&2
        msg "INFO: $CNT_BROKEN total broken read pairs were saved as singletons" >&2
        
        echo -e "!{base}\t${TRIMMO_DISCARD} reads Discarded\t${CNT_BROKEN} reads Singletons" \
        > !{base}_trimmo.tsv

        cat !{base}_R1.unpaired.fq !{base}_R2.unpaired.fq > !{base}_single.fq

        rm -f !{base}_R1.unpaired.fq !{base}_R2.unpaired.fq

        minimum_size=$(( !{size}/120 ))
        for suff in R1.paired.fq R2.paired.fq ; do
            verify_file_minimum_size "!{base}_${suff}" 'cleaned read' ${minimum_size}c
        done

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            trimmomatic: $(trimmomatic -version)
        END_VERSIONS

        '''
}
