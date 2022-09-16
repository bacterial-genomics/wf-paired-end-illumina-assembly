process EXTRACT_SINGLETONS {

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.gz"
    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    label "process_low"

    container "snads/flash@sha256:363b2f44d040c669191efbc3d3ba99caf5efd3fdef370af8f00f3328932143a6"

    input:
        path input
        path R1_paired
        path R2_paired
        val base
        val size

    output:
        path "*R1.paired.fq.gz", emit: R1_paired_gz
        path "*R2.paired.fq.gz", emit: R2_paired_gz
        path "*single.fq.gz", emit: single_gz
        path "*overlap.tsv"
        path "*clean-reads.tsv"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''

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
            flash -m ${OVERLAP_LEN} -M ${READ_LEN} -o flash -t !{task.cpus} !{R1_paired} !{R2_paired}

            minimum_size=$(( !{size}/150 ))
            for suff in flash.notCombined_1.fastq flash.notCombined_2.fastq ; do
                verify_file_minimum_size "${suff}" 'cleaned non-overlapping read' ${minimum_size}c
            done

            rm !{R1_paired} !{R2_paired}
            mv flash.notCombined_1.fastq !{base}_R1.paired.fq
            mv flash.notCombined_2.fastq !{base}_R2.paired.fq

            if [ -f  flash.extendedFrags.fastq ] && \
                [ -s  flash.extendedFrags.fastq ]; then
                CNT_READS_OVERLAPPED=$(awk '{lines++} END{print lines/4}' \
                flash.extendedFrags.fastq)

                cat flash.extendedFrags.fastq >> !{base}_single.fq
                rm flash.extendedFrags.fastq
            fi

            msg "INFO: ${CNT_READS_OVERLAPPED:-0} pairs overlapped into singleton reads" >&2
            echo -e "!{base}\t${CNT_READS_OVERLAPPED:-0} reads Overlapped" \
            > !{base}_overlap.tsv
        fi

        # Summarize final read set and compress

        count_R1=$(echo $(cat !{base}_R1.paired.fq | wc -l))
        CNT_CLEANED_PAIRS=$(echo $((${count_R1}/4)))
        msg "INFO: CNT_CLEANED_PAIRS ${CNT_CLEANED_PAIRS}"

        count_single=$(echo $(cat !{base}_single.fq | wc -l))
        CNT_CLEANED_SINGLETON=$(echo $((${count_single}/4)))
        msg "INFO: CNT_CLEANED_SINGLETON ${CNT_CLEANED_SINGLETON}"
        

        echo -e "!{base}\t${CNT_CLEANED_PAIRS} cleaned pairs\t${CNT_CLEANED_SINGLETON} cleaned singletons" \
        > !{base}_clean-reads.tsv

        gzip !{base}_single.fq\
        !{base}_R1.paired.fq\
        !{base}_R2.paired.fq

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            flash: $(flash --version | head -n 1 | awk 'NF>1{print $NF}')
        END_VERSIONS

        '''
}