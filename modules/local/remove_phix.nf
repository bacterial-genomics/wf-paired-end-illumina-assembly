process REMOVE_PHIX {

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    label "process_low"
    tag { "${base}" }
    
    container "snads/bbtools@sha256:9f2a9b08563839cec87d856f0fc7607c235f464296fd71e15906ea1d15254695"
    
    input:
        tuple val(base), val(size), path(input)

    output:
        tuple val(base), val(size), path("${base}_noPhiX-R1.fsq"), path("${base}_noPhiX-R2.fsq"), emit: nophix
        path "${base}.raw.tsv"
        path "${base}.phix.tsv"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        source bash_functions.sh

        # Get PhiX, check if it exists, and verify file size
        PHIX="${DIR}/PhiX_NC_001422.1.fasta"
        check_if_file_exists_allow_seconds ${PHIX} '60'
        verify_file_minimum_size ${PHIX} 'PhiX genome' '5k' '100'

        # Remove PhiX
        msg "INFO: Running bbduk with !{task.cpus} threads"
        
        bbduk.sh threads=!{task.cpus} k=31 hdist=1\
        ref="${PHIX}" in="!{input[0]}" in2="!{input[1]}"\
        out=!{base}_noPhiX-R1.fsq out2=!{base}_noPhiX-R2.fsq\
        qin=auto qout=33 overwrite=t

        for suff in R1.fsq R2.fsq; do
            verify_file_minimum_size "!{base}_noPhiX-${suff}" 'PhiX cleaned read' "!{size}" '0.8'
        done

        TOT_READS=$(grep '^Input: ' .command.err \
        | awk '{print $2}')
        TOT_BASES=$(grep '^Input: ' .command.err \
        | awk '{print $4}')

        if [[ -z "${TOT_READS}" || -z "${TOT_BASES}" ]]; then
            msg 'ERROR: unable to parse input counts from bbduk log' >&2
            exit 1
        fi

        PHIX_READS=$(grep '^Contaminants: ' .command.err \
        | awk '{print $2}' | sed 's/,//g')
        PHIX_BASES=$(grep '^Contaminants: ' .command.err \
        | awk '{print $5}' | sed 's/,//g')

        msg "INFO: ${TOT_BASES} bp and $TOT_READS reads provided as raw input"
        msg "INFO: ${PHIX_BASES:-0} bp of PhiX were detected and removed in ${PHIX_READS:-0} reads"

        echo -e "!{base}\t${TOT_BASES} bp Raw\t${TOT_READS} reads Raw" \
        > !{base}.raw.tsv
        echo -e "!{base}\t${PHIX_BASES:-0} bp PhiX\t${PHIX_READS:-0} reads PhiX" \
        > !{base}.phix.tsv

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            bbduk: $(bbduk.sh --version 2>&1 | head -n 2 | tail -1 | awk 'NF>1{print $NF}')
        END_VERSIONS
        '''
}