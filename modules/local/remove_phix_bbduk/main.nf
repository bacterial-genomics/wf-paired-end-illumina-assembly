process REMOVE_PHIX_BBDUK {

    publishDir "${params.outdir}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{raw,phix}.tsv"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{PhiX_Genome_File,PhiX-removed_FastQ_Files}.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}" }

    label "process_low"
    tag { "${prefix}" }

    container "snads/bbtools@sha256:9f2a9b08563839cec87d856f0fc7607c235f464296fd71e15906ea1d15254695"

    input:
    tuple val(prefix), path(reads), path(qc_input_filecheck)

    output:
    tuple val(prefix), path("${prefix}_noPhiX-R1.fsq"), path("${prefix}_noPhiX-R2.fsq"), path("*File*.tsv"), emit: phix_removed
    path "${prefix}.PhiX_Genome_File.tsv", emit: qc_phix_genome_filecheck
    path "${prefix}.PhiX-removed_FastQ_Files.tsv", emit: qc_phix_removed_filecheck
    path "${prefix}.raw.tsv"
    path "${prefix}.phix.tsv"
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_input_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Get PhiX, check if it exists, and verify file size
    if ! check_if_file_exists_allow_seconds !{params.phix_reference} '60'; then
      exit 1
    fi

    if verify_minimum_file_size !{params.phix_reference} 'PhiX Genome' "!{params.min_filesize_phix_genome}"; then
      echo -e "!{prefix}\tPhiX Genome\tPASS" >> !{prefix}.PhiX_Genome_File.tsv
    else
      echo -e "!{prefix}\tPhiX Genome\tFAIL" >> !{prefix}.PhiX_Genome_File.tsv
    fi

    # Remove PhiX
    msg "INFO: Running bbduk with !{task.cpus} threads"

    bbduk.sh \
      threads=!{task.cpus} \
      k=31 \
      hdist=1 \
      ref="!{params.phix_reference}" \
      in="!{reads[0]}" \
      in2="!{reads[1]}" \
      out=!{prefix}_noPhiX-R1.fsq \
      out2=!{prefix}_noPhiX-R2.fsq \
      qin=auto \
      qout=33 \
      overwrite=t

    for suff in R1.fsq R2.fsq; do
      if verify_minimum_file_size "!{prefix}_noPhiX-${suff}" 'PhiX-removed FastQ Files' "!{params.min_filesize_fastq_phix_removed}"; then
        echo -e "!{prefix}\tPhiX-removed FastQ ($suff) File\tPASS" \
          >> !{prefix}.PhiX-removed_FastQ_Files.tsv
      else
        echo -e "!{prefix}\tPhiX-removed FastQ ($suff) File\tFAIL" \
          >> !{prefix}.PhiX-removed_FastQ_Files.tsv
      fi
    done

    TOT_READS=$(grep '^Input: ' .command.err | awk '{print $2}')
    TOT_BASES=$(grep '^Input: ' .command.err | awk '{print $4}')

    if [[ -z "${TOT_READS}" || -z "${TOT_BASES}" ]]; then
      msg 'ERROR: unable to parse input counts from bbduk log' >&2
      exit 1
    fi

    PHIX_READS=$(grep '^Contaminants: ' .command.err | awk '{print $2}' | sed 's/,//g')
    PHIX_BASES=$(grep '^Contaminants: ' .command.err | awk '{print $5}' | sed 's/,//g')

    msg "INFO: ${TOT_BASES} bp and $TOT_READS reads provided as raw input"
    msg "INFO: ${PHIX_BASES:-0} bp of PhiX were detected and removed in ${PHIX_READS:-0} reads"

    echo -e "!{prefix}\t${TOT_BASES} bp Raw\t${TOT_READS} reads Raw" \
      > !{prefix}.raw.tsv
    echo -e "!{prefix}\t${PHIX_BASES:-0} bp PhiX\t${PHIX_READS:-0} reads PhiX" \
      > !{prefix}.phix.tsv

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      bbduk: $(bbduk.sh --version 2>&1 | head -n 2 | tail -1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
