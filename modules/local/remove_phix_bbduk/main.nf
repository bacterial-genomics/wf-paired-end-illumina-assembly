process REMOVE_PHIX_BBDUK {

    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{PhiX_Genome_File,PhiX-removed_FastQ_Files}.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${meta.id}.${task.process}${filename}" }

    label "process_low"
    tag { "${meta.id}" }

    container "snads/bbtools@sha256:9f2a9b08563839cec87d856f0fc7607c235f464296fd71e15906ea1d15254695"

    input:
    tuple val(meta), path(reads), path(qc_input_filecheck)
    path phix_reference_file

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                                                                    , emit: versions
    path "${meta.id}.Summary.PhiX.tab"                                                                     , emit: phix_summary
    path "${meta.id}.PhiX_Genome_File.tsv"                                                                 , emit: qc_phix_genome_filecheck
    path "${meta.id}.PhiX-removed_FastQ_Files.tsv"                                                         , emit: qc_phix_removed_filecheck
    tuple val(meta), path("${meta.id}_noPhiX-R1.fsq"), path("${meta.id}_noPhiX-R2.fsq"), path("*File*.tsv"), emit: phix_removed

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

    # Verify PhiX reference file size
    if verify_minimum_file_size !{phix_reference_file} 'PhiX Genome' "!{params.min_filesize_phix_genome}"; then
      echo -e "!{meta.id}\tPhiX Genome\tPASS" >> !{meta.id}.PhiX_Genome_File.tsv
    else
      echo -e "!{meta.id}\tPhiX Genome\tFAIL" >> !{meta.id}.PhiX_Genome_File.tsv
    fi

    # Remove PhiX
    msg "INFO: Removing PhiX using BBDuk"

    run_bbduk() {
      read1=$1
      read2=$2

      bbduk.sh \
        k=31 \
        hdist=1 \
        qout=33 \
        qin=auto \
        overwrite=t \
        in="${read1}" \
        in2="${read2}" \
        threads=!{task.cpus} \
        out=!{meta.id}_noPhiX-R1.fsq \
        out2=!{meta.id}_noPhiX-R2.fsq \
        ref="!{phix_reference_file}"

        echo $?
    }

    # Try reformatting reads if bbduk was unsuccessful
    if [ $(run_bbduk "!{reads[0]}" "!{reads[1]}") == 1 ]; then
      for read in !{reads}; do
        reformat.sh \
          in="${read}" \
          out="reformatted.${read}" \
          tossbrokenreads=t
      done

      # Run bbduk again on reformatted reads
      #  If this fails, input reads are corrupted and will automatically exit
      run_bbduk "reformatted.!{reads[0]}" "reformatted.!{reads[1]}"
    fi

    for suff in R1.fsq R2.fsq; do
      if verify_minimum_file_size "!{meta.id}_noPhiX-${suff}" 'PhiX-removed FastQ Files' "!{params.min_filesize_fastq_phix_removed}"; then
        echo -e "!{meta.id}\tPhiX-removed FastQ ($suff) File\tPASS" \
          >> !{meta.id}.PhiX-removed_FastQ_Files.tsv
      else
        echo -e "!{meta.id}\tPhiX-removed FastQ ($suff) File\tFAIL" \
          >> !{meta.id}.PhiX-removed_FastQ_Files.tsv
      fi
    done

    # Raw input read and bp information
    TOT_READS=$(grep '^Input: ' .command.err | awk '{print $2}')
    TOT_BASES=$(grep '^Input: ' .command.err | awk '{print $4}')

    if [[ -z "${TOT_READS}" || -z "${TOT_BASES}" ]]; then
      msg 'ERROR: unable to parse input counts from bbduk log' >&2
      exit 1
    fi

    # Number of PhiX read/bp contaminants
    NUM_PHIX_READS=$(grep '^Contaminants: ' .command.err | awk '{print $2}' | sed 's/,//g')
    PERCENT_PHIX_READS=$(grep '^Contaminants: ' .command.err | awk '{print $4}' | sed 's/[()]//g')
    NUM_PHIX_BASES=$(grep '^Contaminants: ' .command.err | awk '{print $5}' | sed 's/,//g')
    PERCENT_PHIX_BASES=$(grep '^Contaminants: ' .command.err | awk '{print $7}' | sed 's/[()]//g')

    # Cleaned FastQ file information
    NUM_CLEANED_READS=$(grep '^Result: ' .command.err | awk '{print $2}')
    PERCENT_CLEANED_READS=$(grep '^Result: ' .command.err | awk '{print $4}' | sed 's/[()]//g')
    NUM_CLEANED_BASES=$(grep '^Result: ' .command.err | awk '{print $4}')
    PERCENT_CLEANED_BASES=$(grep '^Result: ' .command.err | awk '{print $7}' | sed 's/[()]//g')

    msg "INFO: Input contains ${TOT_BASES} bp and $TOT_READS reads"
    msg "INFO: ${PHIX_BASES:-0} bp of PhiX were detected and ${PHIX_READS:-0} reads were removed"

    SUMMARY_HEADER="
      Sample name\t
      # Cleaned reads\t
      % Cleaned reads\t
      # Cleaned bp\t
      % Cleaned bp\t
      # PhiX reads\t
      % PhiX reads\t
      # PhiX Bp\t
      % PhiX bp\t
      # Raw reads\t
      # Raw bp
      "

    SUMMARY_OUTPUT="
      !{meta.id}\t
      ${NUM_CLEANED_READS}\t
      ${PERCENT_CLEANED_READS}\t
      ${NUM_CLEANED_BASES}\t
      ${PERCENT_CLEANED_BASES}\t
      ${NUM_PHIX_READS}\t
      ${PERCENT_PHIX_READS}\t
      ${NUM_PHIX_BASES}\t
      ${PERCENT_PHIX_BASES}\t
      ${TOT_READS}\t
      ${TOT_BASES}
      "

    echo -e $SUMMARY_HEADER > !{meta.id}.Summary.PhiX.tab
    echo -e $SUMMARY_OUTPUT >> !{meta.id}.Summary.PhiX.tab

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bbduk: $(bbduk.sh --version 2>&1 | head -n 2 | tail -1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
