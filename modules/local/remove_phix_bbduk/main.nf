process REMOVE_PHIX_BBDUK {

    label "process_low"
    tag { "${meta.id}" }
    container "snads/bbtools@sha256:9f2a9b08563839cec87d856f0fc7607c235f464296fd71e15906ea1d15254695"

    input:
    tuple val(meta), path(reads)
    path phix_reference_file

    output:
    tuple val(meta), path("${meta.id}_noPhiX-R{1,2}.fsq"), emit: fastq_phix_removed
    tuple val(meta), path("${meta.id}.PhiX*_File.tsv")   , emit: qc_filecheck
    path("${meta.id}.Summary.PhiX.tsv")                  , emit: phix_summary
    path(".command.{out,err}")
    path("versions.yml")                                 , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Verify PhiX reference file size
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}.PhiX_Genome_File.tsv"
    if verify_minimum_file_size !{phix_reference_file} 'PhiX Genome' "!{params.min_filesize_phix_genome}"; then
      echo -e "!{meta.id}\tPhiX Genome\tPASS" >> "!{meta.id}.PhiX_Genome_File.tsv"
    else
      echo -e "!{meta.id}\tPhiX Genome\tFAIL" >> "!{meta.id}.PhiX_Genome_File.tsv"
    fi

    # Auto reformat FastQ files
    # msg "INFO: Auto reformatting FastQ files.."
    # for read in !{reads}; do
    #   reformat.sh \
    #     in="${read}" \
    #     out="reformatted.${read}" \
    #     tossbrokenreads=t
    # done

    # Remove PhiX
    msg "INFO: Removing PhiX using BBDuk.."

    bbduk.sh \
      k=31 \
      hdist=1 \
      qout=33 \
      qin=auto \
      overwrite=t \
      in="!{reads[0]}" \
      in2="!{reads[1]}" \
      threads=!{task.cpus} \
      out=!{meta.id}_noPhiX-R1.fsq \
      out2=!{meta.id}_noPhiX-R2.fsq \
      ref="!{phix_reference_file}"

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}.PhiX-removed_FastQ_File.tsv"
    for suff in R1.fsq R2.fsq; do
      if verify_minimum_file_size "!{meta.id}_noPhiX-${suff}" 'PhiX-removed FastQ Files' "!{params.min_filesize_fastq_phix_removed}"; then
        echo -e "!{meta.id}\tPhiX-removed FastQ ($suff) File\tPASS" \
          >> "!{meta.id}.PhiX-removed_FastQ_File.tsv"
      else
        echo -e "!{meta.id}\tPhiX-removed FastQ ($suff) File\tFAIL" \
          >> "!{meta.id}.PhiX-removed_FastQ_File.tsv"
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
    NUM_CLEANED_BASES=$(grep '^Result: ' .command.err | awk '{print $5}')
    PERCENT_CLEANED_BASES=$(grep '^Result: ' .command.err | awk '{print $7}' | sed 's/[()]//g')

    msg "INFO: Input contains ${TOT_BASES} bp and $TOT_READS reads"
    msg "INFO: ${PHIX_BASES:-0} bp of PhiX were detected and ${PHIX_READS:-0} reads were removed"

    SUMMARY_HEADER=(
      "Sample name",
      "# Cleaned reads",
      "% Cleaned reads",
      "# Cleaned bp",
      "% Cleaned bp",
      "# PhiX reads",
      "% PhiX reads",
      "# PhiX Bp",
      "% PhiX bp",
      "# Raw reads",
      "# Raw bp"
    )

    SUMMARY_OUTPUT=(
      "!{meta.id}",
      "${NUM_CLEANED_READS}",
      "${PERCENT_CLEANED_READS}",
      "${NUM_CLEANED_BASES}",
      "${PERCENT_CLEANED_BASES}",
      "${NUM_PHIX_READS}",
      "${PERCENT_PHIX_READS}",
      "${NUM_PHIX_BASES}",
      "${PERCENT_PHIX_BASES}",
      "${TOT_READS}",
      "${TOT_BASES}",
      )

    printf "%s\n" "${SUMMARY_HEADER[@]}" | tr ',' '\t' > "!{meta.id}.Summary.PhiX.tsv"
    printf "%s" "${SUMMARY_OUTPUT[@]}" | tr ',' '\t' >> "!{meta.id}.Summary.PhiX.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bbduk: $(bbduk.sh --version 2>&1 | head -n 2 | tail -1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
