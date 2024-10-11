process REMOVE_PHIX_BBDUK {

    label "process_high"
    tag { "${meta.id}" }
    container "snads/bbtools@sha256:9f2a9b08563839cec87d856f0fc7607c235f464296fd71e15906ea1d15254695"
    // NOTE: "staphb/bbtools@sha256-f7b98063910e2e3b5be12f62076ec5cfdeaa562a01596758feb9a892ce18a363"
    // somtimes gives "bbtools bgzip: error while loading shared libraries: libcurl-gnutls.so.4: cannot open shared object file"
    // error with some samples (e.g., SRR14718846). Need to upgrade but find out what we're doing different.
    // Dockerfile or cmd difference issue? or both?

    input:
    tuple val(meta), path(reads)
    path phix_reference_file

    output:
    tuple val(meta), path("${meta.id}_noPhiX-R{1,2}.fsq")   , emit: phix_removed_reads
    tuple val(meta), path("${meta.id}.PhiX_Genome_File.tsv"), emit: qc_filecheck
    path("${meta.id}.PhiX.tsv")                             , emit: summary
    path("${meta.id}.noPhiX_FastQ.SHA256-checksums.tsv")    , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                    , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Verify PhiX reference file size
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.PhiX_Genome_File.tsv"
    if verify_minimum_file_size !{phix_reference_file} 'PhiX Genome' "!{params.min_filesize_phix_genome}"; then
      echo -e "!{meta.id}\tPhiX Genome FastA File\tPASS" >> "!{meta.id}.PhiX_Genome_File.tsv"
    else
      echo -e "!{meta.id}\tPhiX Genome FastA File\tFAIL" >> "!{meta.id}.PhiX_Genome_File.tsv"
    fi

    # Remove PhiX
    msg "INFO: Removing PhiX from !{meta.id} using BBDuk..."

    # NOTE: With excess sequence reads, it is very normal and possible to see initial error of
    #       "NOTE: Process `ASSEMBLY:REMOVE_PHIX_BBDUK (name)` terminated with an error exit status (140) -- Execution is retried (1)"
    #       But an automatic retry in the workflow with increase RAM should process the bulky sample just fine.
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

    msg "INFO: PhiX removed from !{meta.id} using BBDuk"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.PhiX-removed_FastQ_File.tsv"
    for suff in R1.fsq R2.fsq; do
      if verify_minimum_file_size "!{meta.id}_noPhiX-${suff}" 'PhiX-removed FastQ Files' "!{params.min_filesize_fastq_phix_removed}"; then
        echo -e "!{meta.id}\tPhiX-removed (${suff}) FastQ File\tPASS" \
          >> "!{meta.id}.PhiX-removed_FastQ_File.tsv"
      else
        echo -e "!{meta.id}\tPhiX-removed (${suff}) FastQ File\tFAIL" \
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

    msg "INFO: Input contains ${TOT_BASES} bp and ${TOT_READS} reads"
    msg "INFO: ${PHIX_BASES:-0} bp of PhiX were detected and ${PHIX_READS:-0} reads were removed"

    SUMMARY_HEADER=(
      "Sample_name"
      "Cleaned_reads_(#)"
      "Cleaned_reads_(%)"
      "Cleaned_basepairs_(#)"
      "Cleaned_basepairs_(%)"
      "PhiX_reads_(#)"
      "PhiX_reads_(%)"
      "PhiX_basepairs_(#)"
      "PhiX_basepairs_(%)"
      "Raw_reads_(#)"
      "Raw_basepairs_(#)"
    )

    SUMMARY_OUTPUT=(
      "!{meta.id}"
      "${NUM_CLEANED_READS}"
      "${PERCENT_CLEANED_READS}"
      "${NUM_CLEANED_BASES}"
      "${PERCENT_CLEANED_BASES}"
      "${NUM_PHIX_READS}"
      "${PERCENT_PHIX_READS}"
      "${NUM_PHIX_BASES}"
      "${PERCENT_PHIX_BASES}"
      "${TOT_READS}"
      "${TOT_BASES}"
    )

    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//1')
    SUMMARY_OUTPUT=$(printf "%s\t" "${SUMMARY_OUTPUT[@]}" | sed 's/\t$//1')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.PhiX.tsv"
    echo "${SUMMARY_OUTPUT}" >> "!{meta.id}.PhiX.tsv"

    ### Calculate SHA-256 Checksums of each FastQ file ###
    msg "INFO: Calculating checksums for !{meta.id}_noPhiX-R1.fsq !{meta.id}_noPhiX-R2.fsq ..."

    SUMMARY_HEADER=(
      "Sample_name"
      "Checksum_(SHA-256)"
      "File"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//1')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.noPhiX_FastQ.SHA256-checksums.tsv"

    # Calculate checksums
    for f in "!{meta.id}_noPhiX-R1.fsq" "!{meta.id}_noPhiX-R2.fsq"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.noPhiX_FastQ.SHA256-checksums.tsv"
      awk 'NR%2==0'  "${f}" | paste - - | sort -k1,1 | sha256sum | awk '{print $1 "\t" "'"$f"'"}'
    done >> "!{meta.id}.noPhiX_FastQ.SHA256-checksums.tsv"

    msg "INFO: Calculated checksums for !{meta.id}_noPhiX-R1.fsq !{meta.id}_noPhiX-R2.fsq"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        find: $(find --version | grep "^find" | sed 's/find //1')
        sha256sum: $(sha256sum --version | grep "^sha256sum" | sed 's/sha256sum //1')
        bbduk: $(bbduk.sh --version 2>&1 | grep "^BBMap" | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
