process REMOVE_HOST_SRA_HUMAN_SCRUBBER {

    label "process_medium"
    tag { "${meta.id}" }
    container "quay.io/biocontainers/sra-human-scrubber@sha256:2f6b6635af9ba3190fc2f96640b21f0285483bd1f50d6be229228c52fb747055"

    input:
    tuple val(meta), path(reads)
    path(database)

    output:
    tuple val(meta), path("${meta.id}.SRA_Human_Scrubber_FastQ_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}_R*_scrubbed.fastq.gz")             , emit: host_removed_reads
    path("${meta.id}.SRA_Human_Scrubber_Removal.tsv")                    , emit: summary
    path("${meta.id}.SRA_Human_Scrubber_FastQ.SHA512-checksums.tsv")     , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                                 , emit: versions

    shell:
    // TODO: Use container built on Ubuntu to avoid filesize conversion (`find` is different here)
    minFilesize = params.min_filesize_fastq_sra_human_scrubber_removed
    min_filesize_output_fastq = ( ['c','b','k'].contains(minFilesize[-1])) ? "${minFilesize}" : "25000k"
    '''
    source bash_functions.sh

    # Use a non-default host to remove only if user-specified
    # TO-DO:  implementing this is much more difficult because it needs an
    #         unpackaged software to create the kmer .db file from NCBI's STAT
    #         https://github.com/ncbi/sra-human-scrubber/issues/20#issuecomment-1414392052

    # Remove Host Reads
    msg "INFO: Removing host reads from !{meta.id} using SRA Human Scrubber with !{task.memory} RAM ..."

    # NOTE: input only able to handle FastQ, not compressed
    # NOTE: output .gz filename doesn't compress, so requires gzip
    # NOTE: no handling for PE, do paired files 1-by-1
    # scrub the R1 FastQ file
    if [[ "!{reads[0]}" =~ .gz ]]; then
      msg "INFO: Removing host reads from GZ compressed R1 file "!{reads[0]}" ..."

      zcat "!{reads[0]}" | \
        scrub.sh \
        -d "!{database}" \
        -p "!{task.cpus}" \
        -x \
        -o "!{meta.id}_R1_scrubbed.fastq" \
        2> scrub_R1.stderr.txt

      gzip -f "!{meta.id}_R1_scrubbed.fastq"

    else
      msg "INFO: Removing host reads from uncompressed R1 file "!{reads[0]}" ..."

      scrub.sh \
        -i "!{reads[0]}" \
        -d "!{database}" \
        -p "!{task.cpus}" \
        -x \
        -o "!{meta.id}_R1_scrubbed.fastq" \
        2> scrub_R1.stderr.txt

      gzip -f "!{meta.id}_R1_scrubbed.fastq"
    fi
    msg "INFO: Completed host reads removal from R1 file of !{meta.id}"

    # Parse R1 counts input/output/removed
    R1_COUNT_READS_INPUT=$(grep 'total read count:' scrub_R1.stderr.txt \
      | awk 'BEGIN{FS=OFS="\t"}; {print $2}' | cut -d ':' -f 2 | sed 's/ //g')
    R1_COUNT_READS_REMOVED=$(grep 'spot(s) masked or removed.' scrub_R1.stderr.txt | awk '{print $1}')
    R1_COUNT_READS_OUTPUT=$(("${R1_COUNT_READS_INPUT}"-"${R1_COUNT_READS_REMOVED}"))

    # scrub the R2 FastQ file
    if [[ "!{reads[1]}" =~ .gz ]]; then
      msg "INFO: Removing host reads from GZ compressed R2 file "!{reads[1]}" ..."

      zcat "!{reads[1]}" | \
        scrub.sh \
        -d "!{database}" \
        -p "!{task.cpus}" \
        -x \
        -o "!{meta.id}_R2_scrubbed.fastq" \
        2> scrub_R2.stderr.txt

        gzip -f "!{meta.id}_R2_scrubbed.fastq"

    else
      msg "INFO: Removing host reads from uncompressed R2 file "!{reads[1]}" ..."

      scrub.sh \
        -i "!{reads[1]}" \
        -d "!{database}" \
        -p "!{task.cpus}" \
        -x \
        -o "!{meta.id}_R2_scrubbed.fastq" \
        2> scrub_R2.stderr.txt

      gzip -f "!{meta.id}_R2_scrubbed.fastq"
    fi
    msg "INFO: Completed host reads removal from R2 file of !{meta.id}"

    # Parse R2 counts input/output/removed
    R2_COUNT_READS_INPUT=$(grep 'total read count:' scrub_R2.stderr.txt \
      | awk 'BEGIN{FS=OFS="\t"}; {print $2}' | cut -d ':' -f 2 | sed 's/ //g')
    R2_COUNT_READS_REMOVED=$(grep 'spot(s) masked or removed.' scrub_R2.stderr.txt | awk '{print $1}')
    R2_COUNT_READS_OUTPUT=$(("${R2_COUNT_READS_INPUT}"-"${R2_COUNT_READS_REMOVED}"))

    # Validate output files are sufficient size to continue
    for file in "!{meta.id}_R1_scrubbed.fastq.gz" "!{meta.id}_R2_scrubbed.fastq.gz"; do
      if verify_minimum_file_size "${file}" 'SRA-Human-Scrubber-removed FastQ Files' "!{min_filesize_output_fastq}"; then
        echo -e "!{meta.id}\tSRA-Human-Scrubber-removed FastQ ($file) File\tPASS" \
          >> !{meta.id}.SRA_Human_Scrubber_FastQ_File.tsv
      else
        echo -e "!{meta.id}\tSRA-Human-Scrubber-removed FastQ ($file) File\tFAIL" \
          >> !{meta.id}.SRA_Human_Scrubber_FastQ_File.tsv
      fi
    done

    msg "INFO: Completed QC filesize checks of !{meta.id} for R1 and R2 after SRA Human Scrubbing"

    # Summarize R1 and R2 counts input/output/removed
    COUNT_READS_INPUT=$(("${R1_COUNT_READS_INPUT}"+"${R2_COUNT_READS_INPUT}"))
    COUNT_READS_REMOVED=$(("${R1_COUNT_READS_REMOVED}"+"${R2_COUNT_READS_REMOVED}"))
    COUNT_READS_OUTPUT=$(("${R1_COUNT_READS_OUTPUT}"+"${R2_COUNT_READS_OUTPUT}"))
    PERCENT_REMOVED=$(echo "${COUNT_READS_REMOVED}" "${COUNT_READS_INPUT}" \
      | awk '{proportion=$1/$2} END{printf("%.6f", proportion*100)}')
    PERCENT_OUTPUT=$(echo "${COUNT_READS_REMOVED}" "${COUNT_READS_INPUT}" \
      | awk '{proportion=$1/$2} END{printf("%.6f", 100-(proportion*100))}')

    # Ensure all values parsed properly from stderr output
    for val in "${COUNT_READS_INPUT}" "${COUNT_READS_OUTPUT}" "${COUNT_READS_REMOVED}"; do
      if [[ ! "${val}" =~ [0-9] ]]; then
        msg "ERROR: expected integer parsed from SRA Human Scrubber stderr instead of:${val}" >&2
        exit 1
      fi
    done
    for val in "${PERCENT_REMOVED}" "${PERCENT_OUTPUT}"; do
      if [[ ! "${val}" =~ [0-9.] ]]; then
          msg "ERROR: expected percentage parsed from SRA Human Scrubber stderr instead of:${val}" >&2
          exit 1
      fi
    done

    # Print read counts input/output from this process
    msg "INFO: Input contains ${COUNT_READS_INPUT} reads"
    msg "INFO: ${PERCENT_REMOVED}% of input reads were removed (${COUNT_READS_REMOVED} reads)"
    msg "INFO: ${COUNT_READS_OUTPUT} non-host reads (${PERCENT_OUTPUT}%) were retained"

    SUMMARY_HEADER=(
      "Sample_name"
      "Input_reads_(#)"
      "Removed_reads_(#)"
      "Removed_reads_(%)"
      "Output_reads_(#)"
      "Output_reads_(%)"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    SUMMARY_OUTPUT=(
      "!{meta.id}"
      "${COUNT_READS_INPUT}"
      "${COUNT_READS_REMOVED}"
      "${PERCENT_REMOVED}"
      "${COUNT_READS_OUTPUT}"
      "${PERCENT_OUTPUT}"
    )
    SUMMARY_OUTPUT=$(printf "%s\t" "${SUMMARY_OUTPUT[@]}" | sed 's/\t$//')

    # Store input/output counts
    echo -e "${SUMMARY_HEADER}" > "!{meta.id}.SRA_Human_Scrubber_Removal.tsv"
    echo -e "${SUMMARY_OUTPUT}" >> "!{meta.id}.SRA_Human_Scrubber_Removal.tsv"

    ### Calculate SHA-512 Checksums of each FastQ file ###
    SUMMARY_HEADER=(
      "Sample_name"
      "Checksum_(SHA-512)"
      "File"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.SRA_Human_Scrubber_FastQ.SHA512-checksums.tsv"

    # Calculate checksums
    for f in "!{meta.id}_R1_scrubbed.fastq.gz" "!{meta.id}_R2_scrubbed.fastq.gz"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.SRA_Human_Scrubber_FastQ.SHA512-checksums.tsv"
      BASE="$(basename ${f})"
      HASH=$(zcat "${f}" | awk 'NR%2==0' | paste - - | sort -k1,1 | sha512sum | awk '{print $1}')
      echo -e "${HASH}\t${BASE}"
    done >> "!{meta.id}.SRA_Human_Scrubber_FastQ.SHA512-checksums.tsv"

    # Get process version information
    # NOTE: currently no option to print the software version number, but
    #       track this issue https://github.com/ncbi/sra-human-scrubber/issues/28
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
        sra-human-scrubber: 2.2.1
    END_VERSIONS
    '''
}
