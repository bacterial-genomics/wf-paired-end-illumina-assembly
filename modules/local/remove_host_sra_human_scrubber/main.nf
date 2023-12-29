process REMOVE_HOST_SRA_HUMAN_SCRUBBER {

    label "process_low"
    tag { "${meta.id}" }
    container "ncbi/sra-human-scrubber@sha256:96723dba29c70dcda1e26663c92675c55f9a2340b592db7b1f7013a8b8247fe4"

    input:
    tuple val(meta), path(reads)

    output:
    path(".command.{out,err}")
    path "versions.yml"                                     , emit: versions
    path "${meta.id}.Summary.SRA-Human-Scrubber-Removal.tsv", emit: sra_human_scrubber_summary
    tuple val(meta), path("${meta.id}_R*_scrubbed.fastq.gz"), emit: sra_human_scrubber_removed

    shell:
    '''
    source bash_functions.sh

    # Use a non-default host to remove only if user-specified
    # TO-DO:  implementing this is much more difficult because it needs an
    #         unpackaged software to create the kmer .db file from NCBI's STAT
    #         https://github.com/ncbi/sra-human-scrubber/issues/20#issuecomment-1414392052

    # Remove Host Reads
    msg "INFO: Removing host reads using SRA Human Scrubber"

    # NOTE: input only able to handle FastQ, not compressed
    # NOTE: output .gz filename doesn't compress, so requires gzip
    # NOTE: no handling for PE, do paired files 1-by-1
    # scrub the R1 FastQ file
    OUTFILE_SCRUBBED_R1="!{meta.id}"_R1_scrubbed.fastq.gz
    if [[ "!{reads[0]}" =~ .gz ]]; then
      zcat "!{reads[0]}" | \
        scrub.sh \
        -d /scicomp/groups-pure/OID/NCEZID/DHCPP/BSPB/ZSAL/.databases/sra-human-scrubber/data/human_filter.db \
        -p !{task.cpus} | \
        gzip > "${OUTFILE_SCRUBBED_R1}"
    else
      scrub.sh \
        -i "!{reads[0]}" \
        -p !{task.cpus} | \
        gzip > "${OUTFILE_SCRUBBED_R1}"
    fi

    # Parse R1 counts input/output/removed
    R1_COUNT_READS_INPUT=$(grep 'total read count:' command.err \
     | awk 'BEGIN{FS=OFS="\t"}; {print $2}' | cut -d ':' -f 2 | sed 's/ //g')
    R1_COUNT_READS_REMOVED=$(grep 'spot(s) masked or removed.' command.err | awk '{print $1}')
    R1_COUNT_READS_OUTPUT=$(("${R1_COUNT_READS_INPUT}"-"${R1_COUNT_READS_REMOVED}"))

    # scrub the R2 FastQ file
    OUTFILE_SCRUBBED_R2="!{meta.id}"_R2_scrubbed.fastq.gz
    if [[ "!{reads[1]}" =~ .gz ]]; then
      zcat "!{reads[1]}" | \
        scrub.sh \
        -d /scicomp/groups-pure/OID/NCEZID/DHCPP/BSPB/ZSAL/.databases/sra-human-scrubber/data/human_filter.db \
        -p !{task.cpus} | \
        gzip > "${OUTFILE_SCRUBBED_R2}"
    else
      scrub.sh \
        -i "!{reads[1]}" \
        -p !{task.cpus} | \
        gzip > "${OUTFILE_SCRUBBED_R2}"
    fi

    # Parse R2 counts input/output/removed
    R2_COUNT_READS_INPUT=$(grep 'total read count:' command.err \
     | awk 'BEGIN{FS=OFS="\t"}; {print $2}' | cut -d ':' -f 2 | sed 's/ //g')
    R2_COUNT_READS_REMOVED=$(grep 'spot(s) masked or removed.' command.err | awk '{print $1}')
    R2_COUNT_READS_OUTPUT=$(("${R2_COUNT_READS_INPUT}"-"${R2_COUNT_READS_REMOVED}"))

    # Validate output files are sufficient size to continue
    for file in ${OUTFILE_SCRUBBED_R1} ${OUTFILE_SCRUBBED_R2}; do
      if verify_minimum_file_size "${file}" 'SRA-Human-Scrubber-removed FastQ Files' "!{params.min_filesize_fastq_sra_human_scrubber_removed}"; then
        echo -e "!{meta.id}\tSRA-Human-Scrubber-removed FastQ ($file) File\tPASS" \
          >> !{meta.id}.SRA_Human_Scrubber_FastQ_File.tsv
      else
        echo -e "!{meta.id}\tSRA-Human-Scrubber-removed FastQ ($file) File\tFAIL" \
          >> !{meta.id}.SRA_Human_Scrubber_FastQ_File.tsv
      fi
    done

    # Summarize R1 and R2 counts input/output/removed
    COUNT_READS_INPUT=$(("${R1_COUNT_READS_INPUT}"+"${R2_COUNT_READS_INPUT}"))
    COUNT_READS_REMOVED=$(("${R1_COUNT_READS_REMOVED}"+"${R2_COUNT_READS_REMOVED}"))
    COUNT_READS_OUTPUT=$(("${R1_COUNT_READS_OUTPUT}"+"${R2_COUNT_READS_OUTPUT}"))
    PERCENT_REMOVED=$(echo "${COUNT_READS_REMOVED}" "${COUNT_READS_INPUT}" \
     | awk '{proportion=$1/$2} END{printf("%.6f", proportion*100)}')
    PERCENT_OUTPUT=$(echo "${COUNT_READS_REMOVED}" "${COUNT_READS_INPUT}" \
     | awk '{proportion=$1/$2} END{printf("%.6f", 100-(proportion*100))}')

    # Ensure all values parsed properly from stderr output
    for val in $COUNT_READS_INPUT $COUNT_READS_OUTPUT $COUNT_READS_REMOVED; do
      if [[ ! "${val}" =~ [0-9] ]]; then
        msg "ERROR: expected integer parsed from SRA Human Scrubber stderr instead of:${val}" >&2
        exit 1
      fi
    done
    for val in $PERCENT_REMOVED $PERCENT_OUTPUT; do
      if [[ ! "${val}" =~ [0-9.] ]]; then
          msg "ERROR: expected percentage parsed from SRA Human Scrubber stderr instead of:${val}" >&2
          exit 1
      fi
    done

    # Print read counts input/output from this process
    msg "INFO: Input contains ${COUNT_READS_INPUT} reads"
    msg "INFO: ${PERCENT_REMOVED}% of input reads were removed (${COUNT_READS_REMOVED} reads)"
    msg "INFO: ${COUNT_READS_OUTPUT} non-host reads (${PERCENT_OUTPUT}%) were retained"

    DELIM=$'\t'
    SUMMARY_HEADER=(
      "Sample name"
      "# Input reads"
      "# Output reads"
      "% Output reads"
      "# Removed reads"
      "% Removed reads"
    )
    SUMMARY_HEADER=$(printf "%s${DELIM}" "${SUMMARY_HEADER[@]}")
    SUMMARY_HEADER="${SUMMARY_HEADER%${DELIM}}"

    SUMMARY_OUTPUT=(
      "!{meta.id}"
      "${COUNT_READS_INPUT}"
      "${COUNT_READS_OUTPUT}"
      "${PERCENT_OUTPUT}"
      "${COUNT_READS_REMOVED}"
      "${PERCENT_REMOVED}"
    )
    SUMMARY_OUTPUT=$(printf "%s${DELIM}" "${SUMMARY_OUTPUT[@]}")
    SUMMARY_OUTPUT="${SUMMARY_OUTPUT%${DELIM}}"

    # Store input/output counts
    echo -e "${SUMMARY_HEADER}" > !{meta.id}.Summary.SRA-Human-Scrubber-Removal.tsv
    echo -e "${SUMMARY_OUTPUT}" >> !{meta.id}.Summary.SRA-Human-Scrubber-Removal.tsv

    # Get process version information
    # NOTE: currently no option to print the software version number, but
    #       track this issue https://github.com/ncbi/sra-human-scrubber/issues/28 
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sra-human-scrubber: 2.2.1
    END_VERSIONS
    '''
}
