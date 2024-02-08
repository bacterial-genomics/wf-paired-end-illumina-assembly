process REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR {

    label "process_high"
    tag { "${meta.id}" }
    container "snads/bbtools@sha256:9f2a9b08563839cec87d856f0fc7607c235f464296fd71e15906ea1d15254695"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_repaired-R{1,2}.fastq.gz")              , emit: fastq_removed_broken_pairs
    tuple val(meta), path("${meta.id}.BBTools-Repair-removed_FastQ_Files.tsv"), emit: qc_filecheck
    path("${meta.id}.BBTools-Repair-Removal.tsv")                             , emit: summary
    path(".command.{out,err}")
    path("versions.yml")                                                      , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Remove broken sister read sequences
    msg "INFO: Removing broken sister reads using BBTools' Repair..."

    # NOTE: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/
    # "Repairing (repair flag) arbitrarily disordered files will take a lot of memory"
    # NOTE: "repair.sh requests all available memory by default"
    # NOTE: no CPU flag
    repair.sh \
      overwrite=t \
      in="!{reads[0]}" \
      in2="!{reads[1]}" \
      out=!{meta.id}_repaired-R1.fastq.gz \
      out2=!{meta.id}_repaired-R2.fastq.gz \
      outs=!{meta.id}_discarded_singletons.fastq \
      repair=t

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}.BBTools-Repair-removed_FastQ_Files.tsv"
    for suff in R1.fastq.gz R2.fastq.gz; do
      if verify_minimum_file_size "!{meta.id}_repaired-${suff}" 'Repaired FastQ Files' "!{params.min_filesize_broken_pairs_bbtools_repair_removed}"; then
        echo -e "!{meta.id}\tBBTools-repair-removed FastQ ($suff) File\tPASS" \
          >> "!{meta.id}.BBTools-Repair-removed_FastQ_Files.tsv"
      else
        echo -e "!{meta.id}\tBBTools-repair-removed FastQ ($suff) File\tFAIL" \
          >> "!{meta.id}.BBTools-Repair-removed_FastQ_Files.tsv"
      fi
    done

    # Raw input read and bp information
    INPUT_READS=$(grep '^Input: ' .command.err | awk '{print $2}')
    INPUT_BASES=$(grep '^Input: ' .command.err | awk '{print $4}')

    # Number of read/bp removed
    REMOVED_READS_COUNT=$(grep '^Singletons: ' .command.err | awk '{print $2}' | sed 's/,//g')
    REMOVED_READS_PERCENT=$(grep '^Singletons: ' .command.err | awk '{print $4}' | sed -e 's/[()]//g' -e 's/%//')
    REMOVED_BASES_COUNT=$(grep '^Singletons: ' .command.err | awk '{print $5}' | sed 's/,//g')
    REMOVED_BASES_PERCENT=$(grep '^Singletons: ' .command.err | awk '{print $7}' | sed -e 's/[()]//g' -e 's/%//')

    # Cleaned FastQ file information
    OUTPUT_READS_COUNT=$(grep '^Result: ' .command.err | awk '{print $2}')
    OUTPUT_READS_PERCENT=$(grep '^Result: ' .command.err | awk '{print $4}' | sed -e 's/[()]//g' -e 's/%//')
    OUTPUT_BASES_COUNT=$(grep '^Result: ' .command.err | awk '{print $5}')
    OUTPUT_BASES_PERCENT=$(grep '^Result: ' .command.err | awk '{print $7}' | sed -e 's/[()]//g' -e 's/%//')

    # Ensure all values parsed properly from stderr output
    for val in $INPUT_READS $INPUT_BASES $REMOVED_READS_COUNT $REMOVED_BASES_COUNT $OUTPUT_READS_COUNT $OUTPUT_BASES_COUNT; do
      if [[ ! "${val}" =~ [0-9] ]]; then
        msg "ERROR: expected integer parsed from bbtools repair stderr instead of:${val}" >&2
        exit 1
      fi
    done
    for val in $REMOVED_READS_PERCENT $REMOVED_BASES_PERCENT $OUTPUT_READS_PERCENT $OUTPUT_BASES_PERCENT; do
      if [[ ! "${val}" =~ [0-9.] ]]; then
          msg "ERROR: expected percentage parsed from SRA Human Scrubber stderr instead of:${val}" >&2
          exit 1
      fi
    done

    msg "INFO: Input contains ${INPUT_BASES} bp and ${INPUT_READS} reads"
    msg "INFO: ${REMOVED_READS_COUNT} (${REMOVED_READS_PERCENT}) reads were removed"
    msg "INFO: Output contains ${OUTPUT_BASES_COUNT} bp and ${OUTPUT_READS_COUNT} reads"

    DELIM=$'\t'
    SUMMARY_HEADER=(
      "Sample name"
      "# Input reads"
      "# Input bases"
      "# Output reads"
      "% Output reads"
      "# Output bases"
      "% Output bases"
      "# Removed reads"
      "% Removed reads"
      "# Removed bases"
      "% Removed bases"
    )
    SUMMARY_HEADER=$(printf "%s${DELIM}" "${SUMMARY_HEADER[@]}")
    SUMMARY_HEADER="${SUMMARY_HEADER%${DELIM}}"

    SUMMARY_OUTPUT=(
      "!{meta.id}"
      "${INPUT_READS}"
      "${INPUT_BASES}"
      "${OUTPUT_READS_COUNT}"
      "${OUTPUT_READS_PERCENT}"
      "${OUTPUT_BASES_COUNT}"
      "${OUTPUT_BASES_PERCENT}"
      "${REMOVED_READS_COUNT}"
      "${REMOVED_READS_PERCENT}"
      "${REMOVED_BASES_COUNT}"
      "${REMOVED_BASES_PERCENT}"
    )
    SUMMARY_OUTPUT=$(printf "%s${DELIM}" "${SUMMARY_OUTPUT[@]}")
    SUMMARY_OUTPUT="${SUMMARY_OUTPUT%${DELIM}}"

    # Store input/output counts
    echo -e "${SUMMARY_HEADER}" > !{meta.id}.BBTools-Repair-Removal.tsv
    echo -e "${SUMMARY_OUTPUT}" >> !{meta.id}.BBTools-Repair-Removal.tsv

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        repair.sh: $(repair.sh --version 2>&1 | head -n 2 | tail -1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
