process REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR {

    label "process_high"
    tag { "${meta.id}" }
    container "staphb/bbtools@sha256:f7b98063910e2e3b5be12f62076ec5cfdeaa562a01596758feb9a892ce18a363"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.BBTools-Repair-removed_FastQ_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}_repaired-R{1,2}.fastq.gz")             , emit: repaired_reads
    path("${meta.id}.BBTools_Repair_Removal.tsv")                            , emit: summary
    path("${meta.id}.BBTools_Repair_Removed_FastQ.SHA512-checksums.tsv")     , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                                     , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Remove broken sister read sequences
    msg "INFO: Removing broken sister reads for !{meta.id} using BBTools' Repair with !{task.memory} RAM ..."

    # NOTE: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/
    # "Repairing (repair flag) arbitrarily disordered files will take a lot of memory"
    # NOTE: "repair.sh requests all available memory by default"
    # NOTE: no CPU flag
    repair.sh \
      overwrite=t \
      in=!{reads[0]} \
      in2=!{reads[1]} \
      out=!{meta.id}_repaired-R1.fastq \
      out2=!{meta.id}_repaired-R2.fastq \
      outs=!{meta.id}_discarded_singletons.fastq \
      repair=t

    msg "INFO: Completed removal of broken sister reads for !{meta.id} using BBTools' Repair"

    # NOTE: repair.sh handles .gz outfile extension but it can get stuck hanging
    #       when there's errors like:
    #       "bgzip: error while loading shared libraries: libcurl-gnutls.so.4: cannot open shared object file: No such file or directory"
    #       "Caused by: java.io.IOException: Stream closed"
    gzip -f !{meta.id}_repaired-R1.fastq !{meta.id}_repaired-R2.fastq

    msg "INFO: Completed FastQ compression of repaired sister reads for !{meta.id}"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.BBTools-Repair-removed_FastQ_File.tsv"
    for suff in R1.fastq.gz R2.fastq.gz; do
      if verify_minimum_file_size "!{meta.id}_repaired-${suff}" 'Repaired FastQ Files' "!{params.min_filesize_broken_pairs_bbtools_repair_removed}"; then
        echo -e "!{meta.id}\tBBTools-repair-removed FastQ ($suff) File\tPASS" \
          >> "!{meta.id}.BBTools-Repair-removed_FastQ_File.tsv"
      else
        echo -e "!{meta.id}\tBBTools-repair-removed FastQ ($suff) File\tFAIL" \
          >> "!{meta.id}.BBTools-Repair-removed_FastQ_File.tsv"
          exit 1
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
    for val in "${INPUT_READS}" "${INPUT_BASES}" "${REMOVED_READS_COUNT}" "${REMOVED_BASES_COUNT}" "${OUTPUT_READS_COUNT}" "${OUTPUT_BASES_COUNT}"; do
      if [[ ! "${val}" =~ [0-9] ]]; then
        msg "ERROR: expected integer parsed from bbtools repair stderr instead of:${val}" >&2
        exit 1
      fi
    done
    for val in "${REMOVED_READS_PERCENT}" "${REMOVED_BASES_PERCENT}" "${OUTPUT_READS_PERCENT}" "${OUTPUT_BASES_PERCENT}"; do
      if [[ ! "${val}" =~ [0-9.] ]]; then
        msg "ERROR: expected percentage parsed from SRA Human Scrubber stderr instead of:${val}" >&2
        exit 1
      fi
    done

    msg "INFO: Input contains ${INPUT_BASES} bp and ${INPUT_READS} reads"
    msg "INFO: ${REMOVED_READS_COUNT} (${REMOVED_READS_PERCENT}) reads were removed"
    msg "INFO: Output contains ${OUTPUT_BASES_COUNT} bp and ${OUTPUT_READS_COUNT} reads"

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
      "${INPUT_READS}"
      "${REMOVED_READS_COUNT}"
      "${REMOVED_READS_PERCENT}"
      "${OUTPUT_READS_COUNT}"
      "${OUTPUT_READS_PERCENT}"
    )
    SUMMARY_OUTPUT=$(printf "%s\t" "${SUMMARY_OUTPUT[@]}" | sed 's/\t$//')

    # Store input/output counts
    echo -e "${SUMMARY_HEADER}" > "!{meta.id}.BBTools_Repair_Removal.tsv"
    echo -e "${SUMMARY_OUTPUT}" >> "!{meta.id}.BBTools_Repair_Removal.tsv"

    # Calculate checksums
    for f in "!{meta.id}_repaired-R1.fastq.gz" "!{meta.id}_repaired-R2.fastq.gz"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.BBTools_Repair_Removed_FastQ.SHA512-checksums.tsv"
      BASE="$(basename ${f})"
      HASH=$(zcat "${f}" | awk 'NR%2==0' | paste - - | sort -k1,1 | sha512sum | awk '{print $1}')
      echo -e "${HASH}\t${BASE}"
    done >> "!{meta.id}.BBTools_Repair_Removed_FastQ.SHA512-checksums.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
        repair.sh: $(repair.sh --version 2>&1 | head -n 2 | tail -1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
