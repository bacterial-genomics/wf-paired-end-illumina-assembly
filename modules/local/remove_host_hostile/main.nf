process REMOVE_HOST_HOSTILE {

    label "process_high"
    tag { "${meta.id}" }
    container "quay.io/biocontainers/hostile:0.2.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.Hostile-removed_FastQ_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("hostile/${meta.id}*.clean_*")              , emit: host_removed_reads
    path("${meta.id}.Hostile_Removal.tsv")                            , emit: summary
    path("${meta.id}.Hostile_FastQ.SHA512-checksums.tsv")             , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                              , emit: versions

    shell:
    // TODO: Use container built on Ubuntu
    minFilesize = params.min_filesize_fastq_hostile_removed
    min_filesize_output_fastq = ( ['c','b','k'].contains(minFilesize[-1])) ? "${minFilesize}" : "25000k"
    '''
    source bash_functions.sh

    # Use a non-default host to remove only if user-specified
    HOST_INDEX_ARGUMENT=''
    if [[ ! -z "!{params.hostile_host_reference_path_prefix}" ]]; then

      msg "INFO: Detected user-specified hostile reference path. Validating its Bowtie2 index files ..."

      for suff in .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2; do
        if verify_minimum_file_size "!{params.hostile_host_reference_path_prefix}${suff}" 'Bowtie2 Index file for Hostile' '1c'; then
          continue
        else
          msg "ERROR: !{params.hostile_host_reference_path_prefix}${suff} absent; required for Illumina host removal by Hostile" >&2
          exit 1
        fi
      done

      bowtie2-inspect \
        --summary \
        !{params.hostile_host_reference_path_prefix} \
        1> !{meta.id}.bowtie2-inspect.stdout.log \
        2> !{meta.id}.bowtie2-inspect.stderr.log
      HOST_INDEX_ARGUMENT="--index !{params.hostile_host_reference_path_prefix}"

      msg "INFO: Completed validation of user-specified hostile Bowtie2 index files"

    fi

    # Remove Host Reads
    msg "INFO: Removing host reads using Hostile for !{meta.id} ..."

    if [[ ! -z "!{params.hostile_host_reference_path_prefix}" ]]; then
      hostile \
        clean \
        --fastq1 "!{reads[0]}" \
        --fastq2 "!{reads[1]}" \
        --out-dir hostile \
        --force \
        "${HOST_INDEX_ARGUMENT}" \
        --threads !{task.cpus}
    else
      hostile \
        clean \
        --fastq1 "!{reads[0]}" \
        --fastq2 "!{reads[1]}" \
        --out-dir hostile \
        --force \
        --threads !{task.cpus}
    fi

    msg "INFO: Completed host removal with Hostile for !{meta.id}"

    # JSON format stdout reports input/output filenames and read counts
    if ! verify_minimum_file_size .command.out 'JSON stdout for Hostile' '300c'; then
      msg "ERROR: JSON stdout missing or empty for reporting Hostile results" >&2
      exit 1
    fi

    msg "INFO: Parsing host removal output files from Hostile for !{meta.id} ..."

    # NOTE: grep used because `jq` absent from package
    RELATIVE_OUTPATH_R1=$(grep '"fastq1_out_path":' .command.out | awk '{print $2}' | sed 's/[",]//g')
    RELATIVE_OUTPATH_R2=$(grep '"fastq2_out_path":' .command.out | awk '{print $2}' | sed 's/[",]//g')

    msg "INFO: Completed parsing host removal output files from Hostile for !{meta.id}"

    # Validate output files are sufficient size to continue
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.Hostile-removed_FastQ_File.tsv"

    i=1
    for file in ${RELATIVE_OUTPATH_R1} ${RELATIVE_OUTPATH_R2}; do
      if verify_minimum_file_size "${file}" 'Hostile-removed FastQ Files' "!{min_filesize_output_fastq}"; then
        echo -e "!{meta.id}\tHostile-removed FastQ (R${i}) File\tPASS" \
          >> !{meta.id}.Hostile-removed_FastQ_File.tsv
      else
        echo -e "!{meta.id}\tHostile-removed FastQ (R${i}) File\tFAIL" \
          >> !{meta.id}.Hostile-removed_FastQ_File.tsv
      fi
      ((i++))
    done

    msg "INFO: Parsing host removal metrics from Hostile for !{meta.id} ..."

    # NOTE: grep used because `jq` absent from package
    COUNT_READS_INPUT=$(grep '"reads_in":' .command.out | awk '{print $2}' | sed 's/,//g')
    COUNT_READS_OUTPUT=$(grep '"reads_out":' .command.out | awk '{print $2}' | sed 's/,//g')
    PERCENT_OUTPUT=$(grep '"reads_removed_proportion":' .command.out | awk '{$2=$2*100; print 100-$2}')
    COUNT_READS_REMOVED=$(grep '"reads_removed":' .command.out | awk '{print $2}' | sed 's/,//g')
    PERCENT_REMOVED=$(grep '"reads_removed_proportion":' .command.out | awk '{$2=$2*100; print $2}')

    msg "INFO: Completed parsing host removal metrics from Hostile for !{meta.id}"

    # Ensure all values parsed properly from JSON output report
    for val in "${COUNT_READS_INPUT}" "${COUNT_READS_OUTPUT}" "${COUNT_READS_REMOVED}"; do
      if [[ ! "${val}" =~ [0-9] ]]; then
        msg "ERROR: expected integer parsed from Hostile JSON instead of:${val}" >&2
        exit 1
      fi
    done
    for val in "${PERCENT_REMOVED}" "${PERCENT_OUTPUT}"; do
      if [[ ! "${val}" =~ [0-9.] ]]; then
          msg "ERROR: expected percentage parsed from Hostile JSON instead of:${val}" >&2
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
    echo -e "${SUMMARY_HEADER}" > !{meta.id}.Hostile_Removal.tsv
    echo -e "${SUMMARY_OUTPUT}" >> !{meta.id}.Hostile_Removal.tsv

    ### Calculate SHA-512 Checksums of each FastQ file ###
    SUMMARY_HEADER=(
      "Sample_name"
      "Checksum_(SHA-512)"
      "File"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.Hostile_FastQ.SHA512-checksums.tsv"

    # Calculate checksums
    for f in "${RELATIVE_OUTPATH_R1}" "${RELATIVE_OUTPATH_R2}"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.Hostile_FastQ.SHA512-checksums.tsv"
      BASE="$(basename ${f})"
      HASH=$(zcat "${f}" | awk 'NR%2==0' | paste - - | sort -k1,1 | sha512sum | awk '{print $1}')
      echo -e "${HASH}\t${BASE}"
    done >> "!{meta.id}.Hostile_FastQ.SHA512-checksums.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
        hostile: $(hostile --version)
    END_VERSIONS
    '''
}
