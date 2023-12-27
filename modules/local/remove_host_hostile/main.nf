process REMOVE_HOST_HOSTILE {

    label "process_low"
    tag { "${meta.id}" }
    container "quay.io/biocontainers/hostile:0.2.0--pyhdfd78af_0@sha256:08b38e53f01f78877bdac9638263be0f5b6e02b709c302f56812532d448728b8"

    input:
    tuple val(meta), path(reads), path(qc_input_filecheck)

    output:
    path(".command.{out,err}")
    path "versions.yml"                                                                                    , emit: versions
    path "${meta.id}.Summary.Hostile-Removal.tsv"                                                          , emit: hostile_summary
    tuple val(meta), path("${meta.id}_noPhiX-R1.fsq"), path("${meta.id}_noPhiX-R2.fsq"), path("*File*.tsv"), emit: hostile_removed

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

    # Use a non-default host to remove only if user-specified
    HOST_INDEX_ARGUMENT=''
    if [[ ! -z "!{params.hostile_host_reference_path_prefix}" ]]; then
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
      fi
    fi

    # Remove Host Reads
    msg "INFO: Removing host reads using Hostile"

    hostile \
      clean \
      --fastq1 "!{reads[0]}" \
      --fastq2 "!{reads[0]}" \
      --out-dir hostile \
      "${HOST_INDEX_ARGUMENT}" \
      --threads !{task.cpus}

    # JSON format stdout reports input/output filenames and read counts
    if ! verify_minimum_file_size command.err 'JSON stdout for Hostile' '1k'; then
      msg "ERROR: JSON stdout missing or empty for reporting Hostile results" >&2
      exit 1
    fi

    # NOTE: grep used because `jq` absent from package
    RELATIVE_OUTPATH_R1=$(grep '"fastq1_out_path":' command.err | awk '{print $2}' | sed 's/[",]//g')
    RELATIVE_OUTPATH_R2=$(grep '"fastq2_out_path":' command.err | awk '{print $2}' | sed 's/[",]//g')

    # Validate output files are sufficient size to continue
    for file in ${RELATIVE_OUTPATH_R1} ${RELATIVE_OUTPATH_R2}; do
      if verify_minimum_file_size "${file}" 'Hostile-removed FastQ Files' "!{params.min_filesize_fastq_hostile_removed}"; then
        echo -e "!{meta.id}\tHostile-removed FastQ ($file) File\tPASS" \
          >> !{meta.id}.Hostile-removed_FastQ_Files.tsv
      else
        echo -e "!{meta.id}\tHostile-removed FastQ ($file) File\tFAIL" \
          >> !{meta.id}.Hostile-removed_FastQ_Files.tsv
      fi
    done

    # NOTE: grep used because `jq` absent from package
    COUNT_READS_INPUT=$(grep '"reads_in":' command.err | awk '{print $2}' | sed 's/,//g')
    COUNT_READS_OUTPUT=$(grep '"reads_out":' command.err | awk '{print $2}' | sed 's/,//g')
    COUNT_READS_REMOVED=$(grep '"reads_removed":' command.err | awk '{print $2}' | sed 's/,//g')
    COUNT_PERCENT_REMOVED=$(grep '"reads_removed_proportion":' command.err | awk '{$2=$2*100; print $2}')

    # Ensure all values parsed properly from JSON output report
    for val in $COUNT_READS_INPUT $COUNT_READS_OUTPUT $COUNT_READS_REMOVED; do
      if [[ ! "${val}" =~ [0-9] ]]; then
        msg "ERROR: expected integer parsed from Hostile JSON instead of:${val}" >&2
        exit 1
      fi
    done
    if [[ ! "${COUNT_PERCENT_REMOVED}" =~ [0-9.] ]]; then
        msg "ERROR: expected percentage parsed from Hostile JSON instead of:${COUNT_PERCENT_REMOVED}" >&2
        exit 1
    fi

    # Print read counts input/output from this process
    msg "INFO: Input contains ${COUNT_READS_INPUT} reads"
    msg "INFO: ${COUNT_PERCENT_REMOVED}% of input reads were removed (${COUNT_READS_REMOVED} reads)"
    msg "INFO: ${COUNT_READS_OUTPUT} non-host reads were retained"

    DELIM=$'\t'
    SUMMARY_HEADER=(
      "Sample name"
      "# Input reads"
      "# Output reads"
      "# Removed reads"
      "% Removed reads"
    )
    SUMMARY_HEADER=$(printf "%s${DELIM}" "${SUMMARY_HEADER[@]}")
    SUMMARY_HEADER="${SUMMARY_HEADER%${DELIM}}"

    SUMMARY_OUTPUT=(
      "!{meta.id}"
      "${COUNT_READS_INPUT}"
      "${COUNT_READS_OUTPUT}"
      "${COUNT_READS_REMOVED}"
      "${COUNT_PERCENT_REMOVED}"
    )
    SUMMARY_OUTPUT=$(printf "%s${DELIM}" "${SUMMARY_OUTPUT[@]}"
    SUMMARY_OUTPUT="${SUMMARY_OUTPUT%${DELIM}}"

    # Store input/output counts
    echo -e "${SUMMARY_HEADER}" > !{meta.id}.Summary.Hostile-Removal.tsv
    echo -e "${SUMMARY_OUTPUT}" >> !{meta.id}.Summary.Hostile-Removal.tsv

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        hostile: $(hostile --version)
    END_VERSIONS
    '''
}
