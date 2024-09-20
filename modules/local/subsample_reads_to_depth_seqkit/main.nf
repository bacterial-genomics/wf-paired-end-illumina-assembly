process SUBSAMPLE_READS_TO_DEPTH_SEQKIT {

    label "process_medium"
    tag { "${meta.id}" }
    container "staphb/seqkit@sha256:8eb09a52ae932f7c25cfbb8db0df7110567087a187c7e90d46f499962d1c82c9"

    input:
    tuple val(meta), path(reads), path(depth), path(fraction_of_reads)

    output:
    tuple val(meta), path("*.{fastq,fq}.gz", includeInputs: true), emit: reads
    path("${meta.id}.Subsampled_FastQ.SHA512-checksums.tsv")     , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                         , emit: versions

    shell:
    seqkit_seed = (params.seqkit_seed >= 1)? params.seqkit_seed : 947266746

    '''
    source bash_functions.sh

    fraction_of_reads_to_use=$(cat !{fraction_of_reads})
    initial_depth=$(cat !{depth})

    depth="!{params.depth}"

    echo "!{params.seqkit_seed}" > seed-value.txt

    if ! [[ ${fraction_of_reads_to_use} =~ ^[0-9.]+$ ]]; then
      msg "ERROR: Unable to calculate a fraction to subsample; ${fraction_of_reads_to_use} not a floating point value" >&2
      exit 1
    fi
    if [ ${depth%.*} -gt 0 ] && [ ${initial_depth%.*} -gt ${depth%.*} ]; then
      msg "INFO: Subsampling !{meta.id} R1 with seqkit using seed:!{params.seqkit_seed} ..."

      seqkit sample \
        !{reads[0]} \
        --threads !{task.cpus} \
        --proportion ${fraction_of_reads_to_use} \
        --rand-seed "!{params.seqkit_seed}" \
        --out-file "!{meta.id}_R1.subsampled.fastq.gz" \
        2> seqkit.R1.stderr.txt

      msg "INFO: Subsampling !{meta.id} R2 with seqkit using seed:!{params.seqkit_seed} ..."

      seqkit sample \
        !{reads[1]} \
        --threads !{task.cpus} \
        --proportion ${fraction_of_reads_to_use} \
        --rand-seed "!{params.seqkit_seed}" \
        --out-file "!{meta.id}_R2.subsampled.fastq.gz" \
        2> seqkit.R2.stderr.txt

      msg "INFO: Completed subsampling of R1 and R2 by seqkit"

      # Discard symlink infiles to avoid them being passed as outfiles when
      #   subsampling occurred.
      rm -f !{reads[0]} !{reads[1]}

      number_output_R1_sequences=$(grep 'sequences outputted' seqkit.R1.stderr.txt | awk '{print $2}')
      number_output_R2_sequences=$(grep 'sequences outputted' seqkit.R2.stderr.txt | awk '{print $2}')

      msg "INFO: Subsampled reads contain: ${number_output_R1_sequences} and ${number_output_R2_sequences} sequences"

    else
      # The input FastQ files that were never subsampled will get passed on
      #   as outputs here with the 'includeInputs: true'
      msg "INFO: Subsampling not requested or required for !{meta.id}"
      touch "!{meta.id}.Subsampled_FastQ.SHA512-checksums.tsv" versions.yml
      exit 0
    fi

    ### Calculate SHA-512 Checksums of each FastQ file ###
    SUMMARY_HEADER=(
      "Sample_name"
      "Checksum_(SHA-512)"
      "File"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.Subsampled_FastQ.SHA512-checksums.tsv"

    # Calculate checksums
    for f in "!{meta.id}_R1.subsampled.fastq.gz" "!{meta.id}_R2.subsampled.fastq.gz"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.Subsampled_FastQ.SHA512-checksums.tsv"
      zcat "${f}" | awk 'NR%2==0' | paste - - | sort -k1,1 | sha512sum | awk '{print $1 "\t" "'"${f}"'"}'
    done >> "!{meta.id}.Subsampled_FastQ.SHA512-checksums.tsv"

    msg "INFO: calculated checksums for !{meta.id}_R1.subsampled.fastq.gz !{meta.id}_R2.subsampled.fastq.gz"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
        seqkit: $(seqkit 2>&1 | grep "^Version: " | sed 's/^Version: //1')
    END_VERSIONS
    '''
}
