process SUBSAMPLE_READS_TO_DEPTH_SEQTK {

    label "process_medium"
    tag { "${meta.id}" }
    container "staphb/seqtk@sha256:e3105ea1c7375e6bfe0603f6e031b022068b3d4d529f295c5fa24e0a6709dd2c"

    input:
    tuple val(meta), path(reads), path(depth), path(fraction_of_reads)

    output:
    tuple val(meta), path("*.{fastq,fq}.gz", includeInputs: true), emit: reads
    path("${meta.id}.Subsampled_FastQ.SHA256-checksums.tsv")     , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                         , emit: versions

    shell:
    seqtk_seed = (params.seqtk_seed >= 1)? params.seqtk_seed : 947266746

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
      msg "INFO: Subsampling !{meta.id} R1 with seqtk using seed:!{params.seqtk_seed} ..."

      seqtk sample \
        -s "!{params.seqtk_seed}" \
        !{reads[0]} \
        ${fraction_of_reads_to_use} \
        > "!{meta.id}_R1.subsampled.fastq"

      msg "INFO: Subsampling !{meta.id} R2 with seqtk using seed:!{params.seqtk_seed} ..."

      seqtk sample \
        -s "!{params.seqtk_seed}" \
        !{reads[1]} \
        ${fraction_of_reads_to_use} \
        > "!{meta.id}_R2.subsampled.fastq"

      # Discard symlink infiles to avoid them being passed as outfiles when
      #   subsampling occurred.
      rm -f !{reads[0]} !{reads[1]}

      gzip -9f \
        "!{meta.id}_R1.subsampled.fastq" \
        "!{meta.id}_R2.subsampled.fastq"

      msg "INFO: Subsampled !{meta.id} R1 and R2 with seqtk"

    else
      # The input FastQ files that were never subsampled will get passed on
      #   as outputs here with the 'includeInputs: true'
      msg "INFO: Subsampling not requested or required"
      touch "!{meta.id}.Subsampled_FastQ.SHA256-checksums.tsv" versions.yml
      exit 0
    fi

    ### number of contigs and repeats elements with Lander-Waterman statistics
    ###  https://pubmed.ncbi.nlm.nih.gov/7497129/   ???

    ### Calculate SHA-256 Checksums of each FastQ file ###
    SUMMARY_HEADER=(
      "Sample_name"
      "Checksum_(SHA-256)"
      "File"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    echo "${SUMMARY_HEADER}" > "!{meta.id}.Subsampled_FastQ.SHA256-checksums.tsv"

    # Calculate checksums
    for f in "!{meta.id}_R1.subsampled.fastq.gz" "!{meta.id}_R2.subsampled.fastq.gz"; do
      echo -ne "!{meta.id}\t" >> "!{meta.id}.Subsampled_FastQ.SHA256-checksums.tsv"
      zcat "${f}" | awk 'NR%2==0' | paste - - | sort -k1,1 | sha256sum | awk '{print $1 "\t" "'"${f}"'"}'
    done >> "!{meta.id}.Subsampled_FastQ.SHA256-checksums.tsv"

    msg "INFO: calculated checksums for !{meta.id}_R1.subsampled.fastq.gz !{meta.id}_R2.subsampled.fastq.gz"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sha256sum: $(sha256sum --version | grep ^sha256sum | sed 's/sha256sum //1')
        seqtk: $(seqtk 2>&1 | grep "^Version: " | sed 's/^Version: //1')
    END_VERSIONS
    '''
}
