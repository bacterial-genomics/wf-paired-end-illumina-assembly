process SUBSAMPLE_READS_TO_DEPTH_SEQTK {

    tag { "${meta.id}" }
    container "staphb/seqtk@sha256:82797114adb664ba939b4f6dfcb822483a4af827def4288e5207be559055f2cc"

    input:
    tuple val(meta), path(reads), path(depth), path(fraction_of_reads)

    output:
    tuple val(meta), path("*.{fastq,fq}.gz", includeInputs: true), emit: reads
    path(".command.{out,err}")
    path("versions.yml")                                         , emit: versions

    shell:
    '''
    source bash_functions.sh

    fraction_of_reads_to_use=$(cat !{fraction_of_reads})
    initial_depth=$(cat !{depth})

    depth="!{params.depth}"

    if ! [[ ${fraction_of_reads_to_use} =~ ^[0-9.]+$ ]]; then
      msg "ERROR: Unable to calculate a fraction to subsample; ${fraction_of_reads_to_use} not a floating point value" >&2
      exit 1
    fi
    if [ ${depth%.*} -gt 0 ] && [ ${initial_depth%.*} -gt ${depth%.*} ]; then
      seqtk sample !{reads[0]} ${fraction_of_reads_to_use} > "!{meta.id}_R1.subsampled.fastq"
      seqtk sample !{reads[1]} ${fraction_of_reads_to_use} > "!{meta.id}_R2.subsampled.fastq"

      rm -f !{reads[0]} !{reads[1]}

      gzip -9f "!{meta.id}_R1.subsampled.fastq" \
        "!{meta.id}_R2.subsampled.fastq"

    else
      msg "INFO: Subsampling not requested or required"
      ###   NOTE:
      ###       this gonna be tricky?!
      ###       pass onto phix-remove-bbduk either subsampled reads() or initial input reads() tuple
    fi

    ### number of contigs and repeats elements with Lander-Waterman statistics
    ###  https://pubmed.ncbi.nlm.nih.gov/7497129/   ???

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        seqtk: $(seqtk 2>&1 | grep "^Version: " | sed 's/^Version: //1')
    END_VERSIONS
    '''
}
