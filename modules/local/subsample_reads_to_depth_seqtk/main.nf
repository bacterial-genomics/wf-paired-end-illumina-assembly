process SUBSAMPLE_READS_TO_DEPTH_SEQTK {

    tag { "${meta.id}" }
    container "staphb/seqtk@sha256:e3105ea1c7375e6bfe0603f6e031b022068b3d4d529f295c5fa24e0a6709dd2c"

    input:
    tuple val(meta), path(reads), path(depth), path(fraction_of_reads)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                   , emit: versions
    tuple val(meta), path("*.fastq*", includeInputs: true), emit: reads

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
      seqtk sample !{reads[0]} ${fraction_of_reads_to_use} > !{meta.id}.R1.subsampled.fastq
      seqtk sample !{reads[1]} ${fraction_of_reads_to_use} > !{meta.id}.R2.subsampled.fastq

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
