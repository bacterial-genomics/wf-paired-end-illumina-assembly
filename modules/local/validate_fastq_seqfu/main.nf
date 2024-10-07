process VALIDATE_FASTQ_SEQFU {

    tag { "${meta.id}" }
    container "staphb/seqfu@sha256:20831d2727d0f613f753eb301e19b345f5c9ea82c23762cb78a0c273539a3647"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.Raw_Initial_FastQ_Format_Validation_File.tsv"), emit: qc_filecheck
    tuple val(meta), path(reads)                                                    , emit: input
    path(".command.{out,err}")
    path("versions.yml")                                                            , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Validating !{meta.id} FastQ input with SeqFu..."

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.Raw_Initial_FastQ_Format_Validation_File.tsv"

    msg "INFO: Checking for FastQ valid format in R1: !{reads[0]} and R2: !{reads[1]}"

    # https://telatin.github.io/seqfu2/tools/check.html#integrity-check
    # A single FASTQ file is considered valid if:
    #     1 - each record has the same sequence and quality length
    #     2 - only A,C,G,T,N characters are present in the sequence
    #
    # A paired-end set of FASTQ files is considered valid if:
    #       - each file is individually valid
    #     3 - the two files have the same number of sequences
    #     4 - the first and last sequence of both files has the same name (the last three characters are ignored if the remaining - sequence name is greater than 4 characters)
    #     5 - the first and last sequence of the two files are not identical (R1 != R2)
    # Deep check
    #     If you are parsing NGS files, i.e. FASTQ files, with four lines per record and you expect them to be accepted by any program, use --deep.
    seqfu check \
      --deep \
      --verbose \
      !{reads[0]} !{reads[1]}

    # Retain the exit code status by exiting the exit value after error message
    retVal=$?
    if [ $retVal -ne 0 ]; then
      msg "ERROR: FastQ format validation tests with SeqFu failed for: !{meta.id} with exit status code: ${retVal}" >&2
      echo -e "!{meta.id}\tRaw Initial FastQ (R1 and R2) Valid Format\tFAIL" >> "!{meta.id}.Raw_Initial_FastQ_Format_Validation_File.tsv"
      exit $retVal
    fi

    msg "INFO: SeqFu check on !{reads[0]} !{reads[1]} completed without errors, suggesting the pair is a valid read set."
    echo -e "!{meta.id}\tRaw Initial FastQ (R1 and R2) Valid Format\tPASS" >> "!{meta.id}.Raw_Initial_FastQ_Format_Validation_File.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        seqfu: $(seqfu --version)
    END_VERSIONS
    '''
}
