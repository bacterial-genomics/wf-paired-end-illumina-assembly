process CLASSIFY_16S_RDP {

    tag { "${meta.id}" }
    container "tpaisie/rdp@sha256:ee388dff2e17c567946b7f2bf326765586d30f4ea0a203800616c44f599d53cc"

    input:
    tuple val(meta), path(barnapp_extracted_rna)

    output:
    path(".command.{out,err}")
    path "versions.yml",                             emit: versions
    path "${meta.id}.rdp.tsv",                       emit: qc_rdp_filecheck
    tuple val(meta), path("${meta.id}.rdp.tsv"),     emit: rdp_tsv

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Performing RDP 16S Classification"

###         ##############################################################
###         ##############################################################
###  NOTE: pin this to just 1 CPU, low RAM, < 1 min runtime
###         ##############################################################
###         ##############################################################

    classifier \
      classify \
      --format fixrank \
      --outputFile "!{meta.id}.rdp.tsv" \
      "!{barnapp_extracted_rna}"

    if verify_minimum_file_size "!{meta.id}.rdp.tsv" '16S Classification Output File' "!{params.min_filesize_rdp_output}"; then
      echo -e "!{meta.id}\t16S RDP Output File\tPASS" > !{meta.id}.rdp.tsv
    else
      echo -e "!{meta.id}\t16S RDP Output File\tFAIL" > !{meta.id}.rdp.tsv
    fi

    # Get process version information
    #cat <<-END_VERSIONS > versions.yml
    #"!{task.process}":
        #rdp: $(rdp_classifier version')
    #END_VERSIONS

    '''
}
