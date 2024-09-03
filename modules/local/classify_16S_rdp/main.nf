process CLASSIFY_16S_RDP {

    label "process_medium"  // single CPU but needs the RAM boost
    tag { "${meta.id}" }
    container "tpaisie/rdp@sha256:ee388dff2e17c567946b7f2bf326765586d30f4ea0a203800616c44f599d53cc"

    input:
    tuple val(meta), path(barnapp_extracted_rna)

    output:
    tuple val(meta), path("${meta.id}.RDP_Classification_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}.RDP.tsv")                    , emit: rdp_tsv
    path(".command.{out,err}")
    path("versions.yml")                                           , emit: versions

    shell:
    // WARN: RDP does not report version information. This variable must be updated when container is updated.
    VERSION='2.14'
    '''
    source bash_functions.sh

    msg "INFO: Performing RDP 16S Classification of !{meta.id} ..."

    classifier \
      classify \
      --format "!{params.rdp_output_format}" \
      --gene "!{params.rdp_phylomarker}" \
      --outputFile "!{meta.id}.RDP.tsv" \
      "!{barnapp_extracted_rna}"

    msg "INFO: Completed RDP 16S Classification of !{meta.id}"

    if [[ "!{params.rdp_output_format}" == "fixrank" ]]; then
      # Drop unnecessary columns
      awk -F $'\t' 'BEGIN{OFS = FS} {print $1,$3,$5,$6,$8,$9,$11,$12,$14,$15,$17,$18,$20}' \
        "!{meta.id}.RDP.tsv" \
        > "!{meta.id}.RDP_tmp.tsv"

      mv -f "!{meta.id}.RDP_tmp.tsv" "!{meta.id}.RDP.tsv"

      # Add header
      sed -i \
        '1i Sample_name\tDomain\tDomain_result\tPhylum\tPhylum_result\tClass\tClass_result\tOrder\tOrder_result\tFamily\tFamily_result\tGenus\tGenus_result' \
        "!{meta.id}.RDP.tsv"
    else
      # Other `--format <arg>` options have varying numbers and names for header, so avoid adding any for now
      msg "WARN: RDP Classifier with `--format !{params.rdp_output_format}` unknown header column names might prevent downstream XLSX summary conversion"
    fi

    if verify_minimum_file_size "!{meta.id}.RDP.tsv" '16S Classification Output File' "!{params.min_filesize_rdp_output}"; then
      echo -e "!{meta.id}\t16S RDP Output File\tPASS" >> !{meta.id}.RDP_Classification_File.tsv
    else
      echo -e "!{meta.id}\t16S RDP Output File\tFAIL" >> !{meta.id}.RDP_Classification_File.tsv
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        rdp: $(echo !{VERSION})
    END_VERSIONS
    '''
}
