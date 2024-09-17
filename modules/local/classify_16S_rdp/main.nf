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
      --outputFile "!{meta.id}.RDP-raw.tsv" \
      "!{barnapp_extracted_rna}"

    msg "INFO: Completed RDP 16S Classification of !{meta.id}"

  if [[ "!{params.rdp_output_format}" == "fixrank" ]]; then
      # Split up the first column "<Sample_name>_<int>"; add header;
      #   discard some columns that are now stored as headers
      #   (e.g., "domain", "phylum", "class", "order", "family", "genus")
      awk 'BEGIN {
        FS=OFS="\t";
        print "Sample_name\tUnique_16S_rRNA_extraction_count\tDomain_result\tPhylum_result\tClass_result\tOrder_result\tFamily_result\tGenus_result"
      }
      {
        split($1, a, "_");
        print a[1], a[2], $3, $6, $9, $12, $15, $18
      }' "!{meta.id}.RDP-raw.tsv" \
      > "!{meta.id}.RDP.tsv"

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
