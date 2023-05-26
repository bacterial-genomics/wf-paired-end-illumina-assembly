process EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    tag { "${prefix}" }

    container "snads/bedtools@sha256:9b80fb5c5ef1b6f4a4a211d8739fa3fe107da34d1fb6609d6b70ddc7afdce12c"

    input:
    tuple val(prefix), path(paired_bam), path(single_bam), path(qc_assembly_filecheck)

    output:
    tuple val(prefix), path("${prefix}.Summary.Illumina.CleanedReads-AlnStats.tab"), emit: summary_stats
    path "${prefix}.Summary.Illumina.CleanedReads-AlnStats.tab", emit: summary_alnstats
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_assembly_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Calculate and report coverage of paired-reads and singleton reads separately
    msg "INFO: Running bedtools"

    single_cov='0 bp TooFewToMap Singleton Reads (0.0x)\t'
    if [ -s !{single_bam} ]; then
      single_cov=$(bedtools genomecov -d -split -ibam !{single_bam} |\
        awk '{sum+=$3} END{print sum " bp Singleton Reads Mapped (" sum/NR "x)\t"}')
    fi

    cov_nfo=$(bedtools genomecov -d -split -ibam !{prefix}.paired.bam |\
      awk -v SEcov="${single_cov}" 'BEGIN{sum=0} {sum+=$3} END{
      print sum " bp Paired Reads Mapped (" sum/NR "x)\t" SEcov NR " bp Genome"}')

    echo -e "!{prefix}\t${cov_nfo}" \
      >> !{prefix}.Summary.Illumina.CleanedReads-AlnStats.tab

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      bedtools: $(bedtools --version | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
