process EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS {

    tag { "${meta.id}-${meta.assembler}" }
    container "snads/bedtools@sha256:9b80fb5c5ef1b6f4a4a211d8739fa3fe107da34d1fb6609d6b70ddc7afdce12c"

    input:
    tuple val(meta), path(bam_files)

    output:
    path(".command.{out,err}")
    path("versions.yml")                                                           , emit: versions
    tuple val(meta), path("${meta.id}-${meta.assembler}.CleanedReads-AlnStats.tsv"), emit: summary_alignment_stats

    shell:
    '''
    source bash_functions.sh

    # Calculate and report coverage of paired-reads and singleton reads separately
    msg "INFO: Extracting read alignment depths using bedtools"

    single_cov='0 bp TooFewToMap Singleton Reads (0.0x)\t'
    if [ -s !{bam_files[1]} ]; then
      single_cov=$(bedtools genomecov -d -split -ibam !{bam_files[1]} |\
        awk '{sum+=$3} END{print sum " bp Singleton Reads Mapped (" sum/NR "x)\t"}')
    fi

    cov_info=$(bedtools genomecov -d -split -ibam "!{bam_files[0]}" |\
      awk -v OFS='\t' -v SEcov="${single_cov}" 'BEGIN{sum=0} {sum+=$3} END{
      print sum " bp Paired Reads Mapped (" sum/NR "x)\t" SEcov NR " bp Genome"}')

    echo -e "!{meta.id}\t${cov_info}" \
      > "!{meta.id}-!{meta.assembler}.CleanedReads-AlnStats.tsv"

    sed -i \
      '1i Sample name\tCoverage of paired reads\tCoverage of singleton reads\tGenome size' \
      "!{meta.id}-!{meta.assembler}.CleanedReads-AlnStats.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bedtools: $(bedtools --version | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
