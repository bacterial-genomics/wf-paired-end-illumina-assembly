process EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS {

    tag { "${meta.id}-${meta.assembler}" }
    container "staphb/bedtools@sha256:52d4a9359d3adaa6ac8f8ebbdc5596bae791a738974e9a85f72892486a43336e"

    input:
    tuple val(meta), path(bam_files)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Clean_Reads-AlnStats.tsv"), emit: summary
    path(".command.{out,err}")
    path("versions.yml")                                                           , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Calculate and report coverage of paired-reads and singleton reads separately, in addition to combined

    SINGLE_STATS=(0 0 0 0)
    echo -n '' > single_coverage.tsv
    if [ -s !{bam_files[1]} ]; then
      msg "INFO: Extracting singleton read alignment depths using bedtools from !{bam_files[1]}"
      bedtools genomecov -d -split -ibam !{bam_files[1]} > single_coverage.tsv

      total_sites=$(wc -l < single_coverage.tsv)
      sum_single_coverage=$(awk '{sum+=$3} END{print sum}' single_coverage.tsv)
      mean_single_coverage=$(awk -v sum="$sum_single_coverage" -v total="$total_sites" 'BEGIN{printf "%.1f", sum/total}')
      stdev_single_coverage=$(awk -v mean="$mean_single_coverage" '{sum+=($3-mean)^2} END{printf "%.1f", sqrt(sum/NR)}' single_coverage.tsv)

      SINGLE_STATS=($total_sites $sum_single_coverage $mean_single_coverage $stdev_single_coverage)
    fi

    msg "INFO: Extracting paired-end read alignment depths using bedtools from !{bam_files[0]}"
    bedtools genomecov -d -split -ibam !{bam_files[0]} > paired_coverage.tsv

    total_sites=$(wc -l < paired_coverage.tsv)
    sum_paired_coverage=$(awk '{sum+=$3} END{print sum}' paired_coverage.tsv)
    mean_paired_coverage=$(awk -v sum="$sum_paired_coverage" -v total="$total_sites" 'BEGIN {printf "%.1f", sum/total}')
    stdev_paired_coverage=$(awk -v mean="$mean_paired_coverage" '{sum+=($3-mean)^2} END{printf "%.1f", sqrt(sum/NR)}' paired_coverage.tsv)

    PAIRED_STATS=($total_sites $sum_paired_coverage $mean_paired_coverage $stdev_paired_coverage)

    # Combine coverage values to report total depth of coverage statistics
    sum_total_coverage=$(awk '{sum+=$3} END{print sum}' {paired,single}_coverage.tsv)
    mean_total_coverage=$(awk -v sum="$sum_paired_coverage" -v total="$total_sites" 'BEGIN {printf "%.1f", sum/total}')
    stdev_total_coverage=$(awk -v mean="$mean_paired_coverage" '{sum+=($3-mean)^2} END{printf "%.1f", sqrt(sum/NR)}' {paired,single}_coverage.tsv)

    TOTAL_STATS=($total_sites $sum_total_coverage $mean_total_coverage $stdev_total_coverage)

    COVERAGE_DATA=(
      "!{meta.id}"
      "${PAIRED_STATS[1]}"
      "${PAIRED_STATS[2]}"
      "${PAIRED_STATS[3]}"
      "${SINGLE_STATS[1]}"
      "${SINGLE_STATS[2]}"
      "${SINGLE_STATS[3]}"
      "${TOTAL_STATS[1]}"
      "${TOTAL_STATS[2]}"
      "${TOTAL_STATS[3]}"
      "${TOTAL_STATS[0]}"
    )
    COVERAGE_DATA=$(printf "%s\t" "${COVERAGE_DATA[@]}" | sed 's/\t$//')

    SUMMARY_HEADER=(
      "Sample_name"
      "Total_mapped_paired_reads_(bp)"
      "Mean_coverage_of_paired_reads_(x)"
      "Stdev_coverage_of_paired_reads_(x)"
      "Total_mapped_singleton_reads_(bp)"
      "Mean_coverage_of_singleton_reads_(x)"
      "Stdev_coverage_of_singleton_reads_(x)"
      "Total_mapped_paired_and_singleton_reads_(bp)"
      "Mean_coverage_of_paired_and_singleton_reads_(x)"
      "Stdev_coverage_of_paired_and_singleton_reads_(x)"
      "Genome_assembly_length_(bp)"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    echo -e "${SUMMARY_HEADER}\n${COVERAGE_DATA}" > "!{meta.id}-!{meta.assembler}.Clean_Reads-AlnStats.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bedtools: $(bedtools --version | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
