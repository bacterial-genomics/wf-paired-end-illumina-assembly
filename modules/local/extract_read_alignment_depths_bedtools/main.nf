process EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS {

    tag { "${meta.id}-${meta.assembler}" }
    container "snads/bedtools@sha256:9b80fb5c5ef1b6f4a4a211d8739fa3fe107da34d1fb6609d6b70ddc7afdce12c"

    input:
    tuple val(meta), path(bam_files)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.CleanedReads-AlnStats.tsv"), emit: summary
    path(".command.{out,err}")
    path("versions.yml")                                                           , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Calculate and report coverage of paired-reads and singleton reads separately
    msg "INFO: Extracting read alignment depths using bedtools"

    single_cov='0 bp TooFewToMap Singleton Reads (0.0x)\t'
    if [ -s !{bam_files[1]} ]; then
      single_cov=$(bedtools genomecov -d -split -ibam !{bam_files[1]} | \
                    awk '
                    BEGIN {
                        # Initialize the sum and count
                        sum = 0
                        count = 0
                    }
                    {
                        # Accumulate the sum of the third column (depth at each position)
                        sum += $3

                        # Store each value in an array for standard deviation calculation
                        values[NR] = $3

                        # Increment the count of records
                        count += 1
                    }
                    END {
                        # Calculate the mean
                        mean = sum / count

                        # Calculate the sum of squared differences from the mean for standard deviation
                        sum_of_squares = 0
                        for (i = 1; i <= count; i++) {
                            sum_of_squares += (values[i] - mean) * (values[i] - mean)
                        }
                        stddev = sqrt(sum_of_squares / count)

                        # Print total sum, mean, and standard deviation on one line
                        printf "%d bp Singleton Reads Mapped (%.2fx ± %.2fx)\n", sum, mean, stddev
                    }'
                   )
    fi

    cov_info=$(bedtools genomecov -d -split -ibam "!{bam_files[0]}" |\
                awk -v OFS='\t' -v SEcov="${single_cov}" '
                BEGIN {
                    sum = 0
                    count = 0
                }
                {
                    sum += $3  
                    values[NR] = $3    
                    count += 1
                }
                END {
                    mean = sum / count

                    sum_of_squares = 0
                    for (i = 1; i <= count; i++) {
                        sum_of_squares += (values[i] - mean) * (values[i] - mean)
                    }
                    stddev = sqrt(sum_of_squares / count)

                    printf("%d bp Paired Reads Mapped (%.2fx ± %.2fx)\t%s\t%d bp Genome Assembly Length\n", sum, mean, stddev, SEcov, count)
                }'
              )

    echo -e "!{meta.id}\t${cov_info}" \
      > "!{meta.id}-!{meta.assembler}.CleanedReads-AlnStats.tsv"

    sed -i \
      '1i Sample_name\tCoverage_of_paired_reads_(mean_±_stdev)\tCoverage_of_singleton_reads_(mean_±_stdev)\tGenome_assembly_length_(bp)' \
      "!{meta.id}-!{meta.assembler}.CleanedReads-AlnStats.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bedtools: $(bedtools --version | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
