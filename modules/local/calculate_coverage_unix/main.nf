process CALCULATE_COVERAGE_UNIX {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    tag { "${prefix}" }

    container "ubuntu:jammy"

    input:
    tuple val(prefix), path(summary_assemblies), path(summary_reads), path(summary_stats)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                 , emit: versions
    path "${prefix}.Summary.Illumina.GenomeCoverage.tab", emit: genome_coverage

    shell:
    '''
    source bash_functions.sh

    # Report coverage
    echo -n '' > !{prefix}.Summary.Illumina.GenomeCoverage.tab
    i=0
    while IFS=$'\t' read -r -a ln; do
      if grep -q -e "skesa_" -e "unicyc_" -e ".uncorrected" <<< "${ln[0]}"; then
        continue
      fi

      basepairs=$(grep ${ln[0]} !{summary_stats} 2> /dev/null \
        | awk 'BEGIN{FS="\t"}; {print $2}' | awk '{print $1}' | sort -u)

      if [[ "${basepairs}" =~ ^[0-9]+$ ]]; then
        msg "INFO: read alignment data for ${ln[0]} used for coverage" >&2
      else
        basepairs=$(grep ${ln[0]} !{summary_reads} | cut -f 2)
        msg "INFO: read alignment data absent for ${ln[0]}, so cleaned bases" >&2
        msg "      given to the assembler were used to calculate coverage" >&2
      fi

      genomelen=${ln[7]}
      cov=$(echo | awk -v x=${basepairs} -v y=${genomelen} '{printf ("%0.1f", x/y)}')
      msg "INFO: cov = $cov"

      if [[ "${cov}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo -e "${ln[0]}\t${cov}x" >> !{prefix}.Summary.Illumina.GenomeCoverage.tab
        ((i=i+1))
      fi
    done < <(grep -v 'Total length' !{summary_assemblies})

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
