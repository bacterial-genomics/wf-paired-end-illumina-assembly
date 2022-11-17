process GENOME_COVERAGE {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    container "ubuntu:focal"
    
    input:
        path summary_stats
        path summary_assemblies
        path summary_bases
        val base

    output:
        path "*.Summary.Illumina.GenomeCoverage.tab", emit: genome_coverage
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        source bash_functions.sh

        # Report coverage
        echo -n '' > !{base}.Summary.Illumina.GenomeCoverage.tab
        i=0
        while IFS=$'\t' read -r -a ln; do
            if grep -q -e "skesa_" -e "unicyc_" -e ".uncorrected" <<< "${ln[0]}"; then
                continue
            fi

            basepairs=$(grep ${ln[0]} !{summary_stats} \
            2> /dev/null | awk 'BEGIN{FS="\t"}; {print $2}' | awk '{print $1}' | sort -u)

            if [[ "${basepairs}" =~ ^[0-9]+$ ]]; then
                msg "INFO: read alignment data for ${ln[0]} used for coverage" >&2
            else
                basepairs=$(grep ${ln[0]} !{summary_bases} | cut -f 2)
                msg "INFO: read alignment data absent for ${ln[0]}, so cleaned bases" >&2
                msg "      given to the assembler were used to calculate coverage" >&2
            fi

            genomelen=${ln[7]}
            cov=$(echo | awk -v x=${basepairs} -v y=${genomelen} '{printf ("%0.1f", x/y)}')
            msg "INFO: cov = $cov"
            
            if [[ "${cov}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                echo -e "${ln[0]}\t${cov}x" >> !{base}.Summary.Illumina.GenomeCoverage.tab
                ((i=i+1))
            fi
        done < <(grep -v 'Total length' !{summary_assemblies})

        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            ubuntu: $(awk -F ' ' '{print $1,$2,$3}' /etc/issue | tr -d '\\n')
        END_VERSIONS
        '''
}