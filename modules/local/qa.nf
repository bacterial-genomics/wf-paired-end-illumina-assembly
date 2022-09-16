process QA {

    publishDir "${params.outpath}/qa",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tab"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}
    
    label "process_low"
    
    container "snads/quast@sha256:c8147a279feafbc88bafeeda3817ff32d43db87d31dd0978df1cd2f8022d324c"

    input:
        path base_fna
        path R1_paired_gz
        path R2_paired_gz
        path single_gz
        val base

    output:
        path "Summary.Assemblies.tab", emit: summary_assemblies
        path "Summary.Illumina.CleanedReads-Bases.tab", emit: summary_bases
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''

        source bash_functions.sh

        # Run Quast
        msg "INFO: Running QUAST with !{task.cpus} threads"

        quast.py --output-dir quast --min-contig 100 --threads !{task.cpus} \
        --no-html --gene-finding --gene-thresholds 300 --contig-thresholds 500,1000 \
        --ambiguity-usage one --strict-NA --silent "!{base_fna}" >&2

        mv -f quast/transposed_report.tsv Summary.Assemblies.tab

        # Quast modifies basename. Need to check and modify if needed.
        assemblies_name=$(awk '{print $1}' Summary.Assemblies.tab | awk 'NR!=1 {print}')
        if [ ${assemblies_name} != !{base} ]; then
            sed -i "s|${assemblies_name}|!{base}|g" Summary.Assemblies.tab
        fi

        # Count nucleotides per read set
        echo -n '' > Summary.Illumina.CleanedReads-Bases.tab
        for (( i=0; i<3; i+=3 )); do
            R1=$(basename "!{R1_paired_gz}" _R1.paired.fq.gz)
            R2=$(basename "!{R2_paired_gz}" _R2.paired.fq.gz)
            single=$(basename "!{single_gz}" _single.fq.gz)

            # Verify each set of reads groups properly
            nr_uniq_str=$(echo -e "${R1}\n${R2}\n${single}" | sort -u | wc -l)
            if [ "${nr_uniq_str}" -ne 1 ]; then
                msg "ERROR: improperly grouped ${R1} ${R2} ${single}" >&2
                exit 1
            fi
            echo -ne "${R1}\t" >> Summary.Illumina.CleanedReads-Bases.tab
            zcat "!{R1_paired_gz}" "!{R2_paired_gz}" "!{single_gz}" | \
            awk 'BEGIN{SUM=0} {if(NR%4==2){SUM+=length($0)}} END{print SUM}' \
            >> Summary.Illumina.CleanedReads-Bases.tab
        done

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            quast: $(quast.py --version | awk 'NF>1{print $NF}')
        END_VERSIONS

        '''
}