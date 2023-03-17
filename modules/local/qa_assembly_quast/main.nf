process QA_ASSEMBLY_QUAST {

    // errorStrategy 'terminate'

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}
    
    label "process_low"
    tag { "${base}" }
    
    container "snads/quast@sha256:c8147a279feafbc88bafeeda3817ff32d43db87d31dd0978df1cd2f8022d324c"

    input:
        tuple val(base), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck), path(base_fna)

    output:
        tuple val(base), path("${base}.Summary.Assemblies.tab"), path("${base}.Summary.Illumina.CleanedReads-Bases.tab"), emit: qa_summaries
        path "${base}.Summary.Assemblies.tab", emit: summary_assemblies
        path "${base}.Summary.Illumina.CleanedReads-Bases.tab", emit: summary_reads
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        source bash_functions.sh

        # Exit if previous process fails qc filecheck
        for filecheck in !{qc_nonoverlap_filecheck}; do
          if [[ $(grep "FAIL" ${filecheck}) ]]; then
            error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
            msg "FAILURE: ${error_message} Check FAILED" >&2
            exit 1
          else
            rm ${filecheck}
          fi
        done

        # Run Quast
        msg "INFO: Running QUAST with !{task.cpus} threads"

        quast.py \
         --output-dir quast \
         --min-contig 100 \
         --threads !{task.cpus} \
         --no-html \
         --gene-finding \
         --gene-thresholds 300 \
         --contig-thresholds 500,1000 \
         --ambiguity-usage one \
         --strict-NA \
         --silent \
         "!{base_fna}" >&2

        mv -f quast/transposed_report.tsv !{base}.Summary.Assemblies.tab

        # Quast modifies basename. Need to check and modify if needed.
        assemblies_name=$(awk '{print $1}' !{base}.Summary.Assemblies.tab | awk 'NR!=1 {print}')
        if [ ${assemblies_name} != !{base} ]; then
          sed -i "s|${assemblies_name}|!{base}|g" !{base}.Summary.Assemblies.tab
        fi

        # TO-DO: move this unix-only component to separate QA_READS_BASEPAIR_COUNT_UNIX
        # Count nucleotides per read set
        echo -n '' > Summary.Illumina.CleanedReads-Bases.tab
        for (( i=0; i<3; i+=3 )); do
          R1=$(basename "!{paired_R1_gz}" _R1.paired.fq.gz)
          R2=$(basename "!{paired_R2_gz}" _R2.paired.fq.gz)
          single=$(basename "!{single_gz}" .single.fq.gz)

          # Verify each set of reads groups properly
          nr_uniq_str=$(echo -e "${R1}\n${R2}\n${single}" | sort -u | wc -l)
          if [ "${nr_uniq_str}" -ne 1 ]; then
            msg "ERROR: improperly grouped ${R1} ${R2} ${single}" >&2
            exit 1
          fi
          echo -ne "${R1}\t" >> !{base}.Summary.Illumina.CleanedReads-Bases.tab
          zcat "!{paired_R1_gz}" "!{paired_R2_gz}" "!{single_gz}" | \
           awk 'BEGIN{SUM=0} {if(NR%4==2){SUM+=length($0)}} END{print SUM}' \
           >> !{base}.Summary.Illumina.CleanedReads-Bases.tab
        done

        # Get process version
        echo -e '"!{task.process} (!{base})":' > versions.yml
        echo -e "    quast: $(quast.py --version | awk 'NF>1{print $NF}')" >> versions.yml
        '''
}