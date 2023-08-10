process QA_ASSEMBLY_QUAST {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${meta.id}.${task.process}${filename}"}

    label "process_low"
    tag { "${meta.id}" }

    container "snads/quast@sha256:c8147a279feafbc88bafeeda3817ff32d43db87d31dd0978df1cd2f8022d324c"

    input:
      tuple val(meta), path(R1), path(R2), path(single), path(qc_nonoverlap_filecheck), path(assembly)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                                                                                   , emit: versions
    path "${meta.id}.Summary.Assemblies.tab"                                                                              , emit: summary_assemblies
    path "${meta.id}.Summary.Illumina.CleanedReads-Bases.tab"                                                             , emit: summary_reads
    tuple val(meta), path("${meta.id}.Summary.Assemblies.tab"), path("${meta.id}.Summary.Illumina.CleanedReads-Bases.tab"), emit: qa_summaries

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_nonoverlap_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Run Quast
    msg "INFO: Running QUAST with !{task.cpus} threads"

    quast.py \
      --silent \
      --no-html \
      --strict-NA \
      --gene-finding \
      --min-contig 100 \
      --output-dir quast \
      --gene-thresholds 300 \
      --ambiguity-usage one \
      --threads !{task.cpus} \
      --contig-thresholds 500,1000 \
      "!{assembly}" >&2

    mv -f quast/transposed_report.tsv !{meta.id}.Summary.Assemblies.tab

    # Quast modifies basename. Need to check and modify if needed.
    assemblies_name=$(awk '{print $1}' !{meta.id}.Summary.Assemblies.tab | awk 'NR!=1 {print}')
    if [ ${assemblies_name} != !{meta.id} ]; then
      sed -i "s|${assemblies_name}|!{meta.id}|g" !{meta.id}.Summary.Assemblies.tab
    fi

    # TO-DO: move this unix-only component to separate QA_READS_BASEPAIR_COUNT_UNIX
    # Count nucleotides per read set
    echo -n '' > Summary.Illumina.CleanedReads-Bases.tab
    for (( i=0; i<3; i+=3 )); do
      R1=$(basename "!{R1}" _R1.paired.fq.gz)
      R2=$(basename "!{R2}" _R2.paired.fq.gz)
      single=$(basename "!{single}" .single.fq.gz)

      # Verify each set of reads groups properly
      nr_uniq_str=$(echo -e "${R1}\n${R2}\n${single}" | sort -u | wc -l)
      if [ "${nr_uniq_str}" -ne 1 ]; then
        msg "ERROR: improperly grouped ${R1} ${R2} ${single}" >&2
        exit 1
      fi
      echo -ne "${R1}\t" >> !{meta.id}.Summary.Illumina.CleanedReads-Bases.tab
      zcat "!{R1}" "!{R2}" "!{single}" | \
        awk 'BEGIN{SUM=0} {if(NR%4==2){SUM+=length($0)}} END{print SUM}' \
          >> !{meta.id}.Summary.Illumina.CleanedReads-Bases.tab
    done

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      quast: $(quast.py --version | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
