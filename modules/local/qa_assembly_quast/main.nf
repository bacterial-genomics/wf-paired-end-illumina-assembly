process QA_ASSEMBLY_QUAST {

    label "process_low"
    tag { "${meta.id}-${meta.assembler}" }
    container "snads/quast@sha256:c8147a279feafbc88bafeeda3817ff32d43db87d31dd0978df1cd2f8022d324c"

    input:
    tuple val(meta), path(cleaned_fastq_files), path(assembly)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.QuastSummary.tsv"), path("${meta.id}-${meta.assembler}.CleanedReads-Bases.tsv"), emit: qa_summaries
    path("${meta.id}-${meta.assembler}.QuastSummary.tsv")                                                                              , emit: summary_assemblies
    path("${meta.id}-${meta.assembler}.CleanedReads-Bases.tsv")                                                                        , emit: summary_reads
    path(".command.{out,err}")
    path("versions.yml")                                                                                                               , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Run Quast
    msg "INFO: Evaluating assembly using QUAST"

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

    mv -f quast/transposed_report.tsv "!{meta.id}-!{meta.assembler}.QuastSummary.tsv"

    # Quast modifies basename. Need to check and modify if needed.
    assemblies_name=$(awk '{print $1}' "!{meta.id}-!{meta.assembler}.QuastSummary.tsv" | awk 'NR!=1 {print}')
    if [ ${assemblies_name} != !{meta.id} ]; then
      sed -i "s|${assemblies_name}|!{meta.id}|g" "!{meta.id}-!{meta.assembler}.QuastSummary.tsv"
    fi

    # TO-DO: move this unix-only component to separate QA_READS_BASEPAIR_COUNT_UNIX
    # Count nucleotides per read set
    echo -n '' > "!{meta.id}-!{meta.assembler}.CleanedReads-Bases.tsv"
    for (( i=0; i<3; i+=3 )); do
      R1=$(basename "!{cleaned_fastq_files[0]}" _R1.paired.fq.gz)
      R2=$(basename "!{cleaned_fastq_files[1]}" _R2.paired.fq.gz)
      single=$(basename "!{cleaned_fastq_files[2]}" _single.fq.gz)

      # Verify each set of reads groups properly
      nr_uniq_str=$(echo -e "${R1}\\n${R2}\\n${single}" | sort -u | wc -l)
      if [ "${nr_uniq_str}" -ne 1 ]; then
        msg "ERROR: improperly grouped ${R1} ${R2} ${single}" >&2
        exit 1
      fi
      echo -ne "${R1}\t" >> "!{meta.id}-!{meta.assembler}.CleanedReads-Bases.tsv"
      zcat "!{cleaned_fastq_files[0]}" "!{cleaned_fastq_files[1]}" "!{cleaned_fastq_files[2]}" | \
        awk 'BEGIN{SUM=0} {if(NR%4==2){SUM+=length($0)}} END{OFMT="%f"; print SUM}' \
          >> "!{meta.id}-!{meta.assembler}.CleanedReads-Bases.tsv"

      sed -i '1i Sample name\t# cleaned bases' "!{meta.id}-!{meta.assembler}.CleanedReads-Bases.tsv"
    done

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        quast: $(quast.py --version | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
