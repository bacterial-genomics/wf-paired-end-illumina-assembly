process QA_ASSEMBLY_QUAST {

    label "process_low"
    tag { "${meta.id}-${meta.assembler}" }
    container "staphb/quast@sha256:83ea0fd6c28ca01508fd7a93c0942b19089e6de25c9b8a496a34e138d240e0e8"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.QuastSummary.tsv"), emit: qa_summaries
    path("${meta.id}-${meta.assembler}.QuastSummary.tsv")                 , emit: summary_assemblies
    path(".command.{out,err}")
    path("versions.yml")                                                  , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Run Quast
    msg "INFO: Evaluating !{meta.id} assembly using QUAST ..."

    quast.py \
      --silent \
      --no-html \
      --no-plots \
      --min-contig 100 \
      --output-dir quast \
      --threads !{task.cpus} \
      --contig-thresholds 500,1000 \
      "!{assembly}"

    msg "INFO: Completed QUAST evaluation of !{meta.id} assembly"

    mv -f quast/transposed_report.tsv "!{meta.id}-!{meta.assembler}.QuastSummary.tsv"

    # Quast modifies basename. Need to check and modify if needed.
    assemblies_name=$(awk '{print $1}' "!{meta.id}-!{meta.assembler}.QuastSummary.tsv" | awk 'NR!=1 {print}')
    if [ ${assemblies_name} != !{meta.id} ]; then
      sed -i "s|${assemblies_name}|!{meta.id}|1" "!{meta.id}-!{meta.assembler}.QuastSummary.tsv"
    fi

    # Keep same first column header column name as all others -- "Sample_name"
    sed -i '1s/^Assembly/Sample_name/1' "!{meta.id}-!{meta.assembler}.QuastSummary.tsv"

    # Replace space characters in header line with underscores
    sed -i '1s/ /_/g' "!{meta.id}-!{meta.assembler}.QuastSummary.tsv"

    msg "INFO: Completed QUAST output renaming for !{meta.id}"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        quast: $(quast.py --version | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
