process MLST {

    publishDir "${params.outpath}/qa",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tab"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/mlst@sha256:27f290753760c44204d6e04b6ead7935d03b48d5f0a5ccce068def9ce33babe6"

    input:
        path base_fna

    output:
        path "Summary.MLST.tab"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''

        source bash_functions.sh

        # MLST for each assembly
        msg "INFO: Running MLST with !{task.cpus} threads"

        if [ -s !{base_fna} ]; then
            mlst --threads !{task.cpus} "!{base_fna}" \
            >> Summary.MLST.tab
        fi

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            mlst: $(mlst --version | awk 'NF>1{print $NF}')
        END_VERSIONS

        '''
}