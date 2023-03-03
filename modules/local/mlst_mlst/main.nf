process MLST_MLST {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    tag { "${base}" }
    
    container "snads/mlst@sha256:27f290753760c44204d6e04b6ead7935d03b48d5f0a5ccce068def9ce33babe6"

    input:
        tuple val(base), path(paired_bam), path(single_bam), path(base_fna)

    output:
        path "${base}.Summary.MLST.tab", emit: summary_mlst
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        source bash_functions.sh

        # MLST for each assembly
        msg "INFO: Running MLST with !{task.cpus} threads"

        if [ -s !{base_fna} ]; then
          mlst \
           --threads !{task.cpus} \
           "!{base_fna}" \
           >> !{base}.Summary.MLST.tab
        fi

        # Get process version
        cat <<-END_VERSIONS > versions.yml
       "!{task.process} (!{base})":
            mlst: $(mlst --version | awk 'NF>1{print $NF}')
        END_VERSIONS
        '''
}