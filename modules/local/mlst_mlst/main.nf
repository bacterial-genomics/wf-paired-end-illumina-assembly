process MLST_MLST {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    tag { "${base}" }
    
    container "gregorysprenger/mlst@sha256:69c8c8027474b8f361ef4a579df171702f3ed52f45e3fb388a41ccbf4542706f"

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