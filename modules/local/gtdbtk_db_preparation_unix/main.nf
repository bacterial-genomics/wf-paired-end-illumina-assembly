process GTDBTK_DB_PREPARATION_UNIX {

    publishDir   "${params.process_log_dir}",
        mode:    "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs:  { filename -> "${task.process}${filename}"}

    label "process_medium"
    tag { "${database.getSimpleName()}" }

    container "ubuntu:jammy"

    input:
    path database

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                         , emit: versions
    tuple val("${database.getSimpleName()}"), path("database/*"), emit: db

    shell:
    '''
    mkdir database
    tar -xzf !{database} -C database --strip 1

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
