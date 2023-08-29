process GTDBTK_DB_PREPARATION_UNIX {

    publishDir   "${params.process_log_dir}",
        mode:    "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs:  { filename -> "${task.process}${filename}"}

    label "process_medium"
    tag { "${db_name}" }

    container "ubuntu:jammy"

    input:
    path(database)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                        , emit: versions
    tuple val("${db_name}"), path("database/*"), emit: gtdb_db

    shell:
    db_name = database.toString().split('\\.')[0]
    '''
    if [ -d !{database} ]; then
      ln -sf !{database}/ database
    else
      mkdir database
      tar -xzf !{database} -C database --strip 1
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
