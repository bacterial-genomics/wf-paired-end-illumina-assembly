process GTDBTK_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${meta.id}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(database)

    output:
    tuple val("${meta.id}"), path("database/*"), emit: db
    path(".command.{out,err}")
    path("versions.yml")                       , emit: versions

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
