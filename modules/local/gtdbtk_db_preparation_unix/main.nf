process GTDBTK_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${database.getSimpleName()}" }
    container "ubuntu:jammy"

    input:
    path(database)

    output:
    path(".command.{out,err}")
    path("versions.yml")                                        , emit: versions
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
