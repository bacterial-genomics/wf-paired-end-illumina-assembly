process CHECKM2_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${meta.id}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(database)

    output:
    path("database/*.dmnd")   , emit: db
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    mkdir -p tmp_db database
    tar -xzf !{database} --strip 1 -C tmp_db

    mv `find tmp_db -name "*.dmnd"` database/

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
