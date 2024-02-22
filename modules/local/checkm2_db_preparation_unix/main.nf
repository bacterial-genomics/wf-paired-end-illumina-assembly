process CHECKM2_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${meta.id}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(database)

    output:
    path("*.dmnd")            , emit: db
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    tar -xzf !{database} --strip 1

    mv -fv $(find . -name "*.dmnd") .

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
