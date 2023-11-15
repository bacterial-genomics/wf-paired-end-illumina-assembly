process BLAST_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${db_name}" }
    container "ubuntu:jammy"

    input:
    path(database)

    output:
    path(".command.{out,err}")
    path "versions.yml"                                     , emit: versions
    tuple val(db_name), path("database/16S_ribosomal_RNA.*"), emit: db

    shell:
    db_name = "16S_ribosomal_RNA"
    '''
    mkdir -p db_tmp
    tar -xzf !{database} -C db_tmp/

    mkdir -p database
    mv `find db_tmp/ -name "16S_ribosomal_RNA*"` database/

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
