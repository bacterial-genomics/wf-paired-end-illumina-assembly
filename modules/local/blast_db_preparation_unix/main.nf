process BLAST_DB_PREPARATION_UNIX {

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
    path "versions.yml"                , emit: versions
    path "database/16S_ribosomal_RNA.*", emit: db

    shell:
    '''
    mkdir db_tmp
    tar -xzf !{database} -C db_tmp

    mkdir database
    mv `find db_tmp/ -name "16S_ribosomal_RNA*"` database/

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
