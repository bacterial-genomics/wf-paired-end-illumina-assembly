process BUSCO_DB_PREPARATION_UNIX {

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
    path "versions.yml", emit: versions
    path "database/"   , emit: db

    shell:
    '''
    mkdir db_tmp
    tar -xzf !{database} -C db_tmp/

    # Place lineage datasets in correct format for BUSCO v5
    mkdir -p database/lineages
    mv \
      `find db_tmp/ -maxdepth 2 -type d -name "*_odb10"` \
      database/lineages

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
