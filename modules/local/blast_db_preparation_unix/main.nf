process BLAST_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${db_name}" }
    container "ubuntu:jammy"

    input:
    path(database)

    output:
    tuple val(db_name), path("database/*"), emit: db
    path(".command.{out,err}")
    path("versions.yml")                  , emit: versions

    shell:
    db_name = "16S_ribosomal_RNA"
    '''
    mkdir -p database
    tar -xzf !{database} -C database/

    # Make sure database contains 16S_ribosomal_RNA files
    if [[ $(find database/ -type f -name "16S_ribosomal_RNA*" | wc -l) -lt 1 ]]; then
        msg "ERROR: Missing 16S ribosomal RNA database files from NCBI BLAST!"
        exit 1
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
