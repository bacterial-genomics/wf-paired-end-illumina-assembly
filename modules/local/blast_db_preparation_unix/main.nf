process BLAST_DB_PREPARATION_UNIX {

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
    mkdir -p database
    tar -xzf !{database} -C database/

    # Make sure database contains 16S_ribosomal_RNA files
    if [[ $(find database/ -type f -name "16S_ribosomal_RNA*" | wc -l) -lt 1 ]]; then
        msg "ERROR: Missing 16S ribosomal RNA database files from NCBI BLAST!" >&2
        exit 1
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
