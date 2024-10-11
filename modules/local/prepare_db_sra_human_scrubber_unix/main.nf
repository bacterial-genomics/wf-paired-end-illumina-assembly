process PREPARE_DB_SRA_HUMAN_SCRUBBER {

    label "process_low"
    tag { "${database.getSimpleName()}" }
    container "ubuntu:jammy"

    input:
    path(database)

    output:
    path("${database.getSimpleName()}*"), emit: db
    path(".command.{out,err}")
    path("versions.yml")                , emit: versions

    shell:
    '''
    gunzip !{database}

    # Make sure the decompression worked
    if [[ $(find . -type f -name !{database%.gz} | wc -l) -lt 1 ]]; then
        msg "ERROR: Missing decompressed SRA Human Scrubber db file" >&2
        exit 1
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
