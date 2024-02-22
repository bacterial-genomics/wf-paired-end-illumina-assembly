process CAT_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${meta.id}" }
    container "denolehov/curl:latest"

    input:
    tuple val(meta), path(database)

    output:
    path("database/{tax,db}") , emit: db
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    mkdir -p database
    # tar -xzf !{database} -C database/

    # Download webpage and parse HTML
    #  sort and take newest files (.tar.gz and .tar.gz.md5)
    curl -L https://tbb.bio.uu.nl/tina/CAT_prepare/ \
    | grep -o "<a href=['"'"'"][^'"'"'"]*" \
    | grep "tar.gz" \
    | sed 's|<a href="||g' \
    | sort -nr \
    | head -n 2 \
    | xargs \
        -0 \
        -I '{}' \
        sh -c 'wget -c https://tbb.bio.uu.nl/tina/CAT_prepare/{}'

    # Verify filename in md5 file matches downloaded database
    if [[ ! -f $(find . -name $(awk '{print $2}' *.tar.gz.md5 )) ]]; then
    echo "ERROR: Database file not found!"
    exit 1
    else
    [[ $(md5sum *.tar.gz | awk '{print $1}') = $(awk '{print $1}' *.tar.gz.md5) ]] \
        && echo "md5sum of files match!"; tar -xvzf CAT*.tar.gz \
        || echo "md5sum of files DO NOT MATCH!"; exit 1
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
