process DOWNLOAD_CAT_DB_UNIX {

    label "process_medium"
    tag { "${meta.id}" }
    container "quay.io/biocontainers/aria2@sha256:acb0c86334ea0b2ba9454cc1b4d08f30c5a6ec7852159fd28fe34698154798d6"

    input:
    tuple val(meta)

    output:
    tuple val(meta), path("*.tar.gz"), emit: db
    path(".command.{out,err}")
    path("versions.yml")             , emit: versions

    shell:
    '''
    # Download webpage and parse HTML
    wget -qSO - https://tbb.bio.uu.nl/tina/CAT_prepare/ \
      | grep -o "<a href=['"'"'"][^'"'"'"]*" \
      | sed 's|<a href="||g' \
      | grep "tar.gz" \
      | grep "nr" \
      | while read fname; do
          aria2c \
            --continue \
            --max-connection-per-server=16 \
            --split=16 \
            --check-integrity \
            https://tbb.bio.uu.nl/tina/CAT_prepare/${fname}
        done

    ## IGNORE md5 check due to files not matching on database website!
    # Verify filename in md5 file matches downloaded database
    # if [[ ! -f $(find . -name $(awk '{print $2}' *.tar.gz.md5 )) ]]; then
    #   echo "ERROR: Database file not found!"
    #   exit 1
    # else
    #   [[ $(md5sum *.tar.gz | awk '{print $1}') = $(awk '{print $1}' *.tar.gz.md5) ]] \
    #     && echo "md5sum of files match!" \
    #     || echo "md5sum of files DO NOT MATCH!"; exit 1
    # fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
