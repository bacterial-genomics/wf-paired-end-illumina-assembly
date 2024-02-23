process DOWNLOAD_CAT_DB_UNIX {

    label "process_medium"
    tag { "${meta.id}" }
    container "quay.io/biocontainers/gnu-wget@sha256:a7253d598eb7aa4a826f56b5da6079f1176d2c8cad584ee1a6c06f9f9b25d8b0"

    input:
    tuple val(meta)

    output:
    path("database/{db,tax}") , emit: db
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    # Download webpage and parse HTML
    wget -qSO - https://tbb.bio.uu.nl/tina/CAT_prepare/ \
      | grep -o "<a href=['"'"'"][^'"'"'"]*" \
      | sed 's|<a href="||g' \
      | grep "tar.gz" \
      | grep "nr" \
      | while read fname; do
          wget -c https://tbb.bio.uu.nl/tina/CAT_prepare/${fname}
        done

    # Verify filename in md5 file matches downloaded database
    if [[ ! -f $(find . -name $(awk '{print $2}' *.tar.gz.md5 )) ]]; then
      echo "ERROR: Database file not found!"
      exit 1
    else
      [[ $(md5sum *.tar.gz | awk '{print $1}') = $(awk '{print $1}' *.tar.gz.md5) ]] \
        && echo "md5sum of files match!"; tar -xzf CAT*.tar.gz \
        || echo "md5sum of files DO NOT MATCH!"; exit 1
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
