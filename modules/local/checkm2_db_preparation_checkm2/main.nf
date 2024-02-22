process CHECKM2_DB_PREPARATION_CHECKM2 {

    label "process_medium"
    tag { "${meta.id}" }
    container "tpaisie/checkm2:latest"

    input:
    val(meta)

    output:
    path("database/*.dmnd")   , emit: db
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    mkdir -p database

    #checkm2 database --download
    #mv -fv $(find . -name "*.dmnd") database/

    # NOTE: checkm2 downloads a single *.dmnd file but built-in feature might
    #       fail due error:
    #       "urllib3.exceptions.SSLError: [SSL: CERTIFICATE_VERIFY_FAILED]"
    wget -c https://zenodo.org/records/5571251/files/checkm2_database.tar.gz
    tar -xvzf checkm2_database.tar.gz -C database/
    mv -fv $(find . -name "*.dmnd") database/

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
