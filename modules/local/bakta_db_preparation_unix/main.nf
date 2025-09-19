process BAKTA_DB_PREPARATION_UNIX {

    label "process_low"
    tag { "${meta.id}" }
    // NOTE: bakta version (v1.11.0) here doesn't matter.
    // We just need a container with `xz` library and almost no others have it!
    container 'https://depot.galaxyproject.org/singularity/bakta:1.10.4--pyhdfd78af_0'

    input:
    tuple val(meta), path(database)

    output:
    path("found_db_dir/")     , emit: db
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    source bash_functions.sh

    if command -v xz >/dev/null 2>&1; then
      msg "INFO: xz is installed for bakta database uncompressing"
    else
      msg "ERROR: xz is not installed for bakta database uncompressing" >&2
    fi
    mkdir database
    tar -xf !{database} -C database

    # We need a bunch more files than this, but
    # at least confirm the largest one exists before moving on.
    if [ -s database/*/bakta.db ]; then
      msg "INFO: bakta database uncompressed"

      # Get the parent dirname for the database (e.g., database/db-light/)
      db=$(find database -type f -name "bakta.db" | head -n 1 | xargs dirname)
      mkdir found_db_dir
      cp -r "$db"/* found_db_dir/
      msg "INFO: using ${db} as the bakta database..."
    else
      msg "ERROR: bakta database lacking bakta.db" >&2
      exit 1
    fi

    # We could print the version number from version.json with jq

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
