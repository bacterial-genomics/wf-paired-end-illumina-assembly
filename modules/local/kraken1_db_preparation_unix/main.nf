process KRAKEN1_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${meta.id}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(database)

    output:
    path("database/")         , emit: db
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    source bash_functions.sh

    mkdir -p db_tmp
    tar -xzf !{database} -C db_tmp/ --strip-components 1

    # Place kraken files in correct directory
    mkdir -p database/taxonomy
    mv `find db_tmp/ -name "database.idx" -o -name "database.kdb"` database/
    mv `find db_tmp/ -name "nodes.dmp" -o -name "names.dmp"` database/taxonomy/

    # Verify all 4 files are found
    if [[ $(find database/ -type f | wc -l) != 4 ]]; then
        msg "ERROR: Missing one of the following files: `database.{idx,kdb}, {names,nodes}.dmp`." >&2
        exit 1
    fi

    # ### Calculate SHA-512 Checksum Kraken inspect.txt file ###
    # SUMMARY_HEADER=(
    #   "Sample_name"
    #   "Checksum_(SHA-512)"
    #   "File"
    # )
    # SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')

    # echo "${SUMMARY_HEADER}" > "!{meta.id}.Kraken_Database.SHA512-checksums.tsv"

    # if [ -s "!{database}/inspect.txt" ]; then
    #     msg "INFO: Found pre-calculated inspect.txt Kraken db information"
    # else
    #     msg "INFO: Creating inspect.txt Kraken db information..."
    #     kraken2-inspect --db "!{database}" --threads "!{task.cpus}" > "!{database}/inspect.txt"
    #     msg "INFO: Creating inspect.txt Kraken db information..."
    # fi
    # CHECKSUM=$(sha512sum !{database}/inspect.txt | awk '{print $1}')
    # echo -e "!{meta.id}\t${CHECKSUM}\t!{database}/inspect.txt" >> "!{meta.id}.Kraken_Database.SHA512-checksums.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
