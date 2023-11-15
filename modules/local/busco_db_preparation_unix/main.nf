process BUSCO_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${database.getSimpleName()}" }
    container "ubuntu:jammy"

    input:
    path(database)

    output:
    path(".command.{out,err}")
    path "versions.yml", emit: versions
    path("${dir}")   , emit: db

    shell:
    db_name = database.getSimpleName()
    dir = db_name.contains('odb10') ? "BUSCO/lineages/${db_name}" : db_name
    '''
    source bash_functions.sh

    mkdir -p !{dir}
    tar -xzf !{database} -C !{dir} --strip-components 1

    # Check for `info` and `hmms` directories
    # Check lineage dataset
    if [[ !{db_name} =~ odb10 ]]; then
      for directory in info hmms; do
        if [[ ! -d "!{dir}/${directory}" ]]; then
          msg "ERROR: BUSCO dataset is missing required directory: `${directory}`."
          exit 1
        fi
      done
    else
      # Check if larger BUSCO database
      num_odb10_dirs=$(find !{dir}/ -maxdepth 2 -type d -name "*_odb10" | wc -l)
      num_hmms_dirs=$(find !{dir}/ -maxdepth 3 -type d -name "hmms" | wc -l)
      num_info_dirs=$(find !{dir}/ -maxdepth 3 -type d -name "info" | wc -l)

      if [[ $num_odb10_dirs != $num_hmms_dirs ]] && [[ $num_odb10_dirs != $num_info_dirs ]]; then
        msg "ERROR: BUSCO database does not have the required directories `hmms` and `info` in each lineage dataseet."
        exit 1
      fi
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
