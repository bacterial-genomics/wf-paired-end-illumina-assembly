process BUSCO_DB_PREPARATION_UNIX {

    label "process_medium"
    tag { "${meta.id}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(database)

    output:
    path("${output_dir}")     , emit: db
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    output_dir = meta['id'].contains('odb10') ? "BUSCO/lineages/${meta.id}" : "${meta.id}"
    '''
    source bash_functions.sh

    mkdir -p !{output_dir}
    tar -xzf !{database} -C !{output_dir} --strip-components 1

    # Check for `info` and `hmms` directories
    # Check lineage dataset
    if [[ !{meta.id} =~ odb10 ]]; then
      for directory in info hmms; do
        if [[ ! -d "!{output_dir}/${directory}" ]]; then
          msg "ERROR: BUSCO dataset is missing required directory: `${directory}`." >&2
          exit 1
        fi
      done
    else
      # Check if larger BUSCO database
      num_odb10_dirs=$(find !{output_dir}/ -maxdepth 2 -type d -name "*_odb10" | wc -l)
      num_hmms_dirs=$(find !{output_dir}/ -maxdepth 3 -type d -name "hmms" | wc -l)
      num_info_dirs=$(find !{output_dir}/ -maxdepth 3 -type d -name "info" | wc -l)

      if [[ $num_odb10_dirs != $num_hmms_dirs ]] && [[ $num_odb10_dirs != $num_info_dirs ]]; then
        msg "ERROR: BUSCO database does not have the required directories `hmms` and `info` in each lineage dataseet." >&2
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
