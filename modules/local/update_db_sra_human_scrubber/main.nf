process UPDATE_DB_SRA_HUMAN_SCRUBBER {

    label "process_low"
    tag { "${database}" }
    container "quay.io/biocontainers/sra-human-scrubber@sha256:2f6b6635af9ba3190fc2f96640b21f0285483bd1f50d6be229228c52fb747055"

    input:
    val(database)

    output:
    path("data/human_filter.db"), emit: db
    path(".command.{out,err}")
    path("versions.yml")        , emit: versions

    shell:
    '''
    source bash_functions.sh

    # NOTE: For now, only support the pre-formatted human db file NCBI
    #       provides. To create a DB file requires NCBI's STAT package,
    #       which is not bundled in the scrubber utility.
    #       https://github.com/ncbi/sra-human-scrubber/issues/20#issuecomment-1414392052

    # This shell script by NCBI is a bit weird. It needs the "./data/"
    #   subdir available (won't make it) to download the *.db file.
    # Without any arguments, it fetches the latest human scrub DB file,
    #   which contains a date, then makes a generalized symlink, e.g.,
    #   .
    #   └── data
    #       ├── human_filter.db -> human_filter.db.20240718v2
    #       └── human_filter.db.20240718v2

    mkdir -p data
    init_db.sh

    # Make sure we fetched the *.db file
    if ! verify_minimum_file_size "data/human_filter.db" 'SRA Human Scrubber DB file' "!{params.min_filesize_sra_human_scrubber_db_file}"; then
        msg "ERROR: Missing human_filter.db file for SRA Human Scrubber" >&2
        exit 1
    fi

    # Get process version information
    # NOTE: currently no option to print the software version number, but
    #       track this issue https://github.com/ncbi/sra-human-scrubber/issues/28
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sra-human-scrubber: 2.2.1
    END_VERSIONS
    '''
}
