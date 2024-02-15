process CLASSIFY_ASSEMBLY_CHECKM2 {

    label "process_high"
    tag { "${meta.id}" }
    container "quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0@sha256:sha256:f4adc81bff88ab5a27a2a7c4e7af2cdb0943a7b89e76ef9d2f7ec680a3b95111"

    input:
    tuple val(meta), path(assembly)
    path(database)

    output:
    tuple val(meta), path("${meta.id}.CheckM2_Report_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}.CheckM2.report.tsv")     , emit: summary
    path("${meta.id}.CheckM2.log.gz")
    path(".command.{out,err}")
    path("versions.yml")                                       , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Assess the full FastA assembly with CheckM2
    msg "INFO: Classifying assembly contig set with CheckM2"

    # Run CheckM2
    checkm2 \
      predict \
      --database_path !{database} \
      --input !{assembly} \
      --output-directory checkm2 \
      --force \
      !{params.checkm2_model} \
      --threads !{task.cpus}

    # Move and rename report and log
    mv -f checkm2/quality_report.tsv "!{meta.id}.CheckM2.report.tsv"
    mv -f checkm2/checkm2.log "!{meta.id}.CheckM2.log"

    # Verify output file
    if verify_minimum_file_size "!{meta.id}.CheckM2.report.tsv" 'CheckM2 Report File' "!{params.min_filesize_checkm2_report}"; then
      echo -e "!{meta.id}\tCheckM2 Report File\tPASS" > !{meta.id}.CheckM2_Report_File.tsv
    else
      echo -e "!{meta.id}\tCheckM2 Report File\tFAIL" > !{meta.id}.CheckM2_Report_File.tsv
    fi

    # Compress the logfile for compact storage
    gzip -9f "!{meta.id}.CheckM2.log"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        checkm2: $(checkm2 --version)
    END_VERSIONS
    '''
}
