process ASSESS_ASSEMBLY_CHECKM2 {

    label "process_high"
    tag { "${meta.id}" }
    container "quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0"

    input:
    tuple val(meta), path(assembly)
    path(database)

    output:
    tuple val(meta), path("${meta.id}.CheckM2_Report_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}.CheckM2.report.tsv")     , emit: summary
    path("${meta.id}.CheckM2.log.gz")
    path("${meta.id}.CheckM2.alignments.tsv.gz")
    path(".command.{out,err}")
    path("versions.yml")                                       , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Assess the full FastA assembly with CheckM2
    msg "INFO: Evaluating the assembly contig set for completeness and contamination with CheckM2"

    # Run CheckM2
    checkm2 \
      predict \
      --input !{assembly} \
      --output-directory checkm2 \
      --database_path !{database} \
      --force \
      !{params.checkm2_model} \
      --threads !{task.cpus}

    # Move and rename the report, alignments, and log files
    mv -f checkm2/quality_report.tsv "!{meta.id}.CheckM2.report.tsv"
    mv -f checkm2/diamond_output/DIAMOND_RESULTS.tsv "!{meta.id}.CheckM2.alignments.tsv"
    mv -f checkm2/checkm2.log "!{meta.id}.CheckM2.log"

    # Verify output report file
    if verify_minimum_file_size "!{meta.id}.CheckM2.report.tsv" 'CheckM2 Report File' "!{params.min_filesize_checkm2_report}"; then
      echo -e "!{meta.id}\tCheckM2 Report File\tPASS" > !{meta.id}.CheckM2_Report_File.tsv
    else
      echo -e "!{meta.id}\tCheckM2 Report File\tFAIL" > !{meta.id}.CheckM2_Report_File.tsv
    fi

    # Replace space characters in header line with underscores
    sed -i '1s/ /_/g' !{meta.id}.CheckM2.report.tsv

    # Also replace data content spaces with underscores as workaround issue with pandas converting to XLSX, error message:
    # "pandas.errors.ParserError: Error tokenizing data. C error: Expected 9 fields in line 3, saw 11" where #s can be different, so
    # as a workaround simply replace them with underscores (e.g., "Neural Network (Specific Model)" -> "Neural_Network_(Specific_Model)")
    sed -i '2s/ /_/g' !{meta.id}.CheckM2.report.tsv

    # Compress the log and alignments files for compact storage
    gzip -9f "!{meta.id}.CheckM2.log" "!{meta.id}.CheckM2.alignments.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        checkm2: $(checkm2 --version)
    END_VERSIONS
    '''
}
