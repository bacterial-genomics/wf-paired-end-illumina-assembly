process ASSESS_ASSEMBLY_CHECKM2 {

    label "process_high"  //NOTE: h_vmem=126.6G+ and h_rss=103114M+ normally required or else exit 140 status code
    tag { "${meta.id}" }
    // container "quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0"
    container "quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0"

    input:
    tuple val(meta), path(assembly)
    path(database)

    output:
    tuple val(meta), path("${meta.id}.CheckM2_Report_File.tsv"),               emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.CheckM2.results.tsv"), emit: summary
    path("${meta.id}-${meta.assembler}.CheckM2.log.gz")
    path("${meta.id}-${meta.assembler}.CheckM2.alignments.tsv.gz")
    path(".command.{out,err}")
    path("versions.yml"),                                                      emit: versions

    shell:
    '''
    source bash_functions.sh

    # Assess the full FastA assembly with CheckM2
    msg "INFO: Evaluating the assembly contig set of !{meta.id} for completeness and contamination with CheckM2"

    # Run CheckM2
    # NOTE: h_vmem=126.6G+ and h_rss=103114M+ normally required or else exit 140 status code
    checkm2 \
      predict \
      --input !{assembly} \
      --output-directory checkm2 \
      --database_path !{database} \
      --force \
      !{params.checkm2_model} \
      --threads !{task.cpus}

    msg "INFO: CheckM2 completed for !{meta.id} assembly"

    # Move and rename the report, alignments, and log files
    mv -f checkm2/quality_report.tsv "!{meta.id}-!{meta.assembler}.CheckM2.results.tsv"
    mv -f checkm2/diamond_output/DIAMOND_RESULTS.tsv "!{meta.id}-!{meta.assembler}.CheckM2.alignments.tsv"
    mv -f checkm2/checkm2.log "!{meta.id}-!{meta.assembler}.CheckM2.log"

    # Test/verify paired FastQ outfiles sizes are reasonable to continue
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.CheckM2_Report_File.tsv"
    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.CheckM2.results.tsv" 'CheckM2 Report File' "!{params.min_filesize_checkm2_report}"; then
      echo -e "!{meta.id}\tCheckM2 Report File\tPASS" >> "!{meta.id}.CheckM2_Report_File.tsv"
    else
      echo -e "!{meta.id}\tCheckM2 Report File\tFAIL" >> "!{meta.id}.CheckM2_Report_File.tsv"
    fi

    # Replace space characters in header line with underscores
    sed -i '1s/ /_/g' "!{meta.id}-!{meta.assembler}.CheckM2.results.tsv"
    sed -i '1s/^Name/Sample_name/1' "!{meta.id}-!{meta.assembler}.CheckM2.results.tsv"

    # Replace id-assembler with just id in the data row
    sed -i "2s/!{meta.id}-!{meta.assembler}/!{meta.id}/1" "!{meta.id}-!{meta.assembler}.CheckM2.results.tsv"

    # Compress the log and alignments files for compact storage
    gzip -9f "!{meta.id}-!{meta.assembler}.CheckM2.log" "!{meta.id}-!{meta.assembler}.CheckM2.alignments.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        checkm2: $(checkm2 --version)
    END_VERSIONS
    '''
}
