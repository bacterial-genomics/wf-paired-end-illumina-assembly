process CLASSIFY_CONTIGS_CAT {

    label "process_high"
    tag { "${meta.id}" }
    container "tpaisie/cat:latest"

    input:
    tuple val(meta), path(assembly)
    path(database)

    output:
    tuple val(meta), path("${meta.id}.CAT-Classification*tsv"), emit: summary
    path("${meta.id}.Contigs.tsv")                            , emit: output
    tuple val(meta), path("${meta.id}.CAT_Output_File.tsv")   , emit: qc_filecheck
    path("${meta.id}.CAT-Classification.log.gz")
    path(".command.{out,err}")
    path("versions.yml")                                      , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Classify FastA assembly contigs with CAT
    msg "INFO: Classifying each !{meta.id} assembly contig independently with CAT ..."

    # Run CAT
    CAT \
      contigs \
      --contigs_fasta "!{assembly}" \
      --database_folder db \
      --taxonomy_folder tax \
      --out_prefix "!{meta.id}.CAT-Classification" \
      --nproc !{task.cpus}

    msg "INFO: Completed classification of each !{meta.id} contig with CAT"

    # Verify output files
    msg "INFO: Creating QC summary file of CAT output files for !{meta.id} ..."

    CREATE_SUMMARY=true
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}.CAT_Output_File.tsv"
    for file in !{meta.id}.CAT-Classification.ORF2LCA.txt !{meta.id}.CAT-Classification.contig2classification.txt; do
      cat_file_type=$(echo "${file}" | cut -d '.' -f 3)
      if verify_minimum_file_size "${file}" 'CAT Output File' "!{params.min_filesize_cat_output}"; then
        echo -e "!{meta.id}\tCAT Output File (${cat_file_type})\tPASS" >> !{meta.id}.CAT_Output_File.tsv
      else
        echo -e "!{meta.id}\tCAT Output File (${cat_file_type})\tFAIL" >> !{meta.id}.CAT_Output_File.tsv
        CREATE_SUMMARY=false
      fi
    done

    msg "INFO: Completed QC summary file of CAT output files for !{meta.id}"

    # Add taxonomic names to the CAT output if the GTDB database is not used
    if [[ "${CREATE_SUMMARY}" ]] && \
      [[ ! $(grep "__" *.txt) ]]; then
      msg "INFO: Adding names to CAT ORF2LCA output file for !{meta.id} ..."

      CAT \
        add_names \
        --only_official \
        --input_file !{meta.id}.CAT-Classification.ORF2LCA.txt \
        --output_file !{meta.id}.CAT-Classification.ORF2LCA.names.tsv \
        --taxonomy_folder tax

      msg "INFO: Adding names to CAT contig2classification output file for !{meta.id} ..."

      CAT \
        add_names \
        --only_official \
        --input_file !{meta.id}.CAT-Classification.contig2classification.txt \
        --output_file !{meta.id}.CAT-Classification.names.tsv \
        --taxonomy_folder tax

      msg "INFO: Creating CAT summary file for !{meta.id} ..."
      CAT \
        summarise \
        --input_file !{meta.id}.CAT-Classification.names.tsv \
        --output_file !{meta.id}.CAT-Classification.names.summary.tsv \
        --contigs_fasta "!{assembly}"
      
      # Custom simplified TSV reports
      msg "INFO: Creating custom simpler CAT summary file for !{meta.id} ..."
      grep '^# rank' !{meta.id}.CAT-Classification.names.summary.tsv \
        | sed "s/^# /Sample_name\t/1;s/ /_/g" \
        > Contigs.header_line.tsv
      grep -v -e "no support" -e "^#" !{meta.id}.CAT-Classification.names.summary.tsv \
        | awk -v var=!{meta.id} '{print var "\t" $0}' \
        > Contigs.only-supported-data.tsv
      cat Contigs.header_line.tsv \
        Contigs.only-supported-data.tsv \
        > !{meta.id}.Contigs.tsv

      grep '^# ORF' !{meta.id}.CAT-Classification.ORF2LCA.names.tsv \
        | sed "s/^# /Sample_name\t/1;s/ /_/g" \
        | awk -v var=!{meta.id} '{print var "\t" $0}' \
        > ORF.header_line.unique-lineages.tsv
      grep -v -e $'no support\tno support' -e "^#" !{meta.id}.CAT-Classification.ORF2LCA.names.tsv \
        | sort -u -k3,3 | sed "s/ /_/g" > ORF.only-supported-data.unique-lineages.tsv
      cat ORF.header_line.unique-lineages.tsv \
        ORF.only-supported-data.unique-lineages.tsv \
        > !{meta.id}.unique-lineages.ORF.tsv

      msg "INFO: Completed CAT reports for !{meta.id}"

    else
      msg "WARN: Skipping addition of names to CAT ORF2LCA output file for !{meta.id} ..."

      # Avoid "no output file found" error if name/summary tsv files are not created
      cp "!{meta.id}.CAT-Classification.ORF2LCA.txt" "!{meta.id}.CAT-Classification.ORF2LCA.tsv"
      cp "!{meta.id}.CAT-Classification.contig2classification.txt" "!{meta.id}.CAT-Classification.contig2classification.tsv"
      touch !{meta.id}.Contigs.tsv
    fi

    # Compress the bulky verbose logfile for compact storage
    gzip -9f !{meta.id}.CAT-Classification.log

    msg "INFO: Completed CAT classification process for !{meta.id}"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        cat: $(CAT --version | sed 's/ by F\\. A\\. Bastiaan.*//1' | sed 's/^CAT //1')
    END_VERSIONS
    '''
}
