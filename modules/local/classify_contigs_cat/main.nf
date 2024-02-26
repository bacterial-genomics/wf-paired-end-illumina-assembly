process CLASSIFY_CONTIGS_CAT {

    label "process_high"
    tag { "${meta.id}" }
    // container "quay.io/biocontainers/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-1"
    // container "biocontainers/cat:5.2.3--hdfd78af_1"
    container "tpaisie/cat:latest"

    input:
    tuple val(meta), path(assembly)
    path(database)

    output:
    tuple val(meta), path("${meta.id}.CAT-Classification*tsv"), emit: summary
    tuple val(meta), path("${meta.id}.CAT_Output_File.tsv")   , emit: qc_filecheck
    path("${meta.id}.CAT-Classification.log.gz")
    path(".command.{out,err}")
    path("versions.yml")                                      , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Classify FastA assembly contigs with CAT
    msg "INFO: Classifying each assembly contig independently with CAT"

    # Run CAT
    CAT \
      contigs \
      --contigs_fasta "!{assembly}" \
      --database_folder db \
      --taxonomy_folder tax \
      --out_prefix "!{meta.id}.CAT-Classification" \
      --nproc !{task.cpus}

    # Verify output files
    CREATE_SUMMARY=true
    for file in !{meta.id}.CAT-Classification.ORF2LCA.txt !{meta.id}.CAT-Classification.contig2classification.txt; do
      cat_file_type=$(echo "${file}" | cut -d '.' -f 3)
      if verify_minimum_file_size "${file}" 'CAT Output File' "!{params.min_filesize_cat_output}"; then
        echo -e "!{meta.id}\tCAT Output File (${cat_file_type})\tPASS" >> !{meta.id}.CAT_Output_File.tsv
      else
        echo -e "!{meta.id}\tCAT Output File (${cat_file_type})\tFAIL" >> !{meta.id}.CAT_Output_File.tsv
        CREATE_SUMMARY=false
      fi
    done

    # Add taxonomic names to the CAT output if the GTDB database is not used
    if [[ "${CREATE_SUMMARY}" ]] && \
      [[ ! $(grep "__" *.txt) ]]; then
      msg "INFO: Adding names to CAT ORF2LCA output file"
      CAT \
        add_names \
        --only_official \
        --input_file !{meta.id}.CAT-Classification.ORF2LCA.txt \
        --output_file !{meta.id}.CAT-Classification.ORF2LCA.names.tsv \
        --taxonomy_folder tax

      msg "INFO: Adding names to CAT contig2classification output file"
      CAT \
        add_names \
        --only_official \
        --input_file !{meta.id}.CAT-Classification.contig2classification.txt \
        --output_file !{meta.id}.CAT-Classification.names.tsv \
        --taxonomy_folder tax

      msg "INFO: Creating CAT summary file"
      CAT \
        summarise \
        --input_file !{meta.id}.CAT-Classification.names.tsv \
        --output_file !{meta.id}.CAT-Classification.names.summary.tsv \
        --contigs_fasta "!{assembly}"
    else
      # Avoid "no output file found" error if name/summary tsv files are not created
      cp "!{meta.id}.CAT-Classification.ORF2LCA.txt" "!{meta.id}.CAT-Classification.ORF2LCA.tsv"
      cp "!{meta.id}.CAT-Classification.contig2classification.txt" "!{meta.id}.CAT-Classification.contig2classification.tsv"
    fi

    # Compress the bulky verbose logfile for compact storage
    gzip -9f !{meta.id}.CAT-Classification.log

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        cat: $(CAT --version | sed 's/ by F\\. A\\. Bastiaan.*//1' | sed 's/^CAT //1')
    END_VERSIONS
    '''
}
