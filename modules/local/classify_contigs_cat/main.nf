process CLASSIFY_CONTIGS_CAT {

    label "process_high"
    tag { "${meta.id}" }
    container "biocontainers/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-0"

    input:
    tuple val(meta), path(assembly)
    path(alignment_database)
    path(taxonomy_database)

    output:
    tuple val(meta), path("${meta.id}.CAT-Classification*tsv"), emit: summary
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
      --database_folder !{alignment_database} \
      --taxonomy_folder !{taxonomy_database} \
      --out_prefix "!{meta.id}.CAT-Classification" \
      --nproc !{task.cpus}

    # Verify output files
    CREATE_SUMMARY=true
    for file in !{meta.id}.CAT-Classification.ORF2LCA.txt !{meta.id}.CAT-Classification.contig2classification.txt; do
      cat_file_type=$(echo ${file} | cut -d '.' -f3)
      if verify_minimum_file_size "${file}" 'CAT Output File' "!{params.min_filesize_cat_output}"; then
        echo -e "!{meta.id}\tCAT Output File (${cat_file_type})\tPASS" >> !{meta.id}.CAT_Output_File.tsv
      else
        echo -e "!{meta.id}\tCAT Output File (${cat_file_type})\tFAIL" >> !{meta.id}.CAT_Output_File.tsv
        CREATE_SUMMARY=false
      fi
    done

    # Add taxonomic names to the CAT output
    if [[ $CREATE_SUMMARY ]]; then
      CAT \
        add_names \
        --only_official \
        --input_file !{meta.id}.CAT-Classification.ORF2LCA.txt \
        --output_file !{meta.id}.CAT-Classification.ORF2LCA.names.tsv \
        --taxonomy_folder !{taxonomy_database}

      CAT \
        add_names \
        --only_official \
        --input_file !{meta.id}.CAT-Classification.contig2classification.txt \
        --output_file !{meta.id}.CAT-Classification.names.tsv \
        --taxonomy_folder !{taxonomy_database}

      # Create CAT summary file for contigs
      CAT \
        summarise \
        --input_file !{meta.id}.CAT-Classification.names.tsv \
        --output_file !{meta.id}.CAT-Classification.names.summary.tsv \
        --contigs_fasta "!{assembly}"
    else
      # Avoid "not output file found" error if name/summary tsv files are not created
      cp "!{meta.id}.CAT-Classification.ORF2LCA.txt" "!{meta.id}.CAT-Classification.ORF2LCA.tsv"
      cp "!{meta.id}.CAT-Classification.contig2classification.txt" "!{meta.id}.CAT-Classification.contig2classification.tsv"
    fi

    # Compress the bulky verbose logfile for compact storage
    gzip -9f !{meta.id}.CAT-Classification.log

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        cat: $(CAT --version | sed 's/ by F\. A\. Bastiaan.*//1' | sed 's/^CAT //1')
    END_VERSIONS
    '''
}
