process CONTIG_CLASSIFY_CAT {

    label "process_high"
    tag { "${meta.id}" }
    container "biocontainers/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-0"

    input:
    tuple val(meta), path(assembly)
    path alignment_database
    path taxonomy_database

    output:
    path(".command.{out,err}")
    path "cat-contigs.${meta.id}.log.gz"
    path "versions.yml"                                        , emit: versions
    tuple val(meta), path("cat-contigs.${meta.id}.*.names*tsv"), emit: cat_classified_files

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
      --out_prefix cat-contigs."!{meta.id}" \
      --nproc !{task.cpus}

    # Verify output files
    for file in cat-contigs.!{meta.id}.ORF2LCA.txt cat-contigs.!{meta.id}.contig2classification.txt; do
      if verify_minimum_file_size "${file}" 'CAT Output File' '1c'; then
        continue
      else
        msg "ERROR: ${file} CAT output file missing" >&2
        exit 1
      fi
    done

    # Add taxonomic names to the CAT output
    CAT \
      add_names \
      --only_official \
      --input_file cat-contigs.!{meta.id}.ORF2LCA.txt \
      --output_file cat-contigs.!{meta.id}.ORF2LCA.names.tsv \
      --taxonomy_folder !{taxonomy_database}
    CAT \
      add_names \
      --only_official \
      --input_file cat-contigs.!{meta.id}.contig2classification.txt \
      --output_file cat-contigs.!{meta.id}.contig2classification.names.tab \
      --taxonomy_folder !{taxonomy_database}

    # Create CAT summary file for contigs
    CAT \
      summarise \
      --input_file cat-contigs.!{meta.id}.contig2classification.names.tab \
      --output_file cat-contigs.!{meta.id}.contig2classification.names.summary.tsv \
      --contigs_fasta "!{assembly}"

    # Compress the bulky verbose logfile for compact storage
    gzip -9f cat-contigs.!{meta.id}.log

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        cat: $(CAT --version | sed 's/ by F\. A\. Bastiaan.*//1' | sed 's/^CAT //1')
    END_VERSIONS
    '''
}
