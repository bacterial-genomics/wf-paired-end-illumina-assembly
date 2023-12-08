process ANNOTATE_PROKKA {

    label "process_high"
    tag { "${meta.id}-${meta.assembler}" }
    container "snads/prokka@sha256:ef7ee0835819dbb35cf69d1a2c41c5060691e71f9138288dd79d4922fa6d0050"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Annotated_GenBank_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.gbk")                       , emit: prokka_genbank_file
    path("prokka/${meta.id}-${meta.assembler}.log.gz")
    path(".command.{out,err}")
    path("versions.yml")                                                            , emit: versions

    shell:
    curated_proteins = params.prokka_curated_proteins ? "--proteins ${params.prokka_curated_proteins}" : ""
    '''
    source bash_functions.sh

    # Remove seperator characters from basename for future processes
    short_base=$(echo !{meta.id} | sed 's/[-._].*//g')
    sed -i "s/!{meta.id}/${short_base}/g" !{assembly}

    # Annotate cleaned and corrected assembly
    msg "INFO: Annotating assembly using Prokka"

    # Run Prokka
    prokka \
      --outdir prokka \
      --prefix "!{meta.id}-!{meta.assembler}" \
      --locustag "!{meta.id}-!{meta.assembler}" \
      --evalue !{params.prokka_evalue} \
      --mincontiglen 1 \
      --force \
      --addgenes \
      !{curated_proteins} \
      --cpus !{task.cpus} \
      !{assembly}

    # Regardless of the file extension, unify to GBK extension for GenBank format
    for ext in gb gbf gbff gbk; do
      if [ -s "prokka/!{meta.id}-!{meta.assembler}.${ext}" ]; then
        mv -f "prokka/!{meta.id}-!{meta.assembler}.${ext}" \
          "!{meta.id}-!{meta.assembler}.gbk"
        break
      fi
    done

    # Verify output file
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Annotated_GenBank_File.tsv"
    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.gbk" 'Annotated GenBank File' "!{params.min_filesize_annotated_genbank}"; then
      echo -e "!{meta.id}-!{meta.assembler}\tAnnotated GenBank File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Annotated_GenBank_File.tsv"
    else
      echo -e "!{meta.id}-!{meta.assembler}\tAnnotated GenBank File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Annotated_GenBank_File.tsv"
    fi

    # Compress the bulky verbose logfile for compact storage
    gzip -9f "prokka/!{meta.id}-!{meta.assembler}.log"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        prokka: $(prokka --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
