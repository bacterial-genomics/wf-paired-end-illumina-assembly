process ANNOTATE_PROKKA {

    label "process_high"
    tag { "${meta.id}" }
    container "snads/prokka@sha256:ef7ee0835819dbb35cf69d1a2c41c5060691e71f9138288dd79d4922fa6d0050"

    input:
    tuple val(meta), path(assembly)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                         , emit: versions
    path "${meta.id}.Annotated_GenBank_File.tsv", emit: qc_annotated_filecheck
    tuple val(meta), path("${meta.id}.gbk")     , emit: annotation

    shell:
    '''
    source bash_functions.sh

    # Remove seperator characters from basename for future processes
    short_base=$(echo !{meta.id} | sed 's/[-._].*//g')
    sed -i "s/!{meta.id}/${short_base}/g" !{assembly}

    # Annotate cleaned and corrected assembly
    msg "INFO: Annotating assembly using Prokka"

    prokka \
      --outdir prokka \
      --prefix "!{meta.id}" \
      --force \
      --addgenes \
      --locustag "!{meta.id}" \
      --mincontiglen 1 \
      --evalue 1e-08 \
      --cpus !{task.cpus} \
      !{assembly}

    for ext in gb gbf gbff gbk; do
      if [ -s "prokka/!{meta.id}.${ext}" ]; then
        mv -f prokka/!{meta.id}.${ext} !{meta.id}.gbk
        break
      fi
    done

    if verify_minimum_file_size "!{meta.id}.gbk" 'Annotated GenBank File' "!{params.min_filesize_annotated_genbank}"; then
      echo -e "!{meta.id}\tAnnotated GenBank File\tPASS" \
      > !{meta.id}.Annotated_GenBank_File.tsv
    else
      echo -e "!{meta.id}\tAnnotated GenBank File\tFAIL" \
      > !{meta.id}.Annotated_GenBank_File.tsv
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        prokka: $(prokka --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
