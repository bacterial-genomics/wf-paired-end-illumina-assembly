process ANNOTATE_PROKKA {

    publishDir "${params.outdir}/annot",
        mode: "${params.publish_dir_mode}",
        pattern: "*.gbk"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Annotated_GenBank_File.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    label "process_high"
    tag { "${prefix}" }

    container "snads/prokka@sha256:ef7ee0835819dbb35cf69d1a2c41c5060691e71f9138288dd79d4922fa6d0050"

    input:
    tuple val(prefix), path(paired_bam), path(single_bam), path(qc_assembly_filecheck), path(assembly)

    output:
    tuple val(prefix), path("${prefix}.gbk"), path("*File*.tsv"), emit: annotation
    path "${prefix}.Annotated_GenBank_File.tsv", emit: qc_annotated_filecheck
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_assembly_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done
    
    # Remove seperator characters from basename for future processes
    short_base=$(echo !{prefix} | sed 's/[-._].*//g')
    sed -i "s/!{prefix}/${short_base}/g" !{assembly}

    # Annotate cleaned and corrected assembly
    msg "INFO: Running prokka with !{task.cpus} threads"

    prokka \
      --outdir prokka \
      --prefix "!{prefix}"\
      --force \
      --addgenes \
      --locustag "!{prefix}" \
      --mincontiglen 1 \
      --evalue 1e-08 \
      --cpus !{task.cpus} \
      !{assembly}

    for ext in gb gbf gbff gbk; do
      if [ -s "prokka/!{prefix}.${ext}" ]; then
        mv -f prokka/!{prefix}.${ext} !{prefix}.gbk
        break
      fi
    done

    if verify_minimum_file_size "!{prefix}.gbk" 'Annotated GenBank File' "!{params.min_filesize_annotated_genbank}"; then
      echo -e "!{prefix}\tAnnotated GenBank File\tPASS" \
        > !{prefix}.Annotated_GenBank_File.tsv
    else
      echo -e "!{prefix}\tAnnotated GenBank File\tFAIL" \
        > !{prefix}.Annotated_GenBank_File.tsv
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      prokka: $(prokka --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}