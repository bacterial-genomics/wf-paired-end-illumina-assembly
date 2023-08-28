process EXTRACT_16S_BARRNAP {

    publishDir "${params.outdir}/ssu",
        mode: "${params.publish_dir_mode}",
        pattern: "*.fa"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{SSU_Extracted_File,SSU_Renamed_File}.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${meta.id}.${task.process}${filename}" }

    tag { "${meta.id}" }

    container "snads/barrnap@sha256:e22cbd789c36d5626460feb6c7e5f6f7d55c8628dacae68ba0da30884195a837"

    input:
    tuple val(meta), path(annotation), path(qc_annotated_filecheck), path(assembly), path(extracted_rna)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                           , emit: versions
    path "${meta.id}.SSU_Renamed_File.tsv"                        , emit: qc_ssu_renamed_filecheck
    path "${meta.id}.SSU_Extracted_File.tsv"                      , emit: qc_ssu_extracted_filecheck
    tuple val(meta), path("16S.${meta.id}.fa"), path("*File*.tsv"), emit: extracted_base

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_annotated_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    if [[ ! -f "!{extracted_rna}" ]] || [[ ! -s "!{extracted_rna}" ]]; then
      msg "INFO: absent 16S rRNA gene annotation in !{annotation}" >&2
      msg 'Running barrnap' >&2

      barrnap !{assembly} | grep "Name=16S_rRNA;product=16S" > !{meta.id}.gff

      if [[ $(grep -c "Name=16S_rRNA;product=16S" "!{meta.id}.gff") -eq 0 ]]; then
        msg "INFO: barrnap was unable to locate a 16S rRNA gene sequence in !{assembly}" >&2
        exit 2
      fi

      bedtools getfasta \
        -fi !{assembly} \
        -bed !{meta.id}.gff \
        -fo 16S.!{meta.id}.fa
    fi

    if verify_minimum_file_size "16S.!{meta.id}.fa" 'SSU Extracted File' "!{params.min_filesize_extracted_ssu_file}"; then
      echo -e "!{meta.id}\tSSU Extracted File\tPASS" > !{meta.id}.SSU_Extracted_File.tsv
    else
      echo -e "!{meta.id}\tSSU Extracted File\tFAIL" > !{meta.id}.SSU_Extracted_File.tsv
    fi

    awk -v awk_var="!{meta.id}" \
      '/^>/{print ">" awk_var "_" ++i; next} {print}' \
      16S.!{meta.id}.fa \
      > !{meta.id}.fa-renamed
    rm -f 16S.!{meta.id}.fa
    mv -f !{meta.id}.fa-renamed 16S.!{meta.id}.fa

    if verify_minimum_file_size "16S.!{meta.id}.fa" 'SSU Renamed File' "!{params.min_filesize_renamed_ssu_file}"; then
      echo -e "!{meta.id}\tSSU Renamed File\tPASS" > !{meta.id}.SSU_Renamed_File.tsv
    else
      echo -e "!{meta.id}\tSSU Renamed File\tFAIL" > !{meta.id}.SSU_Renamed_File.tsv
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        barrnap: $(barrnap --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
