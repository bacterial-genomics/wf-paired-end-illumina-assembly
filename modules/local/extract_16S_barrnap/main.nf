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
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    tag { "${prefix}" }

    container "snads/barrnap@sha256:e22cbd789c36d5626460feb6c7e5f6f7d55c8628dacae68ba0da30884195a837"

    input:
    tuple val(prefix), path(annotation), path(qc_annotated_filecheck), path(assembly), path(extracted_rna)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                            , emit: versions
    path "${prefix}.SSU_Renamed_File.tsv"                          , emit: qc_ssu_renamed_filecheck
    path "${prefix}.SSU_Extracted_File.tsv"                        , emit: qc_ssu_extracted_filecheck
    tuple val(prefix), path("16S.${prefix}.fa"), path("*File*.tsv"), emit: extracted_base

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

      barrnap !{assembly} | grep "Name=16S_rRNA;product=16S" > !{prefix}.gff

      if [[ $(grep -c "Name=16S_rRNA;product=16S" "!{prefix}.gff") -eq 0 ]]; then
        msg "INFO: barrnap was unable to locate a 16S rRNA gene sequence in !{assembly}" >&2
        exit 2
      fi

      bedtools getfasta \
        -fi !{assembly} \
        -bed !{prefix}.gff \
        -fo 16S.!{prefix}.fa
    fi

    if verify_minimum_file_size "16S.!{prefix}.fa" 'SSU Extracted File' "!{params.min_filesize_extracted_ssu_file}"; then
      echo -e "!{prefix}\tSSU Extracted File\tPASS" > !{prefix}.SSU_Extracted_File.tsv
    else
      echo -e "!{prefix}\tSSU Extracted File\tFAIL" > !{prefix}.SSU_Extracted_File.tsv
    fi

    awk -v awk_var="!{prefix}" \
      '/^>/{print ">" awk_var "_" ++i; next} {print}' \
      16S.!{prefix}.fa \
      > !{prefix}.fa-renamed
    rm -f 16S.!{prefix}.fa
    mv -f !{prefix}.fa-renamed 16S.!{prefix}.fa

    if verify_minimum_file_size "16S.!{prefix}.fa" 'SSU Renamed File' "!{params.min_filesize_renamed_ssu_file}"; then
      echo -e "!{prefix}\tSSU Renamed File\tPASS" > !{prefix}.SSU_Renamed_File.tsv
    else
      echo -e "!{prefix}\tSSU Renamed File\tFAIL" > !{prefix}.SSU_Renamed_File.tsv
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      barrnap: $(barrnap --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
