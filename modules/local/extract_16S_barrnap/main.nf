process EXTRACT_16S_BARRNAP {

    tag { "${meta.id}" }
    container "snads/barrnap@sha256:e22cbd789c36d5626460feb6c7e5f6f7d55c8628dacae68ba0da30884195a837"

    input:
    tuple val(meta), path(annotation), path(assembly), path(extracted_rna)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                , emit: versions
    path("${meta.id}.SSU_{Renamed,Extracted}_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("16S.${meta.id}.fa")         , emit: extracted_base

    shell:
    '''
    source bash_functions.sh

    if [[ ! -f "!{extracted_rna}" ]] || [[ ! -s "!{extracted_rna}" ]]; then
      msg "INFO: Absent 16S rRNA gene annotation in !{annotation}" >&2
      msg 'Running barrnap' >&2

      barrnap !{assembly} | grep "Name=16S_rRNA;product=16S" > !{meta.id}.gff

      if [[ $(grep -c "Name=16S_rRNA;product=16S" "!{meta.id}.gff") -eq 0 ]]; then
        msg "INFO: Barrnap was unable to locate a 16S rRNA gene sequence in !{assembly}" >&2
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
