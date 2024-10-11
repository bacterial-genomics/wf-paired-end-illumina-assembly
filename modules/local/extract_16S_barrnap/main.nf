process EXTRACT_16S_BARRNAP {

    tag { "${meta.id}-${meta.assembler}" }
    container "snads/barrnap@sha256:e22cbd789c36d5626460feb6c7e5f6f7d55c8628dacae68ba0da30884195a837"

    input:
    tuple val(meta), path(assembly), path(biopython_extracted_rna)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.SSU_{Renamed,Extracted}_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("16S.${meta.id}-${meta.assembler}.fa")                          , emit: extracted_rna
    path(".command.{out,err}")
    path("versions.yml")                                                                  , emit: versions

    shell:
    '''
    source bash_functions.sh

    if [[ ! -s "!{biopython_extracted_rna}" ]]; then
      msg "INFO: Absent 16S rRNA gene annotation in !{meta.id}-!{meta.assembler}.gbk" >&2
      msg 'Running barrnap' >&2

      rm "!{biopython_extracted_rna}"

      barrnap \
        --reject 0.1 \
        !{assembly} \
        | grep "Name=16S_rRNA;product=16S" \
        > "!{meta.id}-!{meta.assembler}.gff" \
        || :

      if [[ $(grep -c "Name=16S_rRNA;product=16S" "!{meta.id}-!{meta.assembler}.gff") -eq 0 ]]; then
        msg "INFO: Barrnap was unable to locate a 16S rRNA gene sequence in !{assembly}" >&2
        touch "16S.!{meta.id}-!{meta.assembler}.fa"
      else
        bedtools getfasta \
          -fi !{assembly} \
          -bed "!{meta.id}-!{meta.assembler}.gff" \
          -fo "16S.!{meta.id}-!{meta.assembler}.fa"
      fi
    fi

    echo -e "Sample _name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.SSU_Extracted_File.tsv"
    if verify_minimum_file_size "16S.!{meta.id}-!{meta.assembler}.fa" 'SSU Extracted FastA File' "!{params.min_filesize_extracted_ssu_file}"; then
      echo -e "!{meta.id}\tSSU Extracted FastA File\tPASS"  \
        >> "!{meta.id}.SSU_Extracted_File.tsv"
    else
      echo -e "!{meta.id}\tSSU Extracted FastA File\tFAIL" \
        >> "!{meta.id}.SSU_Extracted_File.tsv"
    fi

    awk -v awk_var="!{meta.id}" \
      '/^>/{print ">" awk_var "_" ++i; next} {print}' \
      "16S.!{meta.id}-!{meta.assembler}.fa" \
      > "!{meta.id}-!{meta.assembler}.fa-renamed"

    rm -f "16S.!{meta.id}-!{meta.assembler}.fa"
    mv -f "!{meta.id}-!{meta.assembler}.fa-renamed" \
      "16S.!{meta.id}-!{meta.assembler}.fa"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.SSU_Renamed_File.tsv"
    if verify_minimum_file_size "16S.!{meta.id}-!{meta.assembler}.fa" 'SSU Renamed FastA File' "!{params.min_filesize_renamed_ssu_file}"; then
      echo -e "!{meta.id}\tSSU Renamed FastA File\tPASS"  \
        >> "!{meta.id}.SSU_Renamed_File.tsv"
    else
      echo -e "!{meta.id}\tSSU Renamed FastA File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.SSU_Renamed_File.tsv"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        barrnap: $(barrnap --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
