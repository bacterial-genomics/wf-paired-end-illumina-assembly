process BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(blast_output), path(assembly)

    output:
    path(".command.{out,err}")
    path "${meta.id}-${meta.assembler}.blast.tsv.gz"
    path "versions.yml"                                                               , emit: versions
    tuple val(meta), path("${meta.id}-${meta.assembler}.Summary.16S.tab")             , emit: blast_summary
    tuple val(meta), path("${meta.id}-${meta.assembler}.16S-top-species.tsv")         , emit: top_blast_species
    tuple val(meta), path("${meta.id}-${meta.assembler}.Filtered_16S_BLASTn_File.tsv"), emit: qc_filecheck

    shell:
    '''
    source bash_functions.sh

    # Get the top match by bitscore
    filter.blast.py \
      -i "!{blast_output}" \
      -o "!{meta.id}-!{meta.assembler}.blast.tab" \
      -c !{params.filter_blast_column} \
      -s !{params.filter_blast_bitscore}

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.blast.tab" 'Filtered 16S BLASTn File' "!{params.min_filesize_filtered_blastn}"; then
      echo -e "!{meta.id}\tFiltered 16S BLASTn File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

      # Report the top alignment match data: %nucl iden, %query cov aln, taxon
      awk -F $'\t' 'BEGIN{OFS=FS}; {print $1, $3 "% identity", $13 "% alignment", $14}' \
      "!{meta.id}-!{meta.assembler}.blast.tab" \
      > "!{meta.id}-!{meta.assembler}.16S-top-species.tsv"

      sed -i \
        '1i Sample name\tPercent identity\tPercent alignment\tSpecies match' \
        "!{meta.id}-!{meta.assembler}.16S-top-species.tsv"

      cat "!{meta.id}-!{meta.assembler}.16S-top-species.tsv" >> "!{meta.id}-!{meta.assembler}.Summary.16S.tab"

    else
      echo -e "!{meta.id}\tFiltered 16S BLASTn File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

      # Empty files to avoid errors
      touch "!{meta.id}-!{meta.assembler}.16S-top-species.tsv" \
        "!{meta.id}-!{meta.assembler}.Summary.16S.tab"
    fi

    gzip -f !{blast_output}

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
