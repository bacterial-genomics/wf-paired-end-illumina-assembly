process BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(blast_output)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Filtered_16S_BLASTn_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.16S-top-species.tsv")         , emit: top_blast_species
    path("${meta.id}-${meta.assembler}.blast.tsv.gz")
    path(".command.{out,err}")
    path("versions.yml")                                                              , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Get the top match by bitscore
    filter.blast.py \
      -i "!{blast_output}" \
      -o "!{meta.id}-!{meta.assembler}.blast.tsv" \
      -c !{params.filter_blast_column} \
      -s !{params.filter_blast_bitscore}

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.blast.tsv" 'Filtered 16S BLASTn File' "!{params.min_filesize_filtered_blastn}"; then
      echo -e "!{meta.id}-!{meta.assembler}\tFiltered 16S BLASTn File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

      # Report the top alignment match data: %nucl iden, %query cov aln, taxon
      awk -F $'\t' 'BEGIN{OFS=FS}; {print $1, $3 "% identity", $13 "% alignment", $14}' \
      "!{meta.id}-!{meta.assembler}.blast.tsv" \
      > "!{meta.id}-!{meta.assembler}.16S-top-species.tsv"

      sed -i \
        '1i Sample name\tPercent identity\tPercent alignment\tSpecies match' \
        "!{meta.id}-!{meta.assembler}.16S-top-species.tsv"

    else
      echo -e "!{meta.id}-!{meta.assembler}\tFiltered 16S BLASTn File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

      # Empty files to avoid errors - use qcfilecheck function to end workflow
      touch "!{meta.id}-!{meta.assembler}.16S-top-species.tsv"
    fi

    # Add header to BLAST output
    SUMMARY_HEADER=(
      "Query ID"
      "Reference ID"
      "Identity (%)"
      "Alignment length"
      "Number of mismatches"
      "Number of gap openings"
      "Query start position"
      "Query end position"
      "Reference start position"
      "Reference end position"
      "Expect value"
      "Bit score"
      "Query coverage per HSP"
      "Reference scientific name"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}")
    sed -i "1i ${SUMMARY_HEADER}" "!{blast_output}"

    gzip -f !{blast_output}

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
