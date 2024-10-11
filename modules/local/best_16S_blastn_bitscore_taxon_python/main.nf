process BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(blast_output)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Filtered_16S_BLASTn_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.16S-top-species.tsv")         , emit: summary
    path("${meta.id}-${meta.assembler}.blast.tsv.gz")
    path(".command.{out,err}")
    path("versions.yml")                                                              , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: filtering !{blast_output} for highest bitscore value to report top 16S BLASTn species match"

    # Get the top match by bitscore
    filter.blast.py \
      -i "!{blast_output}" \
      -o "!{meta.id}-!{meta.assembler}.top-blast-bitscore.tsv" \
      -c !{params.filter_blast_column} \
      -s !{params.filter_blast_bitscore}

    msg "INFO: top 16S rRNA gene match to species by BLAST created: !{meta.id}-!{meta.assembler}.top-blast-bitscore.tsv"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.top-blast-bitscore.tsv" 'Filtered 16S BLASTn File' "!{params.min_filesize_filtered_blastn}"; then
      echo -e "!{meta.id}-!{meta.assembler}\tFiltered 16S BLASTn File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

      # Report the top alignment match data: %nucl iden, %query cov aln, taxon
      #   and split up the first column "<Sample_name>_<int>"; add header.
      # NOTE: a[length(a)] is used to take the last item in cases where
      #       samplename is Name_S1_L001_1 so it would get the final "1"
      awk 'BEGIN { FS=OFS="\t"; print "Sample_name\tUnique_16S_rRNA_extraction_count\tIdentity_(%)\tAlignment_(%)\tSpecies_match" }
        { split($1, a, "_"); print a[1], a[length(a)], $3, $13, $14 }' \
        "!{meta.id}-!{meta.assembler}.top-blast-bitscore.tsv" \
        > tmp.tsv \
        && mv -f tmp.tsv "!{meta.id}-!{meta.assembler}.16S-top-species.tsv"

    else
      echo -e "!{meta.id}-!{meta.assembler}\tFiltered 16S BLASTn File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_16S_BLASTn_File.tsv"

      # Empty files to avoid errors - use qcfilecheck function to end workflow
      touch "!{meta.id}-!{meta.assembler}.16S-top-species.tsv"
    fi

    # Add header to BLAST output
    SUMMARY_HEADER=(
      "Query_Name"
      "Reference_Name"
      "Identity_(%)"
      "Alignment_length_(bp)"
      "Mismatches_(#)"
      "Gap_openings_(#)"
      "Query_start_position"
      "Query_end_position"
      "Reference_start_position"
      "Reference_end_position"
      "Expect_value_(e-value)"
      "Bit_score_(bits)"
      "Query_coverage_per_HSP_(%)"
      "Reference_scientific_name"
    )
    SUMMARY_HEADER=$(printf "%s\t" "${SUMMARY_HEADER[@]}" | sed 's/\t$//')
    sed -i "1i ${SUMMARY_HEADER}" "!{blast_output}"

    # Compress input TSV full alignment file to be saved in outdir
    gzip -f !{blast_output}

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
