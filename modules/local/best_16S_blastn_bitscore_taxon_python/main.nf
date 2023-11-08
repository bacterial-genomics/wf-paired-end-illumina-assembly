process BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON {

    publishDir "${params.outdir}/ssu",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv.gz"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Filtered_16S_BLASTn_File.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${meta.id}.${task.process}${filename}" }

    tag { "${meta.id}" }

    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(blast_tsv), path(qc_aligned_blast_filecheck), path(assembly)

    output:
    path ".command.out"
    path ".command.err"
    path "${meta.id}.blast.tsv.gz"
    path "versions.yml"                           , emit: versions
    path "${meta.id}.Summary.16S.tab"             , emit: blast_summary
    path "${meta.id}.16S-top-species.tsv"         , emit: ssu_species
    path "${meta.id}.Filtered_16S_BLASTn_File.tsv", emit: qc_filtered_blastn_filecheck

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_aligned_blast_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Get the top match by bitscore
    filter.blast.py \
      -i "!{blast_tsv}" \
      -o "!{meta.id}.blast.tab" \
      -c !{params.filter_blast_column} \
      -s !{params.filter_blast_bitscore}

    if verify_minimum_file_size "!{meta.id}.blast.tab" 'Filtered 16S BLASTn File' "!{params.min_filesize_filtered_blastn}"; then
      echo -e "!{meta.id}\tFiltered 16S BLASTn File\tPASS" \
        > !{meta.id}.Filtered_16S_BLASTn_File.tsv

      # Report the top alignment match data: %nucl iden, %query cov aln, taxon
      awk -F $'\t' 'BEGIN{OFS=FS}; {print $1, $3 "% identity", $13 "% alignment", $14}' \
      "!{meta.id}.blast.tab" \
      > "!{meta.id}.16S-top-species.tsv"

      cat "!{meta.id}.16S-top-species.tsv" >> "!{meta.id}.Summary.16S.tab"
      gzip -f !{blast_tsv}

    else
      echo -e "!{meta.id}\tFiltered 16S BLASTn File\tFAIL" \
        > !{meta.id}.Filtered_16S_BLASTn_File.tsv

      # Empty files to avoid errors
      touch !{meta.id}.16S-top-species.tsv !{meta.id}.Summary.16S.tab !{blast_tsv}.gz
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
