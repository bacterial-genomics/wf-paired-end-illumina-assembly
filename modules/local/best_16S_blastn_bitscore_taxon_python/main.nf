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
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    tag { "${prefix}" }

    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(prefix), path(blast_tsv), path(qc_aligned_blast_filecheck), path(assembly)

    output:
    path "${prefix}.Summary.16S.tab", emit: blast_summary
    path "${prefix}.16S-top-species.tsv", emit: ssu_species
    path "${prefix}.Filtered_16S_BLASTn_File.tsv", emit: qc_filtered_blastn_filecheck
    path "${prefix}.blast.tsv.gz"
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

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

    # Get filter.blast.py and check if it exists
    filter_blast_script="${DIR}/filter.blast.py"
    if ! check_if_file_exists_allow_seconds ${filter_blast_script} '60'; then
      exit 1
    fi

    # Get the top match by bitscore
    python ${filter_blast_script} \
      -i "!{blast_tsv}" \
      -o "!{prefix}.blast.tab" \
      -c !{params.filter_blast_column} \
      -s !{params.filter_blast_bitscore}

    if verify_minimum_file_size "!{prefix}.blast.tab" 'Filtered 16S BLASTn File' "!{params.min_filesize_filtered_blastn}"; then
      echo -e "!{prefix}\tFiltered 16S BLASTn File\tPASS" \
        > !{prefix}.Filtered_16S_BLASTn_File.tsv

      # Report the top alignment match data: %nucl iden, %query cov aln, taxon
      awk -F $'\t' 'BEGIN{OFS=FS}; {print $1, $3 "% identity", $13 "% alignment", $14}' \
      "!{prefix}.blast.tab" \
      > "!{prefix}.16S-top-species.tsv"

      cat "!{prefix}.16S-top-species.tsv" >> "!{prefix}.Summary.16S.tab"
      gzip -f !{blast_tsv}

    else
      echo -e "!{prefix}\tFiltered 16S BLASTn File\tFAIL" \
        > !{prefix}.Filtered_16S_BLASTn_File.tsv

      # Empty files to avoid errors
      touch !{prefix}.16S-top-species.tsv !{prefix}.Summary.16S.tab !{blast_tsv}.gz
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      python: $(python --version 2>&1 | awk '{print $2}')
      biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
