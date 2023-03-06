process BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON {

    // errorStrategy 'terminate'

    publishDir "${params.outpath}/ssu",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv.gz"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Filtered_16S_BLASTn_File.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    tag { "${base}" }

    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
        tuple val(base), path(blast_tsv), path(base_fna)

    output:
        path "${base}.Summary.16S.tab", emit: blast_summary
        path "${base}.16S-top-species.tsv", emit: ssu_species
        path "${base}.Filtered_16S_BLASTn_File.tsv", emit: qc_filecheck
        path "${base}.blast.tsv.gz"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        source bash_functions.sh

        # Get filter.blast.py and check if it exists
        filter_blast_script="${DIR}/filter.blast.py"
        if ! check_if_file_exists_allow_seconds ${filter_blast_script} '60'; then
          exit 1
        fi

        # Get the top match by bitscore
        python ${filter_blast_script} \
         -i "!{blast_tsv}" \
         -o "!{base}.blast.tab" \
         -c !{params.filter_blast_column} \
         -s !{params.filter_blast_bitscore}

        if verify_minimum_file_size "!{base}.blast.tab" 'Filtered 16S BLASTn File' "!{params.min_filesize_filtered_blastn}"; then
          echo -e "!{base}\tFiltered 16S BLASTn File\tPASS" \
           > !{base}.Filtered_16S_BLASTn_File.tsv
        else
          echo -e "!{base}\tFiltered 16S BLASTn File\tFAIL" \
           > !{base}.Filtered_16S_BLASTn_File.tsv
          exit 1
        fi

        # Report the top alignment match data: %nucl iden, %query cov aln, taxon
        awk -F $'\t' 'BEGIN{OFS=FS}; {print $1, $3 "% identity", $13 "% alignment", $14}' \
         "!{base}.blast.tab" \
         > "!{base}.16S-top-species.tsv"

        cat "!{base}.16S-top-species.tsv" >> "!{base}.Summary.16S.tab"
        gzip -f !{blast_tsv}

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            python: $(python --version 2>&1 | awk '{print $2}')
            biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
        END_VERSIONS
        '''
}