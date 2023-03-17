process EXTRACT_16S_BARRNAP {

    // errorStrategy 'terminate'

    publishDir "${params.outpath}/ssu",
        mode: "${params.publish_dir_mode}",
        pattern: "*.fa"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{SSU_Extracted_File,SSU_Renamed_File}.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    tag { "${base}" }
    
    container "snads/barrnap@sha256:e22cbd789c36d5626460feb6c7e5f6f7d55c8628dacae68ba0da30884195a837"

    input:
        tuple val(base), path(annotation), path(qc_annotated_filecheck), path(base_fna), path(extracted_rna)

    output:
        tuple val(base), path("16S.${base}.fa"), path("*File*.tsv"), emit: extracted_base
        path "${base}.SSU_Extracted_File.tsv", emit: qc_ssu_extracted_filecheck
        path "${base}.SSU_Renamed_File.tsv", emit: qc_ssu_renamed_filecheck
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        source bash_functions.sh

        # Exit if previous process fails qc filecheck
        for filecheck in !{qc_annotated_filecheck}; do
          if [[ $(grep "FAIL" ${filecheck}) ]]; then
            error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
            msg "${error_message} Check FAILED" >&2
            exit 1
          else
            rm ${filecheck}
          fi
        done

        if [[ ! -f "!{extracted_rna}" ]] || [[ ! -s "!{extracted_rna}" ]]; then
          msg "INFO: absent 16S rRNA gene annotation in !{annotation}" >&2
          msg 'Running barrnap' >&2
          
          barrnap !{base_fna} | grep "Name=16S_rRNA;product=16S" > !{base}.gff
          
          if [[ $(grep -c "Name=16S_rRNA;product=16S" "!{base}.gff") -eq 0 ]]; then
            msg "INFO: barrnap was unable to locate a 16S rRNA gene sequence in !{base_fna}" >&2
            exit 2
          fi
          
          bedtools getfasta \
           -fi !{base_fna} \
           -bed !{base}.gff \
           -fo 16S.!{base}.fa
        fi

        if verify_minimum_file_size "16S.!{base}.fa" 'SSU Extracted File' "!{params.min_filesize_ssu_file}"; then
          echo -e "!{base}\tSSU Extracted File\tPASS" > !{base}.SSU_Extracted_File.tsv
        else
          echo -e "!{base}\tSSU Extracted File\tFAIL" > !{base}.SSU_Extracted_File.tsv
        fi

        awk -v awk_var="!{base}" \
         '/^>/{print ">" awk_var "_" ++i; next} {print}' \
         16S.!{base}.fa \
         > !{base}.fa-renamed
        rm -f 16S.!{base}.fa
        mv -f !{base}.fa-renamed 16S.!{base}.fa

        if verify_minimum_file_size "16S.!{base}.fa" 'SSU Renamed File' "!{params.min_filesize_ssu_file}"; then
          echo -e "!{base}\tSSU Renamed File\tPASS" > !{base}.SSU_Renamed_File.tsv
        else
          echo -e "!{base}\tSSU Renamed File\tFAIL" > !{base}.SSU_Renamed_File.tsv
        fi

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            barrnap: $(barrnap --version 2>&1 | awk 'NF>1{print $NF}')
        END_VERSIONS
        '''
}