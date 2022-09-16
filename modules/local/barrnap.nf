process BARRNAP {

    publishDir "${params.outpath}/ssu",
        mode: "${params.publish_dir_mode}",
        pattern: "*.fa"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/barrnap@sha256:e22cbd789c36d5626460feb6c7e5f6f7d55c8628dacae68ba0da30884195a837"

    input:
        path extracted_rna
        path base_fna
        path annotation
        val base

    output:
        path "16s.*.fa", emit: extracted_base
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''

        source bash_functions.sh

        if [[ ! -f "!{extracted_rna}" ]] || [[ ! -s "!{extracted_rna}" ]]; then
            msg "INFO: absent 16S rRNA gene annotation in !{annotation};" >&2
            msg 'Running barrnapp' >&2
            barrnap !{base_fna} > !{base}.gff
            bedtools getfasta \
                -fi !{base_fna} \
                -bed !{base}.gff \
                -fo 16S.!{base}.fa

            if [[ $(grep -c '>' "!{extracted_rna}") -eq 0 ]]; then
                msg "INFO: RNAmmer was unable to locate a 16S rRNA gene sequence in !{base_fna}" >&2
                rm "16S.!{base}.fa"
                exit 2
            fi
        fi

        verify_file_minimum_size "16S.!{base}.fa" '16S extracted FastA file' '500c'

        awk -v awk_var="!{base}" '/^>/{print ">" awk_var "_" ++i; next} {print}' \
        16S.!{base}.fa > !{base}.fa-renamed
        rm -f 16S.!{base}.fa
        mv -f !{base}.fa-renamed 16s.!{base}.fa

        verify_file_minimum_size "16s.!{base}.fa" '16S extracted and renamed FastA file' '500c'

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            barrnap: $(barrnap --version 2>&1 | awk 'NF>1{print $NF}')
        END_VERSIONS

        '''
}