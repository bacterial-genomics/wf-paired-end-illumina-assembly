process FILTER_BLAST {

    publishDir "${params.outpath}/ssu",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv*"
    publishDir "${params.outpath}/qa",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tab"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/biopython@sha256:bb041f55fd45d0fb577656e2d1f1a9f477d3ba80878b3b42218adff3322ae06e"

    input:
        path filter_blast
        path blast_tsv
        path base_fna
        path outpath
        val base

    output:
        path "Summary.16S.tab"
        path "*.blast.tsv.gz"
        path "*species.tsv"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''

        source bash_functions.sh

        python3 !{filter_blast} -i "!{blast_tsv}" \
        -o "!{base}.blast.tab"

        verify_file_minimum_size "!{base}.blast.tab" 'filtered 16S blastn nr file' '10c'

        awk -F $'\t' 'BEGIN{OFS=FS}; {print $1, $3 "% identity", $13 "% alignment", $14}' "!{base}.blast.tab" \
        > "!{base}.16S-top-species.tsv"

        cat "!{base}.16S-top-species.tsv" >> "Summary.16S.tab"
        gzip -f !{blast_tsv}

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            biopython: $(grep "version" /usr/local/lib/python*/dist-packages/Bio/__init__.py | head -n 1 | awk 'NF>1{print $NF}' | tr -d '"')
        END_VERSIONS

        '''
}