process EXTRACT_RECORDS {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    container "snads/biopython@sha256:bb041f55fd45d0fb577656e2d1f1a9f477d3ba80878b3b42218adff3322ae06e"

    input:
        path annotation
        val base

    output:
        path "16S.*.fa", emit: extracted_rna
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''

        source bash_functions.sh

        # Get extract.record.from.genbank.py and check if it exists
        extract_record="${DIR}/extract.record.from.genbank.py"
        check_if_file_exists_allow_seconds ${extract_record} '60'

        # 16S extraction
        if [[ -s "!{annotation}" ]]; then
            python3 ${extract_record} -i "!{annotation}" \
            -u product -o "16S.!{base}.fa" -q '16S ribosomal RNA' \
            --search-type any_q_is_rec -f rRNA
        fi

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            biopython: $(grep "version" /usr/local/lib/python*/dist-packages/Bio/__init__.py | head -n 1 | awk 'NF>1{print $NF}' | tr -d '"')
        END_VERSIONS

        '''
}