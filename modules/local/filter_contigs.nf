process FILTER_CONTIGS {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
        path contigs
        path R1_paired_gz
        val base

    output:
        path "*.uncorrected.fna", emit: uncorrected_contigs
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''

        source bash_functions.sh

        # Get filter.contigs.py and check if it exists
        filter_contigs="${DIR}/filter.contigs.py"
        check_if_file_exists_allow_seconds ${filter_contigs} '60'

        python ${filter_contigs} \
        -i !{contigs}\
        -b "!{base}" -l 1000 -o !{base}.uncorrected.fna

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            biopython: $(grep "version" /usr/local/lib/python*/dist-packages/Bio/__init__.py | head -n 1 | awk 'NF>1{print $NF}' | tr -d '"')
        END_VERSIONS

        '''
}