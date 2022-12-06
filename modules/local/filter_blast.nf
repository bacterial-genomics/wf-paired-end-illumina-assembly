process FILTER_BLAST {

    publishDir "${params.outpath}/ssu",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv.gz"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    tag { "${base}" }

    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
        tuple val(base), val(size), path(blast_tsv), path(base_fna)

    output:
        path "${base}.Summary.16S.tab", emit: blast_summary
        path "${base}.blast.tsv.gz"
        path "${base}.16S-top-species.tsv", emit: ssu_species
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        source bash_functions.sh

        # Get filter.blast.py and check if it exists
        filter_blast="${DIR}/filter.blast.py"
        check_if_file_exists_allow_seconds ${filter_blast} '60'

        python ${filter_blast} -i "!{blast_tsv}" \
        -o "!{base}.blast.tab"

        verify_file_minimum_size "!{base}.blast.tab" 'filtered 16S blastn nr file' '10c' '100'

        awk -F $'\t' 'BEGIN{OFS=FS}; {print $1, $3 "% identity", $13 "% alignment", $14}' "!{base}.blast.tab" \
        > "!{base}.16S-top-species.tsv"

        cat "!{base}.16S-top-species.tsv" >> "!{base}.Summary.16S.tab"
        gzip -f !{blast_tsv}

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            python: $(python --version 2>&1 | awk '{print $2}')
            biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
        END_VERSIONS
        '''
}