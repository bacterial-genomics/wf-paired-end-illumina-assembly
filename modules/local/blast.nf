process BLAST {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    tag { "${base}" }

    container "gregorysprenger/ncbi-blast-plus@sha256:2d3e226d2eb31e3e0d5a80d7325b3a2ffd873ad1f2bd81215fd0b43727019279"

    input:
        tuple val(base), val(size), path(extracted_base), path(base_fna)
        path blast_db

    output:
        tuple val(base), val(size), path("${base}.blast.tsv"), emit: blast_tsv
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        '''
        source bash_functions.sh

        # Classify each 16S sequence record
        if [[ -d "!{blast_db}" ]]; then
            database="!{blast_db}"
            msg "INFO: Using user specified BLAST database: !{params.blast_db}"
        else
            database="/db"
            msg "INFO: Using pre-loaded 16S rRNA database for BLAST"
        fi

        # Set BLAST database as an environmental variable
        export BLASTDB=${database}

        msg "INFO: Running blastn with !{task.cpus} threads"

        blastn -word_size 10 -task blastn -db 16S_ribosomal_RNA \
        -num_threads "!{task.cpus}" \
        -query "!{extracted_base}" \
        -out "!{base}.blast.tsv" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ssciname"

        verify_file_minimum_size "!{base}.blast.tsv" '16S blastn nr output file' '10c'

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            blast: $(blastn -version | head -n 1 | awk '{print $2}')
        END_VERSIONS
        '''
}