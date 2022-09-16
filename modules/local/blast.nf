process BLAST {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/ncbi-blast-plus@sha256:9200ea627a96b6640e6fdd1b128c08d44b92b34e51e586d5bbf817cfaf540d10"

    input:
        path extracted_base
        path base_fna
        val base

    output:
        path "*.blast.tsv", emit: blast_tsv
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
    '''

    source bash_functions.sh

    # Classify each 16S sequence record
    if [ !{params.blast_db} == null ]; then
        mkdir db
        cd db
        wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz -O db.tar.gz
        tar -xvf db.tar.gz
        rm db.tar.gz
        cd ..
        database=db
        export BLASTDB=${database}
    else
        database=!{params.blast_db}
        export BLASTDB=${database}
    fi
    echo "BLAST DB = ${database}"
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
        blast: $(blast -version | head -n 1 | awk 'NF>1{print $NF}' | cut -d '+' -f 1)
    END_VERSIONS

    '''
}