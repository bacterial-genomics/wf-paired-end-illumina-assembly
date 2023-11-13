process ALIGN_16S_BLAST {

    tag { "${meta.id}" }
    container "gregorysprenger/ncbi-blast-plus@sha256:2d3e226d2eb31e3e0d5a80d7325b3a2ffd873ad1f2bd81215fd0b43727019279"

    input:
    tuple val(meta), path(extracted_base), path(assembly)
    val database

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                          , emit: versions
    path "${meta.id}.16S_BLASTn_Output_File.tsv" , emit: qc_filecheck
    tuple val(meta), path("${meta.id}.blast.tsv"), emit: blast_tsv

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Performing BLASTn alignments"

    blastn \
      -word_size 10 \
      -task blastn \
      -db !{database} \
      -num_threads "!{task.cpus}" \
      -query "!{extracted_base}" \
      -out "!{meta.id}.blast.tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ssciname"

    if verify_minimum_file_size "!{meta.id}.blast.tsv" '16S BLASTn Output File' "!{params.min_filesize_blastn_output}"; then
      echo -e "!{meta.id}\t16S BLASTn Output File\tPASS" > !{meta.id}.16S_BLASTn_Output_File.tsv
    else
      echo -e "!{meta.id}\t16S BLASTn Output File\tFAIL" > !{meta.id}.16S_BLASTn_Output_File.tsv
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        blast: $(blastn -version | head -n 1 | awk '{print $2}')
    END_VERSIONS
    '''
}
