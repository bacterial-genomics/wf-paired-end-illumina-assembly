process ALIGN_16S_BLAST {

    tag { "${meta.id}" }
    container "gregorysprenger/ncbi-blast-plus@sha256:f187706adb753c44f50e5be82d85c518e9cd0ae090bc30ce5e14bb35565a380a"

    input:
    tuple val(meta), path(barnapp_extracted_rna), path(assembly)
    val database

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                          , emit: versions
    path "${meta.id}.16S_BLASTn_Output_File.tsv" , emit: qc_filecheck
    tuple val(meta), path("${meta.id}.blast.tsv"), emit: blast_output

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Performing BLASTn alignments"

    blastn \
      -word_size 10 \
      -task blastn \
      -db !{database} \
      -num_threads "!{task.cpus}" \
      -query "!{barnapp_extracted_rna}" \
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
