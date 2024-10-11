process ALIGN_16S_BLAST {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/ncbi-blast-plus@sha256:f187706adb753c44f50e5be82d85c518e9cd0ae090bc30ce5e14bb35565a380a"

    input:
    tuple val(meta), path(barnapp_extracted_rna)
    tuple val(db_name), path("database/*")

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.16S_BLASTn_Output_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.blast.tsv")                 , emit: blast_output
    path(".command.{out,err}")
    path("versions.yml")                                                            , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Performing BLASTn alignments"

    export BLASTDB=database

    blastn \
      -word_size 10 \
      -task blastn \
      -db "!{db_name}" \
      -num_threads "!{task.cpus}" \
      -query "!{barnapp_extracted_rna}" \
      -out "!{meta.id}-!{meta.assembler}.blast.tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ssciname"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.16S_BLASTn_Output_File.tsv"
    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.blast.tsv" '16S BLASTn Output TSV File' "!{params.min_filesize_blastn_output}"; then
      echo -e "!{meta.id}\t16S BLASTn Output TSV File\tPASS"  \
        >> "!{meta.id}-!{meta.assembler}.16S_BLASTn_Output_File.tsv"
    else
      echo -e "!{meta.id}\t16S BLASTn Output TSV File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.16S_BLASTn_Output_File.tsv"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        blast: $(blastn -version | head -n 1 | awk '{print $2}')
    END_VERSIONS
    '''
}
