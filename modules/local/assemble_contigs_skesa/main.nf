process ASSEMBLE_CONTIGS_SKESA {

    label "process_high"
    tag { "${meta.id}" }
    container "staphb/skesa@sha256:b520da51cd3929683c5eb94739bcd6c32045863dab16e777a4e02d2ff3802f20"

    input:
    tuple val(meta), path(cleaned_fastq_files)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Raw_Assembly_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-SKESA_contigs.fasta")                    , emit: contigs
    path(".command.{out,err}")
    path("versions.yml")                                                       , emit: versions

    shell:
    allow_snps = params.skesa_allow_snps ? "--allow snps" : ""
    memory     = Math.round(Math.floor(task.memory.toString().replaceAll("[GB]", "").toFloat()))
    '''
    source bash_functions.sh

    msg "INFO: Assembling !{meta.id} contigs using SKESA ..."

    if [[ ! -f "!{meta.id}-SKESA_contigs.fasta" ]]; then
      skesa \
        --reads "!{meta.id}_R1.paired.fq.gz","!{meta.id}_R2.paired.fq.gz" \
        --reads "!{meta.id}_single.fq.gz" \
        !{allow_snps} \
        --cores !{task.cpus} \
        --memory !{memory} \
        --contigs_out "!{meta.id}-SKESA_contigs.fasta" \
        --steps !{params.skesa_steps} \
        --kmer !{params.skesa_kmer_length} \
        --fraction !{params.skesa_fraction} \
        --max_snp_len !{params.skesa_max_snp_length} \
        --min_contig !{params.skesa_min_contig_length} \
        --vector_percent !{params.skesa_vector_percent}
    fi

    msg "INFO: Completed genome assembly for !{meta.id} using SKESA"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    if verify_minimum_file_size "!{meta.id}-SKESA_contigs.fasta" 'Raw Assembly FastA File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{meta.id}\tRaw Assembly FastA File\tPASS"  \
        >> "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    else
      echo -e "!{meta.id}\tRaw Assembly FastA File\tFAIL" > "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    fi

    msg "INFO: Completed QC file checks for !{meta.id} SKESA"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      skesa: $(skesa --version 2>&1 | grep 'SKESA' | cut -d ' ' -f 2)
    END_VERSIONS
    '''
}
