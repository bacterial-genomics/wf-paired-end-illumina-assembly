process ASSEMBLE_CONTIGS_SKESA {

    label "process_high"
    tag { "${meta.id}" }
    container "gregorysprenger/skesa@sha256:4455882b5d0fd968630325428729395422be7340301c31d15874a295904b7f26"

    input:
    tuple val(meta), path(cleaned_fastq_files)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Raw_Assembly_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("contigs.fasta")                                     , emit: contigs
    path(".command.{out,err}")
    path("versions.yml")                                                       , emit: versions

    shell:
    allow_snps = params.skesa_allow_snps ? "--allow snps" : ""
    '''
    source bash_functions.sh

    msg "INFO: Assembling contigs using SKESA"

    if [[ ! -f contigs.fasta ]]; then
      skesa \
        --reads !{cleaned_fastq_files[0]},!{cleaned_fastq_files[1]} \
        --reads !{cleaned_fastq_files[2]} \
        !{allow_snps} \
        --cores !{task.cpus} \
        --memory !{task.memory} \
        --contigs_out contigs.fasta \
        --steps !{params.skesa_steps} \
        --kmer !{params.skesa_kmer_length} \
        --fraction !{params.skesa_fraction} \
        --max_snp_len !{params.skesa_max_snp_length} \
        --min_contig !{params.skesa_min_contig_length} \
        --vector_percent !{params.skesa_vector_percent}
    fi

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    if verify_minimum_file_size "contigs.fasta" 'Raw Assembly File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{meta.id}\tRaw Assembly File\tPASS"  \
        >> "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    else
      echo -e "!{meta.id}\tRaw Assembly File\tFAIL" > "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      skesa: $(skesa --version 2>&1 | grep 'SKESA' | cut -d ' ' -f 2)
    END_VERSIONS
    '''
}
