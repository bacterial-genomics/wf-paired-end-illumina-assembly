process ASSEMBLE_CONTIGS_SPADES {

    label "process_high"
    tag { "${meta.id}" }
    container "staphb/spades@sha256:5df39e8404df2678ccc6c6ed9d7aa0e59b79dfa798aef7fd4fc06cc86ba0b4c0"

    input:
    tuple val(meta), path(cleaned_fastq_files)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Raw_Assembly_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}_contigs.fasta")        , emit: contigs
    path("${meta.id}-${meta.assembler}*.{log,gfa,gz,fasta}")
    path(".command.{out,err}")
    path("versions.yml")                                                       , emit: versions

    shell:
    mode_list  = ["--isolate", "--sc", "--meta", "--plasmid", "--rna", "--metaviral", "--metaplasmid", "--corona"]
    mode       = (params.spades_mode !in mode_list) ? "" : params.spades_mode
    memory     = Math.round(Math.floor(task.memory.toString().replaceAll("[GB]", "").toFloat()))
    '''
    source bash_functions.sh

    # Run SPAdes assembler; try up to 3 times
    msg "INFO: Assembling contigs using SPAdes"

    spades.py \
      -1 "!{meta.id}_R1.paired.fq.gz" \
      -2 "!{meta.id}_R2.paired.fq.gz" \
      -s "!{meta.id}_single.fq.gz" \
      -o SPAdes \
      -k !{params.spades_kmer_sizes} \
      !{mode} \
      --memory !{memory} \
      --threads !{task.cpus}

    # Verify file output
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    if verify_minimum_file_size "SPAdes/contigs.fasta" 'Raw Assembly FastA File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{meta.id}\tRaw Assembly FastA File\tPASS"  \
        >> "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    else
      echo -e "!{meta.id}\tRaw Assembly FastA File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    fi

    if grep -E -q 'N{60}' "SPAdes/contigs.fasta"; then
      # avoid this again: https://github.com/ablab/spades/issues/273
      msg "ERROR: contigs.fasta contains 60+ Ns" >&2
      exit 1
    fi

    # Compress log and paramters files for compact storage
    gzip -9f SPAdes/spades.log \
      SPAdes/params.txt

    # Move and rename files
    mv SPAdes/spades.log.gz "!{meta.id}-!{meta.assembler}.log.gz"
    mv SPAdes/params.txt.gz "!{meta.id}-!{meta.assembler}_params.txt.gz"
    mv SPAdes/contigs.fasta "!{meta.id}-!{meta.assembler}_contigs.fasta"
    mv SPAdes/assembly_graph_with_scaffolds.gfa "!{meta.id}-!{meta.assembler}_graph.gfa"

    # Move extra logfiles if exist
    if [ -f SPAdes/warnings.log ]; then
      mv SPAdes/warnings.log "!{meta.id}-!{meta.assembler}_warnings.log"
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        spades: $(spades.py --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
