process ASSEMBLE_CONTIGS_SPADES {

    label "process_high"
    tag { "${meta.id}" }
    container "gregorysprenger/spades@sha256:3fe1ebda8f5746ca3e3ff79c74c220d2ca75e3120f20441c3e6ae88eff03b4dc"

    input:
    tuple val(meta), path(cleaned_fastq_files)

    output:
    path("${meta.id}-${meta.assembler}*.{log,gfa,gz,fasta}")
    path(".command.{out,err}")
    path("versions.yml")                                                       , emit: versions
    tuple val(meta), path("${meta.id}-${meta.assembler}.Raw_Assembly_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}_contigs.fasta")        , emit: contigs

    shell:
    mode_list = ["--isolate", "--sc", "--meta", "--plasmid", "--rna", "--metaviral", "--metaplasmid", "--corona"]
    mode = (params.spades_mode !in mode_list) ? "" : params.spades_mode
    '''
    source bash_functions.sh

    # Run SPAdes assembler; try up to 3 times
    msg "INFO: Assembling contigs using SPAdes"
    failed=0
    while [[ ! -f SPAdes/contigs.fasta ]] && [ ${failed} -lt 2 ]; do
      RAMSIZE=$(echo !{task.memory} | cut -d ' ' -f 1)

      if [ ${failed} -gt 0 ]; then
        msg "ERROR: assembly file not produced by SPAdes for !{meta.id}" >&2
        mv -f SPAdes/spades.log \
          SPAdes/"${failed}"of3-asm-attempt-failed.spades.log 2> /dev/null
        msg "INFO: SPAdes failure ${failed}; retrying assembly for !{meta.id}" >&2

        spades.py \
          --restart-from last \
          -o SPAdes \
          -t !{task.cpus} >&2

      else

        spades.py \
          -1 !{cleaned_fastq_files[0]} \
          -2 !{cleaned_fastq_files[1]} \
          -s !{cleaned_fastq_files[2]} \
          -o SPAdes \
          -k !{params.spades_kmer_sizes} \
          !{mode} \
          --memory "${RAMSIZE}" \
          --threads !{task.cpus}
      fi
      failed=$(( ${failed}+1 ))
    done

    # Verify file output
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    if verify_minimum_file_size "SPAdes/contigs.fasta" 'Raw Assembly File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{meta.id}\tRaw Assembly File\tPASS"  \
        >> "!{meta.id}-!{meta.assembler}.Raw_Assembly_File.tsv"
    else
      echo -e "!{meta.id}\tRaw Assembly File\tFAIL" \
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
    # if [ -f SPAdes/*of3-asm-attempt-failed.spades.log ]; then
    #   mv SPAdes/*of3-asm-attempt-failed.spades.log "SPAdes/!{meta.id}"
    # fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        spades: $(spades.py --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
