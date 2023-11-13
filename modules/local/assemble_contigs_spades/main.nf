process ASSEMBLE_CONTIGS_SPADES {

    label "process_high"
    tag { "${meta.id}" }
    container "gregorysprenger/spades@sha256:3fe1ebda8f5746ca3e3ff79c74c220d2ca75e3120f20441c3e6ae88eff03b4dc"

    input:
    tuple val(meta), path(R1), path(R2), path(single)

    output:
    path "${meta.id}/"
    path ".command.out"
    path ".command.err"
    path "versions.yml"                              , emit: versions
    path "${meta.id}.Raw_Assembly_File.tsv"          , emit: qc_filecheck
    tuple val(meta), path("${meta.id}/contigs.fasta"), emit: contigs

    shell:
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
          -1 !{R1} \
          -2 !{R2} \
          -s !{single} \
          -o SPAdes \
          --memory "${RAMSIZE}" \
          --threads !{task.cpus}
          # param.mode
          # -k !{params.kmer_sizes} \

      fi
      failed=$(( ${failed}+1 ))
    done

    # Verify file output
    if verify_minimum_file_size "SPAdes/contigs.fasta" 'Raw Assembly File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{meta.id}\tRaw Assembly File\tPASS" > !{meta.id}.Raw_Assembly_File.tsv
    else
      echo -e "!{meta.id}\tRaw Assembly File\tFAIL" > !{meta.id}.Raw_Assembly_File.tsv
    fi

    if grep -E -q 'N{60}' "SPAdes/contigs.fasta"; then
      # avoid this again: https://github.com/ablab/spades/issues/273
      msg "ERROR: contigs.fasta contains 60+ Ns" >&2
      exit 1
    fi

    # Compress log and paramters files for compact storage
    gzip -9f SPAdes/spades.log \
      SPAdes/params.txt

    # Most a few spades files into a new sample name dir for storage
    mkdir -p !{meta.id}
    mv -t !{meta.id}/ \
      SPAdes/spades.log.gz \
      SPAdes/params.txt.gz \
      SPAdes/contigs.fasta \
      SPAdes/assembly_graph_with_scaffolds.gfa

    # Move extra logfiles if exist
    if [ -f SPAdes/warnings.log ]; then
      mv SPAdes/warnings.log !{meta.id}/
    fi
    if [ -f SPAdes/*of3-asm-attempt-failed.spades.log ]; then
      mv SPAdes/*of3-asm-attempt-failed.spades.log !{meta.id}/
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        spades: $(spades.py --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
