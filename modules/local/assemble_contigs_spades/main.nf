process ASSEMBLE_CONTIGS_SPADES {

    label "process_high"
    tag { "${meta.id}" }
    container "gregorysprenger/spades@sha256:3fe1ebda8f5746ca3e3ff79c74c220d2ca75e3120f20441c3e6ae88eff03b4dc"

    input:
    tuple val(meta), path(R1), path(R2), path(single), path(qc_nonoverlap_filecheck)

    output:
    path "${meta.id}/"
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                                  , emit: versions
    path "${meta.id}.Raw_Assembly_File.tsv"                              , emit: qc_raw_assembly_filecheck
    tuple val(meta), path("${meta.id}/contigs.fasta"), path("*File*.tsv"), emit: contigs

    shell:
    '''
    source bash_functions.sh

    # Exit if previous process fails qc filecheck
    for filecheck in !{qc_nonoverlap_filecheck}; do
      if [[ $(grep "FAIL" ${filecheck}) ]]; then
        error_message=$(awk -F '\t' 'END {print $2}' ${filecheck} | sed 's/[(].*[)] //g')
        msg "${error_message} Check failed" >&2
        exit 1
      else
        rm ${filecheck}
      fi
    done

    # Assemble with SPAdes
    failed=0

    msg "INFO: Assembling contigs using SPAdes"

    while [[ ! -f !{meta.id}_tmp/contigs.fasta ]] && [ ${failed} -lt 2 ]; do
      RAMSIZE=$(echo !{task.memory} | cut -d ' ' -f 1)

      if [ ${failed} -gt 0 ]; then
        msg "ERROR: assembly file not produced by SPAdes for !{meta.id}" >&2

        mv -f !{meta.id}_tmp/spades.log \
        !{meta.id}_tmp/"${failed}"of3-asm-attempt-failed.spades.log 2> /dev/null

        msg "INFO: SPAdes failure ${failed}; retrying assembly for !{meta.id}" >&2

        spades.py \
          --restart-from last \
          -o !{meta.id}_tmp \
          -t !{task.cpus} >&2

      else

        spades.py \
          --pe1-1 !{R1} \
          --pe1-2 !{R2} \
          --pe1-s !{single} \
          -t !{task.cpus} \
          -o !{meta.id}_tmp \
          --phred-offset 33 \
          --memory "${RAMSIZE}" \
          --only-assembler >&2

      fi
      failed=$(( ${failed}+1 ))
    done

    if verify_minimum_file_size "!{meta.id}_tmp/contigs.fasta" 'Raw Assembly File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{meta.id}\tRaw Assembly File\tPASS" > !{meta.id}.Raw_Assembly_File.tsv
    else
      echo -e "!{meta.id}\tRaw Assembly File\tFAIL" > !{meta.id}.Raw_Assembly_File.tsv
    fi

    if grep -E -q 'N{60}' "!{meta.id}_tmp/contigs.fasta"; then
      # avoid this again: https://github.com/ablab/spades/issues/273
      msg "ERROR: contigs.fasta contains 60+ Ns" >&2
      exit 1
    fi

    gzip !{meta.id}_tmp/spades.log \
    !{meta.id}_tmp/params.txt

    mkdir !{meta.id}
    mv !{meta.id}_tmp/spades.log.gz \
      !{meta.id}_tmp/params.txt.gz \
      !{meta.id}_tmp/contigs.fasta \
      !{meta.id}_tmp/assembly_graph_with_scaffolds.gfa \
      !{meta.id}/

    if [ -f !{meta.id}_tmp/warnings.log ]; then
      mv !{meta.id}_tmp/warnings.log !{meta.id}/
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        spades: $(spades.py --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
