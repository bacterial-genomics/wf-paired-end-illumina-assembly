process ASSEMBLE_SPADES {

    publishDir "${params.outdir}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "${prefix}"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Raw_Assembly_File.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${prefix}.${task.process}${filename}"}

    label "process_high"
    tag { "${prefix}" }

    container "gregorysprenger/spades@sha256:3fe1ebda8f5746ca3e3ff79c74c220d2ca75e3120f20441c3e6ae88eff03b4dc"

    input:
    tuple val(prefix), path(R1), path(R2), path(single), path(qc_nonoverlap_filecheck)

    output:
    path "${prefix}/"
    tuple val(prefix), path("${prefix}/contigs.fasta"), path("*File*.tsv"), emit: contigs
    path "${prefix}.Raw_Assembly_File.tsv", emit: qc_raw_assembly_filecheck
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

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

    msg "INFO: Running SPAdes with !{task.cpus} threads"

    while [[ ! -f !{prefix}_tmp/contigs.fasta ]] && [ ${failed} -lt 2 ]; do
      RAMSIZE_TOT=$(echo !{task.memory} | cut -d ' ' -f 1)
      msg "INFO: RAMSIZE = ${RAMSIZE_TOT}"
      if [ ${failed} -gt 0 ]; then
        msg "ERROR: assembly file not produced by SPAdes for !{prefix}" >&2
        mv -f !{prefix}_tmp/spades.log \
        !{prefix}_tmp/"${failed}"of3-asm-attempt-failed.spades.log 2> /dev/null
        msg "INFO: SPAdes failure ${failed}; retrying assembly for !{prefix}" >&2
        spades.py \
        --restart-from last \
        -o !{prefix}_tmp \
        -t !{task.cpus} >&2
      else
        spades.py \
        --pe1-1 !{R1} \
        --pe1-2 !{R2} \
        --pe1-s !{single} \
        --memory "${RAMSIZE_TOT}" \
        -o !{prefix}_tmp \
        --phred-offset 33 \
        -t !{task.cpus} \
        --only-assembler >&2
      fi
      failed=$(( ${failed}+1 ))
    done

    if verify_minimum_file_size "!{prefix}_tmp/contigs.fasta" 'Raw Assembly File' "!{params.min_filesize_raw_assembly}"; then
      echo -e "!{prefix}\tRaw Assembly File\tPASS" > !{prefix}.Raw_Assembly_File.tsv
    else
      echo -e "!{prefix}\tRaw Assembly File\tFAIL" > !{prefix}.Raw_Assembly_File.tsv
    fi

    if grep -E -q 'N{60}' "!{prefix}_tmp/contigs.fasta"; then
      # avoid this again: https://github.com/ablab/spades/issues/273
      msg "ERROR: contigs.fasta contains 60+ Ns" >&2
      exit 1
    fi

    gzip !{prefix}_tmp/spades.log \
    !{prefix}_tmp/params.txt

    mkdir !{prefix}
    mv !{prefix}_tmp/spades.log.gz \
    !{prefix}_tmp/params.txt.gz \
    !{prefix}_tmp/contigs.fasta \
    !{prefix}_tmp/assembly_graph_with_scaffolds.gfa \
    !{prefix}/
    
    if [ -f !{prefix}_tmp/warnings.log ]; then
      mv !{prefix}_tmp/warnings.log !{prefix}/
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      spades: $(spades.py --version 2>&1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}