process ASSEMBLE_SPADES {

    // errorStrategy 'terminate'

    publishDir "${params.outpath}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "${base}"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.Raw_Assembly_File.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    label "process_high"
    tag { "${base}" }

    container "staphb/spades@sha256:e9c50ffb4b6f0ce4d3c504dd0ce1cb3381ae942ff4d5bac24dc78119b3bfd0dd"

    input:
        tuple val(base), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck)

    output:
        path "${base}/"
        tuple val(base), path("${base}/contigs.fasta"), path("*File*.tsv"), emit: contigs
        path "${base}.Raw_Assembly_File.tsv", emit: qc_raw_assembly_filecheck
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
            msg "FAILURE: ${error_message} Check FAILED" >&2
            exit 1
          else
            rm ${filecheck}
          fi
        done

        # Assemble with SPAdes
        failed=0

        msg "INFO: Running SPAdes with !{task.cpus} threads"

        while [[ ! -f !{base}_tmp/contigs.fasta ]] && [ ${failed} -lt 2 ]; do
          RAMSIZE_TOT=$(echo !{task.memory} | cut -d ' ' -f 1)
          msg "INFO: RAMSIZE = ${RAMSIZE_TOT}"
          if [ ${failed} -gt 0 ]; then
            msg "ERROR: assembly file not produced by SPAdes for !{base}" >&2
            mv -f !{base}_tmp/spades.log \
            !{base}_tmp/"${failed}"of3-asm-attempt-failed.spades.log 2> /dev/null
            msg "INFO: SPAdes failure ${failed}; retrying assembly for !{base}" >&2
            spades.py \
            --restart-from last \
            -o !{base}_tmp \
            -t !{task.cpus} >&2
          else
            spades.py \
            --pe1-1 !{paired_R1_gz} \
            --pe1-2 !{paired_R2_gz} \
            --pe1-s !{single_gz} \
            --memory "${RAMSIZE_TOT}" \
            -o !{base}_tmp \
            --phred-offset 33 \
            -t !{task.cpus} \
            --only-assembler >&2
          fi
          failed=$(( ${failed}+1 ))
        done

        if verify_minimum_file_size "!{base}_tmp/contigs.fasta" 'Raw Assembly File' "!{params.min_filesize_raw_assembly}"; then
          echo -e "!{base}\tRaw Assembly File\tPASS" > !{base}.Raw_Assembly_File.tsv
        else
          echo -e "!{base}\tRaw Assembly File\tFAIL" > !{base}.Raw_Assembly_File.tsv
        fi

        if grep -E -q 'N{60}' "!{base}_tmp/contigs.fasta"; then
          # avoid this again: https://github.com/ablab/spades/issues/273
          msg "ERROR: contigs.fasta contains 60+ Ns" >&2
          exit 1
        fi

        gzip !{base}_tmp/spades.log \
        !{base}_tmp/params.txt

        mkdir !{base}
        mv !{base}_tmp/spades.log.gz \
        !{base}_tmp/params.txt.gz \
        !{base}_tmp/contigs.fasta \
        !{base}_tmp/assembly_graph_with_scaffolds.gfa \
        !{base}/
        
        if [ -f !{base}_tmp/warnings.log ]; then
          mv !{base}_tmp/warnings.log !{base}/
        fi

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            spades: $(spades.py --version 2>&1 | awk 'NF>1{print $NF}')
        END_VERSIONS
        '''
}