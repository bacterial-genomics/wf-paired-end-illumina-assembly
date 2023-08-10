process MAP_CONTIGS_BWA {

    publishDir "${params.outdir}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "*.fna"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*File*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${meta.id}.${task.process}${filename}"}

    label "process_high"
    tag { "${meta.id}" }

    container "gregorysprenger/bwa-samtools-pilon@sha256:209ac13b381188b4a72fe746d3ff93d1765044cbf73c3957e4e2f843886ca57f"

    input:
    tuple val(meta), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck), path(uncorrected_contigs)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                                                              , emit: versions
    tuple val(meta), path("${meta.id}.fna")                                                          , emit: assembly
    path "${meta.id}.Filtered_Assembly_File.tsv"                                                     , emit: qc_filtered_asm_filecheck
    path "${meta.id}.Binary_PE_Alignment_Map_File.tsv"                                               , emit: qc_pe_alignment_filecheck
    path "${meta.id}.Binary_SE_Alignment_Map_File.tsv"                                               , emit: qc_se_alignment_filecheck
    path "${meta.id}.Final_Corrected_Assembly_FastA_File.tsv"                                        , emit: qc_corrected_asm_filecheck
    tuple val(meta), path("${meta.id}.paired.bam"), path("${meta.id}.single.bam"), path("*File*.tsv"), emit: bam

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

    # Map SKESA contigs with cleaned PE reads
    if verify_minimum_file_size "!{uncorrected_contigs}" 'Filtered Assembly File' "!{params.min_filesize_filtered_assembly}"; then
      echo -e "!{meta.id}\tFiltered Assembly File\tPASS" > !{meta.id}.Filtered_Assembly_File.tsv
    else
      echo -e "!{meta.id}\tFiltered Assembly File\tFAIL" > !{meta.id}.Filtered_Assembly_File.tsv
    fi

    bwa index !{uncorrected_contigs}

    bwa mem \
      -v 2 \
      -x intractg \
      -t !{task.cpus} \
      !{uncorrected_contigs} \
      !{paired_R1_gz} !{paired_R2_gz} \
      | \
      samtools sort \
      -@ !{task.cpus} \
      -l 9 \
      -o !{meta.id}.paired.bam \
      --reference !{uncorrected_contigs}

    if verify_minimum_file_size "!{meta.id}.paired.bam" 'Binary PE Alignment Map File' "!{params.min_filesize_binary_pe_alignment}"; then
      echo -e "!{meta.id}\tBinary PE Alignment Map File\tPASS" \
        >> !{meta.id}.Binary_PE_Alignment_Map_File.tsv
    else
      echo -e "!{meta.id}\tBinary PE Alignment Map File\tFAIL" \
        >> !{meta.id}.Binary_PE_Alignment_Map_File.tsv
    fi

    samtools index !{meta.id}.paired.bam

    mv -f !{meta.id}.uncorrected.fna !{meta.id}.fna

    if verify_minimum_file_size "!{meta.id}.fna" 'Final Corrected Assembly FastA File' "!{params.min_filesize_final_assembly}"; then
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tPASS" \
        > !{meta.id}.Final_Corrected_Assembly_FastA_File.tsv
    else
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tFAIL" \
        > !{meta.id}.Final_Corrected_Assembly_FastA_File.tsv
    fi

    # Single read mapping if available for downstream depth of coverage calculations
    if [[ !{single_gz} ]]; then
      msg "INFO: Single read mapping with !{task.cpus} threads"
      bwa index !{meta.id}.fna

      bwa mem \
        -v 2 \
        -x intractg \
        !{single_gz} \
        !{meta.id}.fna \
        -t !{task.cpus} \
        | \
        samtools sort \
        -l 9 \
        -@ !{task.cpus} \
        -o !{meta.id}.single.bam \
        --reference !{meta.id}.fna

      if verify_minimum_file_size "!{meta.id}.single.bam" 'Binary SE Alignment Map File' '!{params.min_filesize_binary_se_alignment}'; then
        echo -e "!{meta.id}\tBinary SE Alignment Map File\tPASS" \
          > !{meta.id}.Binary_SE_Alignment_Map_File.tsv
      else
        echo -e "!{meta.id}\tBinary SE Alignment Map File\tFAIL" \
          > !{meta.id}.Binary_SE_Alignment_Map_File.tsv
      fi

      samtools index !{meta.id}.single.bam
    fi

    # Get process version
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      bwa: $(bwa 2>&1 | head -n 3 | tail -1 | awk 'NF>1{print $NF}')
      samtools: $(samtools --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
