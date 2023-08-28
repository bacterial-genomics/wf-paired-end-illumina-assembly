process POLISH_ASSEMBLY_BWA_PILON {

    publishDir "${params.outdir}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{txt,fna}"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*File*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${meta.id}.${task.process}${filename}" }

    label "process_high"
    tag { "${meta.id}" }

    container "gregorysprenger/bwa-samtools-pilon@sha256:209ac13b381188b4a72fe746d3ff93d1765044cbf73c3957e4e2f843886ca57f"

    input:
    tuple val(meta), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck), path(uncorrected_contigs)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                                                              , emit: versions
    path "${meta.id}.SNPs-corrected.cnt.txt"
    path "${meta.id}.InDels-corrected.cnt.txt"
    tuple val(meta), path("${meta.id}.fna")                                                          , emit: assembly
    path "${meta.id}.Filtered_Assembly_File.tsv"                                                     , emit: qc_filtered_asm_filecheck
    path "${meta.id}.Polished_Assembly_File.tsv"                                                     , emit: qc_polished_asm_filecheck
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

    # Correct cleaned SPAdes contigs with cleaned PE reads
    if verify_minimum_file_size "!{uncorrected_contigs}" 'Filtered Assembly File' "!{params.min_filesize_filtered_assembly}"; then
      echo -e "!{meta.id}\tFiltered Assembly File\tPASS" > !{meta.id}.Filtered_Assembly_File.tsv
    else
      echo -e "!{meta.id}\tFiltered Assembly File\tFAIL" > !{meta.id}.Filtered_Assembly_File.tsv
    fi

    echo -n '' > !{meta.id}.InDels-corrected.cnt.txt
    echo -n '' > !{meta.id}.SNPs-corrected.cnt.txt

    # Make sure polish_corrections is at least 1, if not set to 1
    polish_corrections=!{params.spades_polish_corrections}
    if [[ ${polish_corrections} -lt 1 ]]; then
      polish_corrections=1
    fi

    msg "INFO: Polishing contigs with paired end reads.."

    for ((i=1;i<=polish_corrections;i++)); do
      msg "INFO: Performing polishing step ${i} of !{params.spades_polish_corrections}"

      bwa index !{uncorrected_contigs}

      bwa mem \
        -v 2 \
        -x intractg \
        -t !{task.cpus} \
        !{uncorrected_contigs} \
        !{paired_R1_gz} !{paired_R2_gz} \
        | \
        samtools sort \
        -l 9 \
        -@ !{task.cpus} \
        -o !{meta.id}.paired.bam \
        --reference !{uncorrected_contigs}

      if verify_minimum_file_size "!{meta.id}.paired.bam" 'Binary PE Alignment Map File' "!{params.min_filesize_binary_pe_alignment}"; then
        echo -e "!{meta.id}\tBinary PE Alignment Map File (${i} of 3)\tPASS" \
          >> !{meta.id}.Binary_PE_Alignment_Map_File.tsv
      else
        echo -e "!{meta.id}\tBinary PE Alignment Map File (${i} of 3)\tFAIL" \
          >> !{meta.id}.Binary_PE_Alignment_Map_File.tsv
      fi

      samtools index !{meta.id}.paired.bam

      pilon \
        --genome !{uncorrected_contigs} \
        --frags !{meta.id}.paired.bam \
        --output "!{meta.id}" \
        --changes \
        --fix snps,indels \
        --mindepth 0.50 \
        --threads !{task.cpus} >&2

      if verify_minimum_file_size "!{uncorrected_contigs}" 'Polished Assembly File' "!{params.min_filesize_polished_assembly}"; then
        echo -e "!{meta.id}\tPolished Assembly File (${i} of 3)\tPASS" \
          >> !{meta.id}.Polished_Assembly_File.tsv
      else
        echo -e "!{meta.id}\tPolished Assembly File (${i} of 3)\tFAIL" \
          >> !{meta.id}.Polished_Assembly_File.tsv
      fi

      echo $(grep -c '-' !{meta.id}.changes >> !{meta.id}.InDels-corrected.cnt.txt)
      echo $(grep -vc '-' !{meta.id}.changes >> !{meta.id}.SNPs-corrected.cnt.txt)

      rm -f !{meta.id}.{changes,uncorrected.fna}
      rm -f "!{meta.id}"Pilon.bed
      mv -f !{meta.id}.fasta !{meta.id}.uncorrected.fna

      sed -i 's/_pilon//1' !{meta.id}.uncorrected.fna
    done

    mv -f !{meta.id}.uncorrected.fna !{meta.id}.fna

    if verify_minimum_file_size "!{meta.id}.fna" 'Final Corrected Assembly FastA File' "!{params.min_filesize_final_assembly}"; then
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tPASS" \
        > !{meta.id}.Final_Corrected_Assembly_FastA_File.tsv
    else
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tFAIL" \
        > !{meta.id}.Final_Corrected_Assembly_FastA_File.tsv
    fi

    # Single read mapping if available for downstream depth of coverage
    #  calculations, not for assembly polishing.
    if [[ !{single_gz} ]]; then
      msg "INFO: Single read mapping"
      bwa index !{meta.id}.fna

      bwa mem \
        -v 2 \
        -x intractg \
        !{meta.id}.fna \
        -t !{task.cpus} \
        !{single_gz} \
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

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        pilon: $(pilon --version | cut -d ' ' -f 3)
        bwa: $(bwa 2>&1 | head -n 3 | tail -1 | awk 'NF>1{print $NF}')
        samtools: $(samtools --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
