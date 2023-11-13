process MAP_CONTIGS_BWA {

    label "process_high"
    tag { "${meta.id}" }
    container "gregorysprenger/bwa-samtools-pilon@sha256:209ac13b381188b4a72fe746d3ff93d1765044cbf73c3957e4e2f843886ca57f"

    input:
    tuple val(meta), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(uncorrected_contigs)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                                          , emit: versions
    tuple val(meta), path("${meta.id}.fna")                                      , emit: assembly
    path("${meta.id}.{Filtered,Binary,Final}*_File.tsv")                         , emit: qc_filecheck
    tuple val(meta), path("${meta.id}.paired.bam"), path("${meta.id}.single.bam"), emit: bam

    shell:
    '''
    source bash_functions.sh

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

    cp -L !{meta.id}.uncorrected.fna !{meta.id}.fna

    if verify_minimum_file_size "!{meta.id}.fna" 'Final Corrected Assembly FastA File' "!{params.min_filesize_final_assembly}"; then
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tPASS" \
        > !{meta.id}.Final_Corrected_Assembly_FastA_File.tsv
    else
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tFAIL" \
        > !{meta.id}.Final_Corrected_Assembly_FastA_File.tsv
    fi

    # Single read mapping if available for downstream depth of coverage calculations
    if [[ !{single_gz} ]]; then
      msg "INFO: Single read mapping"
      bwa index !{meta.id}.fna

      bwa mem \
        -v 2 \
        -x intractg \
        !{meta.id}.fna \
        !{single_gz} \
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

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bwa: $(bwa 2>&1 | head -n 3 | tail -1 | awk 'NF>1{print $NF}')
        samtools: $(samtools --version | head -n 1 | awk 'NF>1{print $NF}')
    END_VERSIONS
    '''
}
