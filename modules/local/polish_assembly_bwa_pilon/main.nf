process POLISH_ASSEMBLY_BWA_PILON {

    label "process_high"
    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/bwa-samtools-pilon@sha256:209ac13b381188b4a72fe746d3ff93d1765044cbf73c3957e4e2f843886ca57f"

    input:
    tuple val(meta), path(cleaned_fastq_files), path(uncorrected_contigs)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.{Filtered,Polished,Binary,Final}*_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.{paired,single}.bam")                       , emit: bam
    path("${meta.id}-${meta.assembler}.{SNPs,InDels}-corrected.cnt.txt")
    tuple val(meta), path("${meta.id}-${meta.assembler}.fna")                                       , emit: assembly
    path(".command.{out,err}")
    path("versions.yml")                                                                            , emit: versions

    shell:
    polish_corrections = (params.spades_polish_corrections >= 1) ? params.spades_polish_corrections : 3
    '''
    source bash_functions.sh

    # Correct cleaned SPAdes contigs with cleaned PE reads
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    if verify_minimum_file_size "!{uncorrected_contigs}" 'Filtered Assembly File' "!{params.min_filesize_filtered_assembly}"; then
      echo -e "!{meta.id}-!{meta.assembler}\tFiltered Assembly File\tPASS"  \
        >> "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    else
      echo -e "!{meta.id}-!{meta.assembler}\tFiltered Assembly File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    fi

    # Set up files to retain InDel and SNP correction information
    echo -e "Correction round\tNumber of InDels corrected" \
      > "!{meta.id}-!{meta.assembler}.InDels-corrected.cnt.txt"
    echo -e "Correction round\tNumber of SNPs corrected" \
      > "!{meta.id}-!{meta.assembler}.SNPs-corrected.cnt.txt"

    msg "INFO: Polishing contigs with paired end reads.."

    for (( i=1; i<=!{polish_corrections}; i++ )); do
      msg "INFO: Performing polishing step ${i} of !{polish_corrections}"

      bwa index !{uncorrected_contigs}

      bwa mem \
        -v 2 \
        -x intractg \
        -t !{task.cpus} \
        !{uncorrected_contigs} \
        !{cleaned_fastq_files[0]} !{cleaned_fastq_files[1]} \
        | \
        samtools sort \
        -l 9 \
        -@ !{task.cpus} \
        -o "!{meta.id}-!{meta.assembler}.paired.bam" \
        --reference !{uncorrected_contigs}

      echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
      if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.paired.bam" 'Binary PE Alignment Map File' "!{params.min_filesize_binary_pe_alignment}"; then
        echo -e "!{meta.id}-!{meta.assembler}\tBinary PE Alignment Map File (${i} of !{polish_corrections})\tPASS" \
          >> "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
      else
        echo -e "!{meta.id}-!{meta.assembler}\tBinary PE Alignment Map File (${i} of !{polish_corrections})\tFAIL" \
          >> "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
      fi

      samtools index "!{meta.id}-!{meta.assembler}.paired.bam"

      pilon \
        --genome !{uncorrected_contigs} \
        --frags "!{meta.id}-!{meta.assembler}.paired.bam" \
        --output "!{meta.id}-!{meta.assembler}" \
        --changes \
        --fix snps,indels \
        --mindepth 0.50 \
        --threads !{task.cpus} >&2

      echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Polished_Assembly_File.tsv"
      if verify_minimum_file_size "!{uncorrected_contigs}" 'Polished Assembly File' "!{params.min_filesize_polished_assembly}"; then
        echo -e "!{meta.id}-!{meta.assembler}\tPolished Assembly File (${i} of !{polish_corrections})\tPASS" \
          >> "!{meta.id}-!{meta.assembler}.Polished_Assembly_File.tsv"
      else
        echo -e "!{meta.id}-!{meta.assembler}\tPolished Assembly File (${i} of !{polish_corrections})\tFAIL" \
          >> "!{meta.id}-!{meta.assembler}.Polished_Assembly_File.tsv"
      fi

      # Place round number and number of InDels/SNPs corrected into respective files
      echo -e "${i}\t$(grep -c '-' !{meta.id}-!{meta.assembler}.changes)" \
        >> "!{meta.id}-!{meta.assembler}.InDels-corrected.cnt.txt"
      echo -e "${i}\t$(grep -vc '-' !{meta.id}-!{meta.assembler}.changes)" \
        >> "!{meta.id}-!{meta.assembler}.SNPs-corrected.cnt.txt"

      rm -f "!{meta.id}-!{meta.assembler}.{changes,uncorrected.fna}"
      rm -f "!{meta.id}-!{meta.assembler}Pilon.bed"
      mv -f "!{meta.id}-!{meta.assembler}.fasta" \
        "!{meta.id}-!{meta.assembler}.uncorrected.fna"

      sed -i 's/_pilon//1' "!{meta.id}-!{meta.assembler}.uncorrected.fna"
    done

    mv -f "!{meta.id}-!{meta.assembler}.uncorrected.fna" "!{meta.id}-!{meta.assembler}.fna"

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.fna" 'Final Corrected Assembly FastA File' "!{params.min_filesize_final_assembly}"; then
      echo -e "!{meta.id}-!{meta.assembler}\tFinal Corrected Assembly FastA File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    else
      echo -e "!{meta.id}-!{meta.assembler}\tFinal Corrected Assembly FastA File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    fi

    # Single read mapping if available for downstream depth of coverage
    #  calculations, not for assembly polishing.
    if [[ !{cleaned_fastq_files[2]} ]]; then
      msg "INFO: Single read mapping"
      bwa index "!{meta.id}-!{meta.assembler}.fna"

      bwa mem \
        -v 2 \
        -x intractg \
        "!{meta.id}-!{meta.assembler}.fna" \
        -t !{task.cpus} \
        !{cleaned_fastq_files[2]} \
        | \
        samtools sort \
        -l 9 \
        -@ !{task.cpus} \
        -o "!{meta.id}-!{meta.assembler}.single.bam" \
        --reference "!{meta.id}-!{meta.assembler}.fna"

      echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Binary_SE_Alignment_Map_File.tsv"
      if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.single.bam" 'Binary SE Alignment Map File' '!{params.min_filesize_binary_se_alignment}'; then
        echo -e "!{meta.id}-!{meta.assembler}\tBinary SE Alignment Map File\tPASS" \
            >> "!{meta.id}-!{meta.assembler}.Binary_SE_Alignment_Map_File.tsv"
      else
        echo -e "!{meta.id}-!{meta.assembler}\tBinary SE Alignment Map File\tFAIL" \
          >> "!{meta.id}-!{meta.assembler}.Binary_SE_Alignment_Map_File.tsv"
      fi

      samtools index "!{meta.id}-!{meta.assembler}.single.bam"
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
