process POLISH_ASSEMBLY_BWA_PILON {

    label "process_high"
    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/bwa-samtools-pilon@sha256:209ac13b381188b4a72fe746d3ff93d1765044cbf73c3957e4e2f843886ca57f"

    input:
    tuple val(meta), path(cleaned_fastq_files), path(uncorrected_contigs)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.{Filtered,Polished,Binary,Final}*_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.{paired,single}.bam")                       , emit: bam
    path("${meta.id}-${meta.assembler}.{SNPs,InDels}-corrected.cnt.tsv")
    tuple val(meta), path("${meta.id}-${meta.assembler}.fna")                                       , emit: assembly
    path("${meta.id}.Assembly_FastA.SHA512-checksums.tsv")                                          , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                                                            , emit: versions

    shell:
    polish_corrections = (params.spades_polish_corrections >= 1) ? params.spades_polish_corrections : 3
    '''
    source bash_functions.sh

    msg "INFO: evaluating input filesize of !{uncorrected_contigs} ..."

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    if verify_minimum_file_size "!{uncorrected_contigs}" 'Filtered Assembly File' "!{params.min_filesize_filtered_assembly}"; then
      echo -e "!{meta.id}\tFiltered Assembly File\tPASS"  \
        >> "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    else
      echo -e "!{meta.id}\tFiltered Assembly File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    fi

    msg "INFO: assembly filesize meets or exceeds !{params.min_filesize_filtered_assembly}"

    msg "INFO: Polishing !{meta.id} contigs with cleaned paired end reads ..."

    # Set up files to retain InDel and SNP correction counts (each line is a subsequent polishing round)
    echo -e "Sample_name\tCorrection_round\tInDels_corrected_[#]" \
      > "!{meta.id}-!{meta.assembler}.InDels-corrected.cnt.tsv"
    echo -e "Sample_name\tCorrection_round\tSNPs_corrected_[#]" \
      > "!{meta.id}-!{meta.assembler}.SNPs-corrected.cnt.tsv"

    # Set up QC File checks for BAM and FastA for each round of polishing
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Polished_Assembly_File.tsv"

    for (( i=1; i<=!{polish_corrections}; i++ )); do
      msg "INFO: Performing polishing step ${i} of !{polish_corrections} total rounds ..."

      bwa index !{uncorrected_contigs}

      msg "INFO: Completed bwa index of !{uncorrected_contigs} FastA assembly file"

      msg "INFO: Cleaned paired read mapping (${i} of !{polish_corrections}) of !{meta.id} ..."

      bwa mem \
        -v 2 \
        -x intractg \
        -t !{task.cpus} \
        !{uncorrected_contigs} \
        "!{meta.id}_R1.paired.fq.gz" "!{meta.id}_R2.paired.fq.gz" \
        | \
        samtools sort \
        -l 9 \
        -@ !{task.cpus} \
        -o "!{meta.id}-!{meta.assembler}.paired.bam" \
        --reference !{uncorrected_contigs}

      msg "INFO: Completed paired-end read mapping (${i} of !{polish_corrections}) of !{meta.id} and formed sorted BAM output file"

      if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.paired.bam" 'Binary PE Alignment Map BAM File' "!{params.min_filesize_binary_pe_alignment}"; then
        echo -e "!{meta.id}\tBinary PE Alignment Map BAM File (${i} of !{polish_corrections} rounds)\tPASS" \
          >> "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
      else
        echo -e "!{meta.id}\tBinary PE Alignment Map BAM File (${i} of !{polish_corrections} rounds)\tFAIL" \
          >> "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
      fi

      samtools index "!{meta.id}-!{meta.assembler}.paired.bam"

      msg "INFO: Completed samtools index of paired-end BAM alignment file for !{meta.id}"

      msg "INFO: Pilon correcting SNPs and InDels (${i} of !{polish_corrections}) of !{meta.id}..."

      pilon \
        --genome !{uncorrected_contigs} \
        --frags "!{meta.id}-!{meta.assembler}.paired.bam" \
        --output "!{meta.id}-!{meta.assembler}" \
        --changes \
        --fix snps,indels \
        --mindepth 0.50 \
        --threads !{task.cpus}

      msg "INFO: Completed Pilon correction of SNPs and InDels (${i} of !{polish_corrections}) for !{meta.id}"

      if verify_minimum_file_size "!{uncorrected_contigs}" 'Polished Assembly File' "!{params.min_filesize_polished_assembly}"; then
        echo -e "!{meta.id}\tPolished Assembly File (${i} of !{polish_corrections}) rounds\tPASS" \
          >> "!{meta.id}-!{meta.assembler}.Polished_Assembly_File.tsv"
      else
        echo -e "!{meta.id}\tPolished Assembly File (${i} of !{polish_corrections} rounds)\tFAIL" \
          >> "!{meta.id}-!{meta.assembler}.Polished_Assembly_File.tsv"
      fi

      # Place round number and number of InDels/SNPs corrected into respective files
      echo -e "!{meta.id}\t${i}\t$(grep -c '-' !{meta.id}-!{meta.assembler}.changes)" \
        >> "!{meta.id}-!{meta.assembler}.InDels-corrected.cnt.tsv"
      echo -e "!{meta.id}\t${i}\t$(grep -vc '-' !{meta.id}-!{meta.assembler}.changes)" \
        >> "!{meta.id}-!{meta.assembler}.SNPs-corrected.cnt.tsv"

      rm -f "!{meta.id}-!{meta.assembler}.{changes,uncorrected.fna}"
      rm -f "!{meta.id}-!{meta.assembler}Pilon.bed"
      mv -f "!{meta.id}-!{meta.assembler}.fasta" \
        "!{meta.id}-!{meta.assembler}.uncorrected.fna"

      sed -i 's/_pilon//1' "!{meta.id}-!{meta.assembler}.uncorrected.fna"
    done

    # Final corrected FastA assembly file handling
    mv -f "!{meta.id}-!{meta.assembler}.uncorrected.fna" "!{meta.id}-!{meta.assembler}.fna"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.fna" 'Final Corrected Assembly FastA File' "!{params.min_filesize_final_assembly}"; then
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    else
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    fi

    # Single read mapping if available for downstream depth of coverage
    #  calculations, not for assembly polishing.
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Binary_SE_Alignment_Map_File.tsv"

    if [[ !{meta.id}_single.fq.gz ]]; then
      msg "INFO: Single read mapping of !{meta.id}..."

      bwa index "!{meta.id}-!{meta.assembler}.fna"

      msg "INFO: Completed bwa index of !{meta.id}-!{meta.assembler}.fna FastA assembly file for single-end read mapping"

      bwa mem \
        -v 2 \
        -x intractg \
        "!{meta.id}-!{meta.assembler}.fna" \
        -t !{task.cpus} \
        "!{meta.id}_single.fq.gz" \
        | \
        samtools sort \
        -l 9 \
        -@ !{task.cpus} \
        -o "!{meta.id}-!{meta.assembler}.single.bam" \
        --reference "!{meta.id}-!{meta.assembler}.fna"

      msg "INFO: Completed single-end read mapping of !{meta.id} and formed sorted BAM output file"

      if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.single.bam" 'Binary SE Alignment Map File' '!{params.min_filesize_binary_se_alignment}'; then
        echo -e "!{meta.id}\tBinary SE Alignment Map File\tPASS" \
            >> "!{meta.id}-!{meta.assembler}.Binary_SE_Alignment_Map_File.tsv"
      else
        echo -e "!{meta.id}\tBinary SE Alignment Map File\tFAIL" \
          >> "!{meta.id}-!{meta.assembler}.Binary_SE_Alignment_Map_File.tsv"
      fi

      samtools index "!{meta.id}-!{meta.assembler}.single.bam"

      msg "INFO: Completed samtools index of single-end BAM alignment file for !{meta.id}"

    fi

    # Calculate checksum
    FILE="!{meta.id}-!{meta.assembler}.fna"
    CHECKSUM=$(awk '/^>/ {print substr($1, 1)} !/^>/ {print}' "${FILE}" | sha512sum | awk '{print $1}')
    echo "${CHECKSUM}" | awk -v sample_id="!{meta.id}" -v file="${FILE}" '
        BEGIN {
            # Print the header once
            print "Sample_name\tChecksum_(SHA-512)\tFile"
        }
        {
            # Print the data row once, using the CHECKSUM from input
            print sample_id "\t" $1 "\t" file
        }' \
        > "!{meta.id}.Assembly_FastA.SHA512-checksums.tsv"

    msg "INFO: Calculated checksum of polished FastA assembly file for !{meta.id}"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bwa: $(bwa 2>&1 | head -n 3 | tail -1 | awk 'NF>1{print $NF}')
        pilon: $(pilon --version | cut -d ' ' -f 3)
        samtools: $(samtools --version | head -n 1 | awk 'NF>1{print $NF}')
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
    END_VERSIONS
    '''
}
