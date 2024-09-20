process MAP_CONTIGS_BWA {

    label "process_high"
    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/bwa-samtools-pilon@sha256:209ac13b381188b4a72fe746d3ff93d1765044cbf73c3957e4e2f843886ca57f"

    input:
    tuple val(meta), path(cleaned_fastq_files), path(uncorrected_contigs)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.{Filtered,Binary,Final}*_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.{paired,single}.bam")              , emit: bam
    tuple val(meta), path("${meta.id}-${meta.assembler}.fna")                              , emit: assembly
    path("${meta.id}.Assembly_FastA.SHA512-checksums.tsv")                                 , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                                                   , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Map SKESA contigs with cleaned PE reads
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    if verify_minimum_file_size "!{uncorrected_contigs}" 'Filtered Assembly FastA File' "!{params.min_filesize_filtered_assembly}"; then
      echo -e "!{meta.id}\tFiltered Assembly FastA File\tPASS"  \
        >> "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    else
      echo -e "!{meta.id}\tFiltered Assembly FastA File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Filtered_Assembly_File.tsv"
    fi

    bwa index !{uncorrected_contigs}

    msg "INFO: Completed bwa index of !{uncorrected_contigs} FastA assembly file"

    msg "INFO: Cleaned paired-end read mapping of !{meta.id}..."

    bwa mem \
      -v 2 \
      -x intractg \
      -t !{task.cpus} \
      !{uncorrected_contigs} \
      "!{meta.id}_R1.paired.fq.gz" "!{meta.id}_R2.paired.fq.gz" \
      | \
      samtools sort \
      -@ !{task.cpus} \
      -l 9 \
      -o "!{meta.id}-!{meta.assembler}.paired.bam" \
      --reference !{uncorrected_contigs}

    msg "INFO: Completed paired read mapping of !{meta.id}"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.paired.bam" 'Binary PE Alignment Map BAM File' "!{params.min_filesize_binary_pe_alignment}"; then
      echo -e "!{meta.id}\tBinary PE Alignment Map BAM File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
    else
      echo -e "!{meta.id}\tBinary PE Alignment Map BAM File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Binary_PE_Alignment_Map_File.tsv"
    fi

    samtools index "!{meta.id}-!{meta.assembler}.paired.bam"

    msg "INFO: Completed samtools index of paired-end BAM alignment file for !{meta.id}"

    cp -L "!{meta.id}-!{meta.assembler}.uncorrected.fna" "!{meta.id}-!{meta.assembler}.fna"

    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.fna" 'Final Corrected Assembly FastA File' "!{params.min_filesize_final_assembly}"; then
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    else
      echo -e "!{meta.id}\tFinal Corrected Assembly FastA File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Final_Corrected_Assembly_FastA_File.tsv"
    fi

    # Single read mapping if available for downstream depth of coverage calculations
    if [[ !{meta.id}_single.fq.gz ]]; then

      msg "INFO: Single read mapping of !{meta.id}..."

      bwa index "!{meta.id}-!{meta.assembler}.fna"

      msg "INFO: Completed bwa index of !{meta.id}-!{meta.assembler}.fna FastA assembly file"

      bwa mem \
        -v 2 \
        -x intractg \
        "!{meta.id}-!{meta.assembler}.fna" \
        "!{meta.id}_single.fq.gz" \
        -t !{task.cpus} \
        | \
        samtools sort \
        -l 9 \
        -@ !{task.cpus} \
        -o "!{meta.id}-!{meta.assembler}.single.bam" \
        --reference "!{meta.id}-!{meta.assembler}.fna"

      msg "INFO: Completed single read mapping of !{meta.id}"

      echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Binary_SE_Alignment_Map_File.tsv"
      if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.single.bam" 'Binary SE Alignment Map BAM File' '!{params.min_filesize_binary_se_alignment}'; then
        echo -e "!{meta.id}\tBinary SE Alignment Map BAM File\tPASS" \
            >> "!{meta.id}-!{meta.assembler}.Binary_SE_Alignment_Map_File.tsv"
      else
        echo -e "!{meta.id}\tBinary SE Alignment Map BAM File\tFAIL" \
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

    msg "INFO: Calculated checksum of FastA assembly file for !{meta.id}"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bwa: $(bwa 2>&1 | head -n 3 | tail -1 | awk 'NF>1{print $NF}')
        samtools: $(samtools --version | head -n 1 | awk 'NF>1{print $NF}')
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
    END_VERSIONS
    '''
}
