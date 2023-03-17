process POLISH_ASSEMBLY_BWA_PILON {

    // errorStrategy 'terminate'

    publishDir "${params.outpath}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "*.{txt,fna}"
    publishDir "${params.qc_filecheck_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: "*File*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    label "process_high"
    tag { "${base}" }

    container "gregorysprenger/bwa-samtools-pilon@sha256:209ac13b381188b4a72fe746d3ff93d1765044cbf73c3957e4e2f843886ca57f"
    
    input:
        tuple val(base), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_nonoverlap_filecheck), path(uncorrected_contigs)

    output:
        tuple val(base), path("${base}.paired.bam"), path("${base}.single.bam"), path("*File*.tsv"), emit: bam
        tuple val(base), path("${base}.fna"), emit: base_fna
        path "${base}.Filtered_Assembly_File.tsv", emit: qc_filtered_asm_filecheck
        path "${base}.Binary_PE_Alignment_Map_File.tsv", emit: qc_pe_alignment_filecheck
        path "${base}.Polished_Assembly_File.tsv", emit: qc_polished_asm_filecheck
        path "${base}.Final_Corrected_Assembly_FastA_File.tsv", emit: qc_corrected_asm_filecheck
        path "${base}.Binary_SE_Alignment_Map_File.tsv", emit: qc_se_alignment_filecheck
        path "${base}.InDels-corrected.cnt.txt"
        path "${base}.SNPs-corrected.cnt.txt"
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

        # Correct cleaned SPAdes contigs with cleaned PE reads
        if verify_minimum_file_size "!{uncorrected_contigs}" 'Filtered Assembly File' "!{params.min_filesize_non_filtered_assembly}"; then
          echo -e "!{base}\tFiltered Assembly File\tPASS" > !{base}.Filtered_Assembly_File.tsv
        else
          echo -e "!{base}\tFiltered Assembly File\tFAIL" > !{base}.Filtered_Assembly_File.tsv
          exit 1
        fi

        echo -n '' > !{base}.InDels-corrected.cnt.txt
        echo -n '' > !{base}.SNPs-corrected.cnt.txt

        msg "INFO: Correcting contigs with PE reads using !{task.cpus} threads"

        for i in {1..3}; do
          bwa index !{uncorrected_contigs}

          bwa mem \
          -t !{task.cpus} \
          -x intractg \
          -v 2 \
          !{uncorrected_contigs} \
          !{paired_R1_gz} !{paired_R2_gz} \
          | \
          samtools sort \
          -@ !{task.cpus} \
          --reference !{uncorrected_contigs} \
          -l 9 \
          -o !{base}.paired.bam

          if verify_minimum_file_size "!{base}.paired.bam" 'Binary PE Alignment Map File' "!{params.min_filesize_binary_pe_alignment}"; then
            echo -e "!{base}\tBinary PE Alignment Map File (${i} of 3)\tPASS" \
            >> !{base}.Binary_PE_Alignment_Map_File.tsv
          else
            echo -e "!{base}\tBinary PE Alignment Map File (${i} of 3)\tFAIL" \
            >> !{base}.Binary_PE_Alignment_Map_File.tsv
            exit 1
          fi

          samtools index !{base}.paired.bam

          pilon \
          --genome !{uncorrected_contigs} \
          --frags !{base}.paired.bam \
          --output "!{base}" \
          --changes \
          --fix snps,indels \
          --mindepth 0.50 \
          --threads !{task.cpus} >&2

          if verify_minimum_file_size "!{uncorrected_contigs}" 'Polished Assembly File' "!{params.min_filesize_polished_assembly}"; then
            echo -e "!{base}\tPolished Assembly File (${i} of 3)\tPASS" \
            >> !{base}.Polished_Assembly_File.tsv
          else
            echo -e "!{base}\tPolished Assembly File (${i} of 3)\tFAIL" \
            >> !{base}.Polished_Assembly_File.tsv
            exit 1
          fi

          echo $(grep -c '-' !{base}.changes >> !{base}.InDels-corrected.cnt.txt)
          echo $(grep -vc '-' !{base}.changes >> !{base}.SNPs-corrected.cnt.txt)

          rm -f !{base}.{changes,uncorrected.fna}
          rm -f "!{base}"Pilon.bed
          mv -f !{base}.fasta !{base}.uncorrected.fna

          sed -i 's/_pilon//1' !{base}.uncorrected.fna
        done

        mv -f !{base}.uncorrected.fna !{base}.fna

        if verify_minimum_file_size "!{base}.fna" 'Final Corrected Assembly FastA File' "!{params.min_filesize_final_assembly}"; then
          echo -e "!{base}\tFinal Corrected Assembly FastA File\tPASS" \
          > !{base}.Final_Corrected_Assembly_FastA_File.tsv
        else
          echo -e "!{base}\tFinal Corrected Assembly FastA File\tFAIL" \
          > !{base}.Final_Corrected_Assembly_FastA_File.tsv
          exit 1
        fi

        # Single read mapping if available for downstream depth of coverage
        #  calculations, not for assembly polishing.
        if [[ !{single_gz} ]]; then
          msg "INFO: Single read mapping with !{task.cpus} threads"
          bwa index !{base}.fna

          bwa mem \
          -t !{task.cpus} \
          -x intractg \
          -v 2 \
          !{base}.fna \
          !{single_gz} \
          | \
          samtools sort \
          -@ !{task.cpus} \
          --reference !{base}.fna \
          -l 9 \
          -o !{base}.single.bam

          if verify_minimum_file_size "!{base}.single.bam" 'Binary SE Alignment Map File' '!{params.min_filesize_binary_se_alignment}'; then
            echo -e "!{base}\tBinary SE Alignment Map File\tPASS" \
            > !{base}.Binary_SE_Alignment_Map_File.tsv
          else
            echo -e "!{base}\tBinary SE Alignment Map File\tFAIL" \
            > !{base}.Binary_SE_Alignment_Map_File.tsv
            exit 1
          fi

          samtools index !{base}.single.bam
        fi

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            bwa: $(bwa 2>&1 | head -n 3 | tail -1 | awk 'NF>1{print $NF}')
            samtools: $(samtools --version | head -n 1 | awk 'NF>1{print $NF}')
            pilon: $(pilon --version | cut -d ' ' -f 3)
        END_VERSIONS
        '''
}