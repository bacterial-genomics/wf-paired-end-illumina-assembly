process ANNOTATE_PROKKA {

    label "process_high"
    tag { "${meta.id}-${meta.assembler}" }
    container "staphb/prokka@sha256:6bb2522e077ef08a8be7a3856fe80372ede3b832becba0728e1ddbe83d89a042"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Annotated_GenBank_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.gbk")                       , emit: prokka_genbank_file
    path("prokka/${meta.id}-${meta.assembler}.log.gz")
    path("${meta.id}.Annotation_GenBank.SHA512-checksums.tsv")                      , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                                            , emit: versions

    shell:
    curated_proteins = params.prokka_curated_proteins ? "--proteins ${params.prokka_curated_proteins}" : ""
    '''
    source bash_functions.sh

    # Remove seperator characters from basename for future processes
    short_base=$(echo !{meta.id} | sed 's/[-._].*//g')
    sed -i "s/!{meta.id}/${short_base}/g" !{assembly}

    # Annotate cleaned and corrected assembly
    msg "INFO: Annotating assembly using Prokka for !{meta.id}..."

    # Run Prokka
    prokka \
      --outdir prokka \
      --prefix "!{meta.id}-!{meta.assembler}" \
      --locustag "!{meta.id}-!{meta.assembler}" \
      --evalue !{params.prokka_evalue} \
      --mincontiglen 1 \
      --force \
      --addgenes \
      !{curated_proteins} \
      --cpus !{task.cpus} \
      !{assembly}

    msg "INFO: Completed annotating the assembly using Prokka for !{meta.id}..."

    # Regardless of the file extension, unify to GBK extension for GenBank format
    for ext in gb gbf gbff gbk; do
      if [ -s "prokka/!{meta.id}-!{meta.assembler}.${ext}" ]; then
        mv -f "prokka/!{meta.id}-!{meta.assembler}.${ext}" \
          "!{meta.id}-!{meta.assembler}.gbk"
        break
      fi
    done

    # Verify output file
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "!{meta.id}-!{meta.assembler}.Annotated_GenBank_File.tsv"
    if verify_minimum_file_size "!{meta.id}-!{meta.assembler}.gbk" 'Annotated GenBank File' "!{params.min_filesize_annotated_genbank}"; then
      echo -e "!{meta.id}-!{meta.assembler}\tAnnotated GenBank File\tPASS" \
        >> "!{meta.id}-!{meta.assembler}.Annotated_GenBank_File.tsv"
    else
      echo -e "!{meta.id}-!{meta.assembler}\tAnnotated GenBank File\tFAIL" \
        >> "!{meta.id}-!{meta.assembler}.Annotated_GenBank_File.tsv"
    fi

    # Compress the bulky verbose logfile for compact storage
    gzip -9f "prokka/!{meta.id}-!{meta.assembler}.log"

    # Calculate checksum
    FILE="!{meta.id}-!{meta.assembler}.gbk"
    CHECKSUM=$(awk '/^LOCUS/ {gsub(/[[:space:]]+[0-9]{2}-[A-Z]{3}-[0-9]{4}/, "", $0); print} !/^LOCUS/ {print}' "${FILE}" | sha512sum | awk '{print $1}')
    echo "${CHECKSUM}" | awk -v sample_id="!{meta.id}" -v file="${FILE}" '
        BEGIN {
            # Print the header once
            print "Sample_name\tChecksum_(SHA-512)\tFile"
        }
        {
            # Print the data row once, using the CHECKSUM from input
            print sample_id "\t" $1 "\t" file
        }' \
        > "!{meta.id}.Annotation_GenBank.SHA512-checksums.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        prokka: $(prokka --version 2>&1 | awk 'NF>1{print $NF}')
        sha512sum: $(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
    END_VERSIONS
    '''
}
