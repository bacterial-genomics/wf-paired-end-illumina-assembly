process POST_BAKTA_QC {
    label "process_high"
    tag { "${meta.id}-${meta.assembler}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(genbank)  // .gbff from bakta

    output:
    tuple val(meta), path("${meta.id}-${meta.assembler}.Annotated_GenBank_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.id}-${meta.assembler}.gbk")                       , emit: bakta_genbank_file
    path("${meta.id}.Annotation_GenBank.SHA512-checksums.tsv")                      , emit: checksums
    path(".command.{out,err}")
    path("versions.yml")                                                            , emit: versions

    script:
    """
    source bash_functions.sh

    # Rename GenBank file from .gbff to .gbk
    cp -f "${genbank}" "${meta.id}-${meta.assembler}.gbk"

    # Perform basic file size check
    echo -e "Sample_name\tQC_step\tOutcome_(Pass/Fail)" > "${meta.id}-${meta.assembler}.Annotated_GenBank_File.tsv"
    if verify_minimum_file_size "${meta.id}-${meta.assembler}.gbk" 'Annotated GenBank File' "${params.min_filesize_annotated_genbank}"; then
        echo -e "${meta.id}\tAnnotated GenBank File\tPASS" >> "${meta.id}-${meta.assembler}.Annotated_GenBank_File.tsv"
    else
        echo -e "${meta.id}\tAnnotated GenBank File\tFAIL" >> "${meta.id}-${meta.assembler}.Annotated_GenBank_File.tsv"
    fi

    # Generate SHA512 checksum of GenBank file with LOCUS line date stripped
    awk '/^LOCUS/ {gsub(/[[:space:]]+[0-9]{2}-[A-Z]{3}-[0-9]{4}/, "", \$0); print} !/^LOCUS/ {print}' "${meta.id}-${meta.assembler}.gbk" \
        | sha512sum \
        | awk -v sample_id="${meta.id}" -v file="${meta.id}.gbk" '
            BEGIN { print "Sample_name\tChecksum_(SHA-512)\tFile" }
            { print sample_id "\t" \$1 "\t" file }
        ' > "${meta.id}.Annotation_GenBank.SHA512-checksums.tsv"

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version | head -n 1)
        sha512sum: \$(sha512sum --version | grep ^sha512sum | sed 's/sha512sum //1')
    END_VERSIONS
    """
}
