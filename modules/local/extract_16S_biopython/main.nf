process EXTRACT_16S_BIOPYTHON {

    tag { "${meta.id}-${meta.assembler}" }
    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(prokka_genbank_file), path(assembly)

    output:
    tuple val(meta), path("16S.${meta.id}-${meta.assembler}.fa"), emit: extracted_rna
    path(".command.{out,err}")
    path("versions.yml")                                        , emit: versions

    shell:
    '''
    source bash_functions.sh

    # 16S extraction
    if [[ -s "!{prokka_genbank_file}" ]]; then
      extract.record.from.genbank.py \
        -i "!{prokka_genbank_file}" \
        -o "16S.!{meta.id}-!{meta.assembler}.fa" \
        -q "!{params.genbank_query}" \
        -f !{params.genbank_query_feature} \
        -u !{params.genbank_query_qualifier} \
        --search-type !{params.genbank_search_type}
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
    END_VERSIONS
    '''
}
