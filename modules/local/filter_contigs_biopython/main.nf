process FILTER_CONTIGS_BIOPYTHON {

    // errorStrategy 'terminate'

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${base}.${task.process}${filename}"}

    tag { "${base}" }

    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
        tuple val(base), path(paired_R1_gz), path(paired_R2_gz), path(single_gz), path(qc_overlap_filecheck), path(contigs), path(qc_assembly_filecheck)

    output:
        tuple val(base), path("${base}.uncorrected.fna"), emit: uncorrected_contigs
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
        filter_contigs_params = "-l ${params.filter_contigs_length}"
        if (params.filter_contigs_discard_file != 'None') {
            filter_contigs_params += " -d ${params.filter_contigs_discard_file}"
        }
        if (params.filter_contigs_coverage != 5) {
            filter_contigs_params += " -c ${params.filter_contigs_coverage}"
        }
        if (params.filter_contigs_gcskew == 'True') {
            filter_contigs_params += " -g"
        }
        if (params.filter_contigs_keep_low_complexity == 'True') {
            filter_contigs_params += " -m"
        }
        if (params.filter_contigs_deflines != 'rename_retain') {
            filter_contigs_params += " --deflines ${params.filter_contigs_deflines}"
        }
        if (params.filter_contigs_no_sort == 'True') {
            filter_contigs_params += " --no-sort"
        }
        '''
        # Exit if previous process fails qc filechecks
        if [ $(grep "FAIL" !{base}*File*.tsv) ]; then
          exit 1
        fi

        source bash_functions.sh

        # Get filter.contigs.py and check if it exists
        filter_contigs_script="${DIR}/filter.contigs.py"
        if ! check_if_file_exists_allow_seconds ${filter_contigs_script} '60'; then
        exit 1
        fi

        # Remove junk contigs
        python ${filter_contigs_script} \
        -i !{contigs} \
        -b "!{base}" \
        -o !{base}.uncorrected.fna \
        !{filter_contigs_params}

        # Get process version
        cat <<-END_VERSIONS > versions.yml
        "!{task.process} (!{base})":
            python: $(python --version 2>&1 | awk '{print $2}')
            biopython: $(python -c 'import Bio; print(Bio.__version__)' 2>&1)
        END_VERSIONS
        '''
}