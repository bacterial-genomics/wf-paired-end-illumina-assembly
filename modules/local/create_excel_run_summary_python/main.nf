process CREATE_EXCEL_RUN_SUMMARY_PYTHON {

    container "gregorysprenger/pandas-excel@sha256:4fad4114df25726e24660d8550da48b926b80ce5b8a32b522b552a2d8e1df156"

    input:
    path(list_of_files)
    path(tab_colors)
    val(wf_version)

    output:
    path("Summary-Report.xlsx"), emit: summary
    path(".command.{out,err}")
    path("versions.yml")       , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Converting summary TSV files: !{list_of_files} into single XLSX workbook..."

    # Default outfile is "Summary-Report.xlsx" in python script
    tsv_to_excel.py \
      !{list_of_files} \
      --color-dict !{tab_colors} \
      -t "wf-paired-end-illumina-assembly_v!{wf_version}" \
      -s "genome_assembly" \
      -c "bacterial-genomics"

    msg "INFO: Converted summary TSV files into single XLSX workbook."

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        pandas: $(python3 -c "import pandas as pd; print(pd.__version__)")
        python: $(python3 --version 2>&1 | awk '{print $2}')
        ubuntu: $(awk -F ' ' '{print $2, $3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
