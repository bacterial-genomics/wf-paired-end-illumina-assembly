process CREATE_EXCEL_RUN_SUMMARY_PYTHON {

    container "gregorysprenger/pandas-excel@sha256:4fad4114df25726e24660d8550da48b926b80ce5b8a32b522b552a2d8e1df156"

    input:
    path(list_of_files)

    output:
    path("Summary-Report_*.xlsx") , emit: summary
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    python3 <<-END_PYTHON
    import glob
    import datetime
    import pandas as pd

    def create_summary_workbook(output_file, tsv_file):
        sheet_name = tsv_file.split(".")[1]
        data = pd.read_csv(tsv_file, sep="\t")
        data.to_excel(output_file, sheet_name=sheet_name, index=False)

    date = datetime.datetime.now()
    date_format = date.strftime("%Y-%b-%d_%H-%M-%S")

    list_of_files = glob.glob("*.tsv")

    with pd.ExcelWriter(f"Summary-Report_{date_format}.xlsx") as output_file:
        for file in list_of_files:
            create_summary_workbook(output_file, file)

    END_PYTHON

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        ubuntu: \$(awk -F ' ' '{print \$2, \$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
