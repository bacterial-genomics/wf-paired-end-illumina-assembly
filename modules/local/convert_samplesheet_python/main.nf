process CONVERT_SAMPLESHEET_PYTHON {

    container "gregorysprenger/pandas-excel@sha256:4fad4114df25726e24660d8550da48b926b80ce5b8a32b522b552a2d8e1df156"

    input:
    path excel_samplesheet

    output:
    path("samplesheet.tsv")   , emit: converted_samplesheets
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    samplesheet="!{excel_samplesheet}"
    sheet_name="!{params.excel_sheet_name}"

    export samplesheet sheet_name

    python3 <<-END_PYTHON
    import os
    import pandas as pd

    # Read Excel/LibreCalc file
    data = pd.read_excel(os.environ['samplesheet'], sheet_name=os.environ['sheet_name'], index_col=None)

    # Replace spaces in columns with underscores
    data.columns = [c.replace(' ', '_') for c in data.columns]

    # Replace fields containing a line break with a space
    df = data.replace('\\n', ' ', regex=True)

    # Write dataframe to TSV
    df.to_csv('samplesheet.tsv', sep='\\t', encoding='utf-8', index=False)
    END_PYTHON

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python3 --version 2>&1 | awk '{print $2}')
        ubuntu: $(awk -F ' ' '{print $2, $3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
