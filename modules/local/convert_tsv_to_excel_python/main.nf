process CONVERT_TSV_TO_EXCEL_PYTHON {

    container "gregorysprenger/pandas-excel@sha256:4fad4114df25726e24660d8550da48b926b80ce5b8a32b522b552a2d8e1df156"

    input:
    path("summary_files/")

    output:
    path(".command.{out,err}")
    path("versions.yml")      , emit: versions

    shell:
    '''
    output_dir="!{projectDir}/!{params.outdir}"
    export output_dir

    python3 <<-END_PYTHON
    import os
    import glob
    import pandas as pd

    def convert_tsv_to_excel(tsv_file):
        path = os.path.dirname(tsv_file)
        basename = os.path.basename(tsv_file).split(".tsv")[0]
        sheet_name = basename.split(".")[-1]

        data = pd.read_csv(tsv_file, sep="\\t")
        data.to_excel(f"{path}/{basename}.xlsx", sheet_name=sheet_name, index=False)

    results_dir = os.environ["output_dir"]
    all_tsv_files = glob.glob(f"{results_dir}/**/*.tsv", recursive=True)

    # Drop all TSV files in pipeline_info directory
    list_of_files = [file for file in all_tsv_files if "pipeline_info" not in file]

    for file in list_of_files:
        convert_tsv_to_excel(file)

    END_PYTHON

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python3 --version 2>&1 | awk '{print $2}')
        ubuntu: $(awk -F ' ' '{print $2, $3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
