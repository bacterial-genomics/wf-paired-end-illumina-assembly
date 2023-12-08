process MERGE_LANE_FILES_PYTHON {

    tag { "${samplesheet.getSimpleName()}" }
    container "gregorysprenger/pandas-excel@sha256:4fad4114df25726e24660d8550da48b926b80ce5b8a32b522b552a2d8e1df156"

    input:
    path samplesheet

    output:
    path("lanes_merged_samplesheet.*"), emit: lanes_merged_samplesheet
    path(".command.{out,err}")
    path("versions.yml")              , emit: versions

    shell:
    '''
    samplesheet=!{samplesheet}

    export samplesheet

    python3 <<-END_PYTHON
    import os
    import sys
    import pandas as pd

    # Save filename extension
    extension = ""

    # Read samplesheet based on file extension
    if os.environ['samplesheet'].endswith('.csv'):
        data = pd.read_csv(os.environ['samplesheet'], index_col=None)
        extension = "csv"
    elif os.environ['samplesheet'].endswith('.tsv'):
        data = pd.read_csv(os.environ['samplesheet'], index_col=None, sep='\\t')
        extension = "tsv"

    # Replace spaces in columns with underscores
    data.columns = [c.replace(' ', '_') for c in data.columns]

    # Replace fields containing a line break with a space
    df = data.replace('\\n', ' ', regex=True)

    # Find all duplicated samples
    sample_names = df[df.duplicated(['sample'])]['sample']

    # Make sure theres at least 1 sample that is duplicated
    if len(sample_names) >= 1:
        for name in sample_names:
            # Grab all fastq_1 and fastq_1 columns for each sample_name
            r1 = df[df['sample']==name]['fastq_1']
            r2 = df[df['sample']==name]['fastq_2']

            # Empty string list to add each fastq path to
            r1_list = ""
            r2_list = ""

            # Add each fastq path to r1
            for path in r1:
                r1 += f'{path} '

            # Add each fastq path to r2
            for path in r2:
                r2 += f'{path} '

            os.system(f"cat {r1_list} > {name}_R1.merged.fastq.gz")
            os.system(f"cat {r2_list} > {name}_R2.merged.fastq.gz")

            # Drop data from current dataframe (df)
            df = df[df['sample'] != name]

            merged_r1 = os.path.abspath(f'{name}_R1.merged.fastq.gz')
            merged_r2 = os.path.abspath(f'{name}_R2.merged.fastq.gz')

            if not os.path.exists(merged_r1) or not os.path.exists(merged_r2):
                sys.stderr.write("ERROR: Merged files do not exist.")
                sys.exit(1)

            # Add merged data to dataframe
            df.loc[len(df)] = [name, merged_r1, merged_r2]

    # Write dataframe to file
    if extension == 'tsv':
        df.to_csv('lanes_merged_samplesheet.tsv', sep='\\t', encoding='utf-8', index=False)
    else:
        df.to_csv('lanes_merged_samplesheet.csv', encoding='utf-8', index=False)

    END_PYTHON

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python3 --version 2>&1 | awk '{print $2}')
        ubuntu: $(awk -F ' ' '{print $2, $3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
