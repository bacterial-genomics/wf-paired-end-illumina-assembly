# wf-paired-end-illumina-assembly: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - June 2, 2023

Converted workflow to nf-core style and added SKESA assembler.

### `Added`
- nf-core Styling
- SKESA
- Allow samplesheet (XLSX, CSV, TSV) OR directory as input

### `Fixed`
- Errors in run_assembly.uge-nextflow scripts

### `Updated`
- Docker containers to ensure they all have built in tests

### `Deprecated`

## v1.0.0 - January 20, 2023

Initial release of wf-paired-end-illumina-assembly.