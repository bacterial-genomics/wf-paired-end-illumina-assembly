# wf-paired-end-illumina-assembly: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Annotated genbank file](#annot) - Genbank annotated file
- [Final assembly file](#asm) - Final, corrected assembly file
- [SPAdes output](#spades-output) - Contigs and SPAdes log files
- [QA](#qa) - QA information on assembly files
- [Small subunit (16S) files](#ssu) - BLAST output files
- [Cleaned reads](#trim-reads) - Cleaned reads
- [Summaries](#summaries) - Output summaries
- [Log files](#log) - Nextflow and HPC logs, software information, and error list if applicable
- [Process logs](#process-logs) - Output and error logs for each process
- [QC file checks](#qc-file-checks) - Process output quality checks to determine pipline continuation

# Output File Structure

_Note: Output file structure is based on the output path given to `--outdir`._

_Note: \<SampleName\> is parsed from input FastQ filenames and excludes everything after \_{R1,R2}._

_Note: \<Assembler\> is the name of the assembler (SPAdes, SKESA) given to `--assembler`. \[Default: SPAdes\]_

| Output Directory                                         | Filename                                               | Explanation                                                                            |
| -------------------------------------------------------- | ------------------------------------------------------ | -------------------------------------------------------------------------------------- |
| <a id="annot">Annotation</a>                             |                                                        | **Annotation files**                                                                   |
| Annotation/Prokka                                        | \<SampleName\>-\<Assembler\>.gbk                       | Genbank file created by Prokka                                                         |
| <a id="asm">Assembly</a>                                 |                                                        | **Assembly files**                                                                     |
|                                                          | \<SampleName\>-\<Assembler\>.fna                       | Corrected assembly                                                                     |
| <a id="spades-output">Assembly/SPAdes/\<SampleName\></a> |                                                        | **SPAdes Assembly files**                                                              |
|                                                          | contigs.fasta                                          | Assembled contigs from SPAdes                                                          |
|                                                          | params.txt.gz                                          | Parameters used with SPAdes                                                            |
|                                                          | spades.log.gz                                          | Log information from SPAdes                                                            |
|                                                          | assembly_graph_with_scaffolds.gfa                      | Contains SPAdes assembly graph and scaffolds paths                                     |
|                                                          | \<SampleName\>-\<Assembler\>.InDels-corrected.cnt.txt  | Each line represents number of corrected InDels (per correction round)                 |
|                                                          | \<SampleName\>-\<Assembler\>.SNPs-corrected.cnt.txt    | Each line represents number of corrected SNPs (per correction round)                   |
|                                                          | \<Num\>of3-asm-attempt-failed.spades.log               | Output log of SPAdes for each attempt tried (up to 3)                                  |
| <a id="qa">Assembly/QA</a>                               |                                                        | **Quality Assurance files**                                                            |
| Assembly/QA/\<SampleName\>                               |                                                        | Output directory for each \<SampleName\>                                               |
|                                                          | \<SampleName\>-\<Assembler\>.CleanedReads-Bases.tsv    | Number of cleaned bases                                                                |
|                                                          | \<SampleName\>-\<Assembler\>.GenomeCoverage.tsv        | Genome Coverage                                                                        |
|                                                          | \<SampleName\>-\<Assembler\>.QuastSummary.tsv          | Quast Summary                                                                          |
| <a id="ssu">SSU</a>                                      |                                                        | **Small Subunit (16S) files**                                                          |
|                                                          | 16S-top-species.tsv                                    | Top BlAST hit results                                                                  |
|                                                          | 16S.\<SampleName\>-\<Assembler\>.fa                    | Top BLAST hit in FastA format                                                          |
| SSU/BLAST                                                |                                                        | BLAST output files                                                                     |
|                                                          | \<SampleName\>-\<Assembler\>.blast.tsv.gz              | BLAST output                                                                           |
| <a id="trim-reads">CleanedReads</a>                      |                                                        | **Trimmed Reads**                                                                      |
|                                                          | \<SampleName\>\_{R1,R2}.paired.fq.gz                   | Cleaned paired reads                                                                   |
|                                                          | \<SampleName\>.single.fq.gz                            | Cleaned single read                                                                    |
| CleanedReads/FLASH                                       |                                                        | FLASH Output                                                                           |
|                                                          | \<SampleName\>.clean-reads.tsv                         | Number of cleaned reads                                                                |
|                                                          | \<SampleName\>.overlap.tsv                             | Number of overlapping reads                                                            |
| CleanedReads/Trimmomatic                                 |                                                        | Trimmomatic Output                                                                     |
|                                                          | \<SampleName\>.trimmomatic.tsv                         | Discarded reads and singletons                                                         |
| <a id="summaries">Summaries</a>                          |                                                        | **Output Summaries**                                                                   |
|                                                          | Summary.16S.tab                                        | Top BLAST hit results                                                                  |
|                                                          | Summary.Assemblies.tab                                 | Contig summary information                                                             |
|                                                          | Summary.BUSCO.tab                                      | Intra-contig gene information using BUSCO                                              |
|                                                          | Summary.GTDB-Tk.tab                                    | Contig taxonomic classification using GTDB-Tk                                          |
|                                                          | Summary.CleanedReads-AlnStats.tab                      | Basepairs of Paired Reads and Singnleton Reads mapped                                  |
|                                                          | Summary.CleanedReads-Bases.tab                         | Number of cleaned bases                                                                |
|                                                          | Summary.GenomeCoverage.tab                             | Genome Coverage                                                                        |
|                                                          | Summary.MLST.tab                                       | MLST result                                                                            |
|                                                          | Summary.PhiX.tsv                                       | Information on the removal of PhiX reads                                               |
|                                                          | Summary.QC_File_Checks.tab                             | QC file checks                                                                         |
| <a id="taxonomy">Taxonomy</a>                            |                                                        | **Taxonomic Classification**                                                           |
| Taxonomy/kraken/\<SampleName\>                           |                                                        | Kraken Output for each \<SampleName\>                                                  |
|                                                          | \<SampleName\>.kraken_output.tab.gz                    | Full Kraken output                                                                     |
|                                                          | \<SampleName\>.kraken_summary.tsv                      | Summarized unclassified, top 2 genus and top 2 species information                     |
| Taxonomy/kraken2/\<SampleName\>                          |                                                        | Kraken2 Output for each \<SampleName\>                                                 |
|                                                          | \<SampleName\>.kraken_output.tab.gz                    | Full Kraken2 output                                                                    |
|                                                          | \<SampleName\>.kraken2_summary.tsv                     | Summarized unclassified, top 2 genus and top 2 species information                     |
| <a id="log">pipeline_info</a>                            |                                                        | **Log files**                                                                          |
|                                                          | ASM\_\<Number of Samples\>.o\<Submission Number\>      | HPC output report                                                                      |
|                                                          | ASM\_\<Number of Samples\>.e\<Submission Number\>      | HPC error report                                                                       |
|                                                          | pipeline_dag.\<YYYY-MM-DD_HH-MM-SS\>.html              | Direct acrylic graph of workflow                                                       |
|                                                          | report.\<YYYY-MM-DD_HH-MM-SS\>.html                    | Nextflow summary report of workflow                                                    |
|                                                          | timeline.\<YYYY-MM-DD_HH-MM-SS\>.html                  | Nextflow execution timeline of each process in workflow                                |
|                                                          | trace.\<YYYY-MM-DD_HH-MM-SS\>.txt                      | Nextflow execution tracing of workflow, which includes percent of CPU and memory usage |
|                                                          | software_versions.yml                                  | Versions of software used in each process                                              |
|                                                          | errors.tsv                                             | Errors file if errors exist and summarizes the errors                                  |
| <a id="process-logs">pipeline_info/process_logs</a>      |                                                        | **Process log files**                                                                  |
|                                                          | \<SampleName\>.\<ProcessName\>.command.out             | Standard output for \<SampleName\> during process \<ProcessName\>                      |
|                                                          | \<SampleName\>.\<ProcessName\>.command.err             | Standard error for \<SampleName\> during process \<ProcessName\>                       |
| <a id="qc-file-checks">pipeline_info/qc_file_checks</a>  |                                                        | **QC file check log files**                                                            |
|                                                          | \<SampleName\>.Raw_Initial_FastQ_Files.tsv             | Raw Initial FastQ File Check                                                           |
|                                                          | \<SampleName\>.PhiX_Genome.tsv                         | PhiX Genome Check                                                                      |
|                                                          | \<SampleName\>.PhiX-removed_FastQ_Files.tsv            | PhiX-removed FastQ File Check                                                          |
|                                                          | \<SampleName\>.Adapters_FastA.tsv                      | Adapters FastA File Check                                                              |
|                                                          | \<SampleName\>.Adapter-removed_FastQ_Files.tsv         | Adapter-removed FastQ File Check                                                       |
|                                                          | \<SampleName\>.Non-overlapping_FastQ_Files.tsv         | Non-overlapping FastQ File Check                                                       |
|                                                          | \<SampleName\>.Raw_Assembly_File.tsv                   | Raw Assembly File Check                                                                |
|                                                          | \<SampleName\>.Filtered_Assembly_File.tsv              | Filtered Assembly File Check                                                           |
|                                                          | \<SampleName\>.Binary_PE_Alignment_Map_File.tsv        | Binary Paired-End Alignment Map File Check                                             |
|                                                          | \<SampleName\>.Polished_Assembly_File.tsv              | Polished Assembly File Check                                                           |
|                                                          | \<SampleName\>.Final_Corrected_Assembly_FastA_File.tsv | Final Corrected Assembly FastA File Check                                              |
|                                                          | \<SampleName\>.Binary_SE_Alignment_Map_File.tsv        | Binary Singletons Alignment Map File Check                                             |
|                                                          | \<SampleName\>.Annotated_GenBank_File.tsv              | Annotated GenBank File Check                                                           |
|                                                          | \<SampleName\>.SSU_Extracted_File.tsv                  | SSU Extracted File Check                                                               |
|                                                          | \<SampleName\>.16S_BLASTn_Output_File.tsv              | 16S BLASTn Output File Check                                                           |
|                                                          | \<SampleName\>.Filtered_16S_BLASTn_File.tsv            | Filtered 16S BLASTn File Check                                                         |
