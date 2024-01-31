# wf-paired-end-illumina-assembly: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and is used to perform *de novo* assembly on raw Illumina paired-end reads from bacterial isolates.

- [Input quality control](#input-quality-control) that includes trimming and contaminant removal
  - [Initial FastQ file check](#initial-fastq-file-check) ensures input FastQ files meet a minimum file size
  - [Host read removal](#host-read-removal)
  - [PhiX read removal](#phix-read-removal)
  - [Adapter clipping and quality trimming](#adapter-clipping-and-quality-trimming)
  - [Merge overlapping sister reads](#merge-overlapping-sister-reads) to create singletons
- [Taxonomic classification of trimmed reads](#taxonomic-classification-of-trimmed-reads)
  - [Kraken](#kraken)
  - [Kraken2](#kraken2)
- [Assembly](#assembly) of trimmed reads
  - [SPAdes](#spades)
  - [SKESA](#skesa)
- [Assembly metrics and classification](#assembly-metrics-and-classification)
  - [Quality assessment and coverage](#quality-assessment-and-coverage)
  - [Multilocus sequence typing (MLST)](#multilocus-sequence-typing-mlst)
  - [Genome annotation](#genome-annotation)
  - [16S ribosomal RNA (rRNA) classification](#16s-ribosomal-rna-rrna-classification)
  - [Assembly taxonomic classification](#assembly-taxonomic-classification) using GTDB-Tk
- [Summaries](#summaries) of the output files generated during the pipeline process
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

> [!NOTE]
> `[sample]` is a unique identifier that is parsed from input FastQ filenames and excludes everything after [R1/R2].
>
> `[assembler]` is the name of the assembler used to assemble contigs. [Default: SPAdes].

## Input quality control

Input files must meet a minimum file size to be processed within this pipeline. If this check passes, the input files go through host removal, down sampling, PhiX read removal, and adapter trimming.

### Initial FastQ file check

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- Input files are checked to ensure that they meet a minimum file size to be processed within this pipeline `[Default: 25M]`. This is to prevent unusually small input sets from wasting compute time processing data that will not yield a usable bacterial genome assembly.

</details>

### Host read removal

Host read removal can be skipped or performed by Hostile and/or NCBI SRA-Human-Scrubber by specifying `--host_remove {both,hostile,sra-human-scrubber,skip}`. For SRA-Human-Scrubber, reads are repaired using BBTools to discard broken sister reads. Information about the number of reads discarded and retained are saved in the output directory. Please see the [host removal using Hostile documentation](../modules/local/remove_host_hostile/README.md) and [host removal using SRA-Human-Scrubber documentation](../modules/local/remove_host_sra_human_scrubber/README.md) for more information.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- FastQ files after host removal are checked to ensure that they meet a minimum file size before continuing to downstream processes `[Default: 25M]`.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `CleanedReads/Hostile/`
  - `[sample].Summary.Hostile-Removal.tsv`: Summary of the number of reads discarded and retained from Hostile.

- `CleanedReads/SRA-Human-Scrubber/`
  - `[sample].Summary.BBTools-Repair-Removal.tsv`: Summary of the number of reads discarded and retained after repairing broken sister reads produced from SRA-Human-Scrubber.
  - `[sample].Summary.SRA-Human-Scrubber-Removal.tsv`: Summary of the number of reads discarded and retained from SRA-Human-Scrubber.

</details>

### PhiX read removal

PhiX reads are commonly used as a positive control for Illumina sequencing. During assembly, PhiX reads are considered contaminants and if retained, a misassembled genome will be formed. Therefore, a PhiX reference file is required and a default [PhiX reference file](../bin/PhiX_NC_001422.1.fasta) is included with this pipeline. Please see the [PhiX removal using  BBDuk documentation](../modules/local/remove_phix_bbduk/README.md) for more information.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- The PhiX reference genome that is used with BBDuk to remove sequence reads must be accessible and meet a minimum file size `[Default: 5k]`.

- FastQ files after PhiX removal are checked to ensure that they meet a minimum file size before continuing to downstream processes `[Default: 25M]`. This is to halt analysis of a sample mostly containing PhiX reads rather than the sample DNA itself.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `CleanedReads/BBDUK/`
  - `[sample].Summary.PhiX.tsv`: Number of reads discarded and retained from BBDuk.

</details>

### Adapter clipping and quality trimming

Illumina instruments can detect and remove adapter sequences, but sometimes adapters can end up in the FastQ output due to sequencing errors. A default [adapters reference file](../bin/adapters_Nextera_NEB_TruSeq_NuGEN_ThruPLEX.fas) is used within this pipeline. Trimmomatic also performs quality trimming, where broken sister reads are retained for downstream processes. Please see the [adapter clipping and quality trimming using Trimmomatic documentation](../modules/local/trim_reads_trimmomatic/README.md) for more information.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- The adapters reference file that is used with Trimmomatic to remove sequence reads must be accessible and meet a minimum file size `[Default: 10k]`.

- FastQ files after removing adapter sequences are checked to ensure that they meet a minimum file size before continuing to downstream processes `[Default: 25M]`. This is to halt analysis of a sample that is primarily contaminated with artificial Illumina sequences.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `CleanedReads/Trimmomatic/`
  - `[sample].trimmomatic.tsv`: Summary of the number of reads discarded and retained from Trimmomatic.

</details>

### Merge overlapping sister reads

Overlapping content between sister reads that are at least 80% similar are collapsed into a singleton read. Please see the [overlapping of paired-end reads documentation](../modules/local/overlap_paired_reads_flash/README.md) for more information.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- After singletons have been formed from overlapping read pairs, the remaining non-overlapping paired-end (R1,R2) read files must meet a minimum file size `[Default: 20M]`. This is to prevent the analysis of paired-end read files that are heavily misconstructed into small fragment sizes.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `CleanedReads/`
  - `[sample]_single.fq.gz`: Final cleaned singleton reads.
  - `[sample]_R[1/2].paired.fq.gz`: Final cleaned paired reads.

- `CleanedReads/FLASH/`
  - `[sample].overlap.tsv`: Number of reads that were overlapped into singleton reads.
  - `[sample].clean-reads.tsv`: Number of non-overlapping reads.

</details>

## Taxonomic classification of trimmed reads

These classifiers perform classifications on a read-by-read basis or through the use of k-mers on the cleaned and trimmed FastQ files. The results that are obtained are heavily dependent on the quality and diversity of the database used. Therefore, the results produced from these classifiers should be used as a quality check to evaluate the possibility of sample contamination.

> [!WARNING]
> Taxonomic classification tools will be skipped if the accompanying database is not specified.

### Kraken

Kraken is a k-mer based classification tool that assigns taxonomic labels using the Lowest Common Ancestor (LCA) algorithm.

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/kraken/[sample]/`
  - `[sample].kraken_summary.tsv`: Summary of the unclassified, top 3 genus and top 3 species classifications from the Kraken report.
  - `[sample].kraken_output.tab.gz`: Taxonomic classification in the Kraken report format.

</details>

### Kraken2

Kraken2 is a k-mer based classification tool that assigns taxonomic labels using the Lowest Common Ancestor (LCA) algorithm.

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/kraken2/[sample]/`
  - `[sample].kraken2_summary.tsv`: Summary of the unclassified, top 3 genus and top 3 species classifications from the Kraken report.
  - `[sample].kraken2_output.tab.gz`: Taxonomic classification in the Kraken report format.

</details>

## Assembly

The cleaned and trimmed reads are used to assemble contigs using SPAdes or SKESA `[Default: SPAdes]`. Contigs that have low compositional complexity are discarded. Please see the [contig filtering documentation](../modules/local/filter_contigs_biopython/README.md) for more information. Contigs from SPAdes require polishing to create the final genome assembly, which is done by using BWA, Pilon, and Samtools. Contigs from SKESA do not require this step.

> [!IMPORTANT]
> Outputs generated by SPAdes and SKESA cannot be compared even when using the same FastQ inputs.

> [!TIP]
> For many input FastQ files, SKESA may be useful in decreasing runtime. For input FastQ files that may be heavily contaminated, SPAdes may help maintain contiguity.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- The contigs produced from an assembler software package must meet a minimum file size criteria `[Default: 1M]`. This is to prevent the analysis of a highly incomplete bacterial genome.

- The resulting contig file after low compositional complexity contigs are discarded must meet a minimum file size `[Default: 1M]`. This is to prevent the analysis of a highly incomplete bacterial genome.

- The cleaned paired-end reads are mapped onto the filtered assembly file in sequential steps (`[Default: 3]`), and the resulting binary paired-end alignment file must meet a minimum file size criteria `[Default: 25M]`. This is to prevent the analysis of an assembly file that has an unusually low read sequence amount.

- The assembly file goes through SNP and InDel corrections in sequential steps (`[Default: 3]`), and the resulting assembly file must meet a minimum file size criteria `[Default: 1M]`. This is to prevent further analysis of an unusually incomplete genome.

- The final error-corrected assembly file must meet a minimum file size criteria `[Default: 1M]`. This is to ensure that the final assembly file is not unexpectedly small or incomplete.

- If singletons (single-end reads) exist after read cleaning, they are mapped onto the assembly file and the resulting binary single-end alignment file must meet a minimum file size criteria `[Default: 1k]`. This is to ensure that read depth calculations can be performed on the single-end reads.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `Assembly/`
  - `[sample]-[assembler]_contigs.fna`: Final genome assembly.

</details>

### SPAdes

SPAdes is a k-mer based software that forms a genome assembly from read sequences. Contigs from SPAdes require polishing to create the final genome assembly, which is done by using BWA, Pilon, and Samtools.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/SPAdes/[sample]/`
  - `[sample]-SPAdes.log.gz`: SPAdes log file.
  - `[sample]-SPAdes_graph.gfa`: Assembly graph in gfa format.
  - `[sample]-SPAdes_warnings.log`: Log file that lists warnings when forming a genome assembly.
  - `[sample]-SPAdes_params.txt.gz`: Command used to perform the SPAdes analysis.
  - `[sample]-SPAdes_contigs.fasta`: Assembled contigs in FastA format.
  - `[sample]-SPAdes.SNPs-corrected.cnt.txt`: Number of SNPs corrected in each round of corrections.
  - `[sample]-SPAdes.InDels-corrected.cnt.txt`: Number of InDels corrected in each round of corrections.

</details>

### SKESA

Strategic K-mer Extension for Scrupulous Assemblies (SKESA) is a software that is based on DeBruijn graphs that forms a genome assembly from read sequences.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/SKESA/[sample]/`
  - `[sample]-SKESA_contigs.fasta`: Assembled contigs in FastA format.

</details>

## Assembly metrics and classification

### Quality assessment and coverage

QUAST is used to perform quality assessment on the assembly file to report metrics such as N50, cumulative length, longest contig length, and GC composition. Using bedtools, the BAM file created with the assembly is used to calculate genome size and the coverage of paired-end and single-end reads.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/QA/[sample]/`
  - `[sample]-[assembler].QuastSummary.tsv`: Assembly metrics such as N50, cumulative length, longest contig length, and GC composition.
  - `[sample]-[assembler].GenomeCoverage.tsv`: Genome coverage information.
  - `[sample]-[assembler].CleanedReads-Bases.tsv`: Number of cleaned bases.

</details>

### Multilocus sequence typing (MLST)

The final assembly file is scanned against PubMLST typing schemes to determine the MLST for each sample. Unless a specific typing scheme is specified, the best typing scheme for each sample is automatically selected.

<details markdown="1">
<summary>Output files</summary>

- `Summaries/`
  - `Summary.MLST.tab`: Summary of the MLST results for all samples.

</details>

<br />

<details markdown="1">
<summary>MLST output interpretation</summary>

Symbol | Meaning                               | Length          | Identity
-------|---------------------------------------|-----------------|---------------
`n`    | exact intact allele                   | 100%            | 100%
`~n`   | novel full length allele similar to n | 100%            | &ge; `--minid`
`n?`   | partial match to known allele         | &ge; `--mincov` | &ge; `--minid`
`-`    | allele missing                        | &lt; `--mincov` | &lt; `--minid`
`n,m`  | multiple alleles                      |                 |

</details>

### Genome annotation

The final assembly file is annotated to identify and label features using Prokka. Please see [genome annotation using Prokka documentation](../modules/local/annotate_prokka/README.md) for more information.

> [!IMPORTANT]
> Gene symbols may not be as update as the product description when using Prokka. It is recommended to use NCBI's [Prokaryotic Genome Annotation Pipeline (PGAP)](https://github.com/ncbi/pgap) to further investigate putatively present genes.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- The resulting annotated GenBank file must meet a minimum file size criteria `[Default: 3M]`. This is to prevent further analysis of a highly incomplete annotation set.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `Annotation/Prokka/`
  - `[sample]-[assembler].gbk`: Annotated genome in GenBank file format.

</details>

### 16S ribosomal RNA (rRNA) classification

The GenBank file is parsed for 16S rRNA gene records. If there are no 16S rRNA gene records, Barrnap is used to predict 16S rRNA genes using the assembly file. BLAST is then used to align these gene records to its database, where the best alignment is filtered out based on bit score.

> [!NOTE]
> Some assembled genomes do not contain identifiable 16S rRNA sequences and therefore 16S is not able to be classified. If the classification of 16S rRNA sequences is required, the sample must be re-sequenced.

> [!IMPORTANT]
> The 16S rRNA classification produced should not be used as a definitive classification as some taxa have 16S sequences that are extremely similar between different species.
>
> For an exact species match, 100% identity and 100% alignment are needed.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- The extracted 16S rRNA gene sequence file from the genome assembly must meet a minimum file size criteria `[Default: 500b]`. This is to prevent the classification of a highly incomplete 16S rRNA gene sequence.

- The sample identifiers in the 16S rRNA gene sequence file gets renamed and this resulting file must meet a minimum file size criteria `[Default: 500b]`. This is to prevent the classification of a highly incomplete 16S rRNA gene sequence.

- The 16S rRNA gene sequence file is classified using BLASTn and the resulting output file must meet a minimum file size criteria `[Default: 10b]`. This is to ensure that the BLASTn output file contains at least one alignment with taxonomic information to be filtered and reported.

- The best 16S BLASTn alignment sequence filtered by bit score is saved to a file in FastA format and must meet a minimum file size criteria `[Default: 10b]`. This is to ensure that the resulting file contains alignment information that can be parsed and reported.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `SSU/`
  - `16S-top-species.tsv`: Summary of the best BLAST alignment for each sample.
  - `16S.[sample]-[assembler].fa`: 16S rRNA gene sequence of the best BLAST alignment in FastA format.

- `SSU/BLAST/`
  - `[sample]-[assembler].blast.tsv.gz`: BLAST output 16S rRNA gene records in tab-separated value (TSV) format.

</details>

### Assembly taxonomic classification

GTDB-Tk is a taxonomic classification tool that uses the Genome Database Taxonomy [GTDB](https://gtdb.ecogenomic.org/). GTDB-Tk can be used as a quality check for the assembly file.

<details markdown="1">
<summary>Output files</summary>

- `Summaries/`
  - `Summary.GTDB-Tk.tab`: Summary of the GTDB-Tk taxonomic classification for each sample.

</details>

## Summaries

Concatenation of output metrics for all samples.

<details markdown="1">
<summary>Output files</summary>

- `Summaries/`
  - `Summary.16S.tab`: Summary of the best BLAST alignment for each sample.
  - `Summary.MLST.tab`: Summary of the MLST results for all samples.
  - `Summary.PhiX.tab`: Number of reads discarded and retained for each sample.
  - `Summary.Assemblies.tab`: Assembly metrics such as N50, cumulative length, longest contig length, and GC composition for each sample.
  - `Summary.GenomeCoverage.tab`: Summary of the overall genome coverage for each sample.
  - `Summary.QC_File_Checks.tab`: Summary of all QC file checks detailing if a sample passes or fails each process.
  - `Summary.CleanedReads-Bases.tab`: Summary of the number of cleaned bases for each sample.
  - `Summary.CleanedReads-AlignmentStats.tab`: Summary of the genome size and coverages of the paired-end and single-end reads for each sample.

</details>

## Pipeline information

Information about the pipeline execution, output logs, error logs, and QC file checks for each sample are stored here.

> [!NOTE]
> Pipeline execution files have the date and time appended to the filename using the following shorthand notation: year (yyyy), month (MM), day (dd), hour (HH), minute (mm), second (ss).

<details markdown="1">
<summary>Pipeline information</summary>

- `pipeline_info/`
  - `software_versions.yml`: Summary of the software packages used in each process and their version information.
  - `nextflow_log.[job_id].txt`: Execution log file produced by Nextflow.
  - `ASM_[num_of_samples].o[job_id]`: Output log file produced by the job scheduler.
  - `ASM_[num_of_samples].e[job_id]`: Error log file produced by the job scheduler.
  - `pipeline_dag_yyyy-MM-dd_HH-mm-ss.html`: Direct acrylic graph (DAG) image of the workflow that gives a visual representation of how each process connects to each other.
  - `execution_trace_yyyy-MM-dd_HH-mm-ss.txt`: Text-based summary report detailing the work directory hash, runtime, CPU usage, memory usage, etc. for each process.
  - `execution_report_yyyy-MM-dd_HH-mm-ss.html`: Summary report of all processes, including processes that passed/failed, resource usage, etc.
  - `execution_timeline_yyyy-MM-dd_HH-mm-ss.html`: Summary report detailing the runtime and memory usage of each process.
  -

</details>

<details markdown="1">
<summary>Process log information</summary>

- `pipeline_info/process_logs/`
  - `[sample].[process].command.out`: Output log file for each sample in each process.
  - `[sample].[process].command.err`: Error log file for each sample in each process.

</details>

<details markdown="1">
<summary>QC file checks</summary>

- `pipeline_info/qc_file_checks/`
  - `[sample].Raw_Initial_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria for the pipeline `[Default: 25M]`.
  - `[sample].Summary.Hostile-Removal.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria for after host removal using Hostile `[Default: 25M]`.
  - `[sample].SRA_Human_Scrubber_FastQ_File.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria for after host removal using SRA-Human-Scrubber `[Default: 25M]`.
  - `[sample].BBTools-Repair-removed_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria after repairing broken sister reads from SRA-Human-Scrubber `[Default: 25M]`.
  - `[sample].PhiX_Genome.tsv`: Details if the input PhiX reference genome meets the minimum file size criteria `[Default: 5k]`.
  - `[sample].PhiX-removed_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria after PhiX reads have been removed `[Default: 25M]`.
  - `[sample].Adapters_FastA.tsv`: Details if the input adapters reference file meets the minimum file size criteria `[Default: 10k]`.
  - `[sample].Adapter-removed_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria after adapter sequences have been removed `[Default: 25M]`.
  - `[sample].Non-overlapping_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria after removing overlapping reads `[Default: 20M]`.
  - `[sample].Raw_Assembly_File.tsv`: Details if the genome assembly file produced by an assembler software package meets the minimum file size criteria `[Default: 1M]`.
  - `[sample].Filtered_Assembly_File.tsv`: Details if the genome assembly file meets the minimum file size criteria after low compositional complexity contigs are discarded `[Default: 1M]`.
  - `[sample].Binary_PE_Alignment_Map_File.tsv`: Details if the binary paired-end (PE) alignment file meets the minimum file size criteria after the cleaned paired-end reads are mapped onto the filtered genome assembly `[Default: 25M]`.
  - `[sample].Polished_Assembly_File.tsv`: Details if the genome assembly file meets the minimum file size criteria after SNP and InDel corrections are performed `[Default: 1M]`.
  - `[sample].Final_Corrected_Assembly_FastA_File.tsv`: Details if the final error-corrected genome assembly file meets the minimum file size criteria `[Default: 1M]`.
  - `[sample].Binary_SE_Alignment_Map_File.tsv`: Details if the single-end (SE) alignment file meets the minimum file size criteria after the cleaned singleton reads are mapped onto the final genome assembly file `[Default: 1k]`.
  - `[sample].Annotated_GenBank_File.tsv`: Details if the annotated GenBank file meets the minimum file size criteria `[Default: 3M]`.
  - `[sample].SSU_Extracted_File.tsv`: Details if the extracted 16S rRNA gene sequence file meets the minimum file size criteria `[Default: 500b]`.
  - `[sample]-[assembler].SSU_Renamed_File.tsv`: Details if the 16S rRNA gene sequence file meets the minimum file size criteria after sample identifiers are added to each sequence `[Default: 500b]`.
  - `[sample].16S_BLASTn_Output_File.tsv`: Details if the BLASTn output file meets the minimum file size criteria `[Default: 10b]`.
  - `[sample].Filtered_16S_BLASTn_File.tsv`: Details if the best BLASTn alignment sequence meets the minimum file size criteria `[Default: 10b]`.

</details>
