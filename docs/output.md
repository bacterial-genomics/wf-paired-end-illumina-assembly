# wf-paired-end-illumina-assembly: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and is used to perform _de novo_ assembly on raw Illumina paired-end reads from bacterial isolates.

- [Input quality control](#input-quality-control) that includes trimming and contaminant removal
  - [Initial FastQ file check](#initial-fastq-file-check) ensures input FastQ files meet a minimum file size
  - [Host read removal](#host-read-removal)
  - [PhiX read removal](#phix-read-removal)
  - [Adapter clipping and quality trimming](#adapter-clipping-and-quality-trimming)
  - [Merge overlapping sister reads](#merge-overlapping-sister-reads) to create singletons
- [Taxonomic classification of cleaned reads](#taxonomic-classification-of-cleaned-reads)
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
>
> `[sample]` is a unique identifier that is parsed from input FastQ filenames and excludes everything after [R1/R2].
>
> `[assembler]` is the name of the assembler used to assemble contigs. [Default: SPAdes].

> [!TIP]
> All tab-separated value (TSV) files can be converted to Excel spreadsheets (XLSX) by using the parameter `--create_excel_outputs` when running the pipeline.
>
> When using this parameter, a summary workbook is created to allow for all summary files to be added to separate worksheets within the workbook.

## Input quality control

Input files are checked for corruption and must meet a minimum file size to be processed within this pipeline. If this check passes, the input files will go through several read cleaning steps before other analyses are performed.

### Initial FastQ file check

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- Input files are checked to ensure that they meet a minimum file size to be processed within this pipeline `[Default: 25M]`. This is to prevent unusually small input sets from wasting compute time processing data that will not yield a usable bacterial genome assembly.

</details>

### Host read removal

Host read removal can be skipped or performed by Hostile and/or NCBI SRA-Human-Scrubber by specifying `--host_remove {both,hostile,sra-human-scrubber,skip}`. When `both` is invoked, they occur in sequential fashion-- first SRA Scrub tool and repair, then hostile occurs. For SRA-Human-Scrubber, reads are repaired using BBTools to discard broken sister reads. Information about the number of reads discarded and retained are saved in the output directory.
Please see the [host removal using Hostile documentation](../modules/local/remove_host_hostile/README.md) and [host removal using SRA-Human-Scrubber documentation](../modules/local/remove_host_sra_human_scrubber/README.md) for more information.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- FastQ files after host removal are checked to ensure that they meet a minimum file size before continuing to downstream processes `[Default: 25M]`.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `Clean_Reads/Hostile/`

  - `[sample].fastq.gz`: Human host-removed FastQ files.
  - `[sample].Summary.Hostile.tsv`: Summary of the number of reads discarded and retained from Hostile.
  - `[sample].Hostile_FastQ.SHA512-checksums.tsv`: Checksum values for each FastQ output from Hostile.

- `Clean_Reads/SRA-Human-Scrubber/`
  - `[sample].Summary.BBTools_Repair_Removal.tsv`: Summary of the number of reads discarded and retained after repairing broken sister reads produced from SRA-Human-Scrubber.
  - `[sample].Summary.SRA_Human_Scrubber_Removal.tsv`: Summary of the number of reads discarded and retained from SRA-Human-Scrubber.

</details>

### PhiX read removal

PhiX reads are commonly used as a positive control for Illumina sequencing. During assembly, PhiX reads are considered contaminants and if retained, a misassembled genome will be formed. Therefore, a PhiX reference file is required and a default [PhiX reference file](../bin/PhiX_NC_001422.1.fasta) is included with this pipeline. Please see the [PhiX removal using BBDuk documentation](../modules/local/remove_phix_bbduk/README.md) for more information.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- The PhiX reference genome that is used with BBDuk to remove sequence reads must be accessible and meet a minimum file size `[Default: 5k]`.

- FastQ files after PhiX removal are checked to ensure that they meet a minimum file size before continuing to downstream processes `[Default: 25M]`. This is to halt analysis of a sample mostly containing PhiX reads rather than the sample DNA itself.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `Clean_Reads/BBDuk/`
  - `[sample].noPhiX_FastQ.SHA512-checksums.tsv`: Checksum values for each PhiX-free FastQ output of BBDuk.
  - `[sample].PhiX_Removed_Reads.metrics_summary.tsv`: Metrics on FastQ reads output including minimum, average, and maximum lengths as well as total counts and total length.
  - `[sample].PhiX.tsv`: Number of reads discarded and retained from BBDuk.

</details>

### Adapter clipping and quality trimming

Illumina instruments can detect and remove adapter sequences, but sometimes adapters can end up in the FastQ output due to sequencing errors.

#### Trimmomatic

A default [adapters reference file](../bin/adapters_Nextera_NEB_TruSeq_NuGEN_ThruPLEX.fas) is used within this pipeline for Trimmomatic. Trimmomatic also performs quality trimming, where broken sister reads are retained for downstream processes.
Please see the [adapter clipping and quality trimming using Trimmomatic documentation](../modules/local/trim_reads_trimmomatic/README.md) for more information.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- The adapters reference file that is used with Trimmomatic to remove sequence reads must be accessible and meet a minimum file size `[Default: 10k]`.

- FastQ files after removing adapters and performing quality trimming are checked to ensure that they meet a minimum file size before continuing to downstream processes `[Default: 25M]`. This is to halt analysis of a sample that is primarily contaminated with artificial Illumina sequences.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `Clean_Reads/Trimmomatic/`
  - `[sample].Adapter_and_Quality_Trimmed_Reads.metrics_summary.tsv`: Metrics on FastQ reads output including minimum, average, and maximum lengths as well as total counts and total length.
  - `[sample].Trim_FastQ.SHA512-checksums.tsv`: Checksum values for each FastQ output from Trimmomatic.
  - `[sample].Trimmomatic.tsv`: Summary of the number of reads discarded and retained from Trimmomatic.

</details>

#### fastp

fastp is able to clip adapters, perform quality trimming, and retain broken sister reads for downstream processes. fastp is able to automatically detect adapter sequences, but the use of a custom list of adapter sequences can be used.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- If used, the adapters reference file that is used to remove sequence reads must be accessible and meet a minimum file size `[Default: 10k]`.

- FastQ files after removing adapters and performing quality trimming are checked to ensure that they meet a minimum file size before continuing to downstream processes `[Default: 25M]`. This is to halt analysis of a sample that is primarily contaminated with artificial Illumina sequences.

</details>

<br />

<details markdown="1">
<summary>Output files</summary>

- `Clean_Reads/fastp/`
  - `[sample].Adapter_and_Quality_Trimmed_Reads.metrics_summary.tsv`: Metrics on FastQ reads output including minimum, average, and maximum lengths as well as total counts and total length.
  - `[sample].Trim_FastQ.SHA512-checksums.tsv`: Checksum values for each FastQ output from Fastp.
  - `[sample].Fastp.tsv`: Summary of the number of reads discarded and retained from Fastp.

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

- `Clean_Reads/`

  - `[sample]_single.fq.gz`: Final cleaned singleton reads.
  - `[sample]_R[1/2].paired.fq.gz`: Final cleaned paired reads.

- `Clean_Reads/FLASh/`
  - `[sample].Clean_Reads_FastQ.metrics_summary.tsv`: Metrics on FastQ reads output including minimum, average, and maximum lengths as well as total counts and total length.
  - `[sample].Clean_Reads_FastQ.SHA512-checksums.tsv`: Checksum values for each FastQ output from FLASh.
  - `[sample].Overlap.tsv`: Number of reads that were overlapped into singleton reads.

> [!NOTE]
> FastQ sequences after overlapping with FLASh are stored in `Clean_Reads/` rather than `Clean_Reads/FLASh` as they are the final outputs for read cleaning.

</details>

## Taxonomic classification of cleaned reads

These classifiers perform classifications on a read-by-read basis or through the use of k-mers on the cleaned FastQ files. The results that are obtained are heavily dependent on the quality and diversity of the database used. Therefore, the results produced from these classifiers should be used as a quality check to evaluate the possibility of sample contamination.

> [!WARNING]
> Taxonomic classification tools will be skipped if the accompanying database is not specified.

### Kraken

Kraken is a k-mer based classification tool that assigns taxonomic labels using the Lowest Common Ancestor (LCA) algorithm. It consumes much more RAM than Kraken2, but Kraken can be valuable for precise exact k-mer matches of target organisms with lots of highly similar near taxonomic neighbors (e.g., occurring within a species complex).

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/Kraken/[sample]/`
  - `[sample].kraken_summary.tsv`: Summary of the unclassified, top 3 genus and top 3 species classifications from the Kraken report.
  - `[sample].kraken_output.tsv.gz`: Full taxonomic k-mer matches, not yet filtered for top hits, in the Kraken report format.

</details>

### Kraken2

Kraken2 is a minimizer based classification tool that assigns taxonomic labels using the Lowest Common Ancestor (LCA) algorithm.

<details markdown="1">
<summary>Output files</summary>

- `Taxonomy/Kraken2/[sample]/`
  - `[sample].kraken2_summary.tsv`: Summary of the unclassified, top 3 genus and top 3 species classifications from the Kraken2 report.
  - `[sample].kraken2_output.tsv.gz`: Full taxonomic minimizer matches, not yet filtered for top hits, in the Kraken2 report format.

> [!TIP]
> Some pure isolates consistently give near neighbor matches, but you might wonder if an assembly used in database creation was contaminated. If you're curious which specific assemblies were used in your Kraken or Kraken2 database creation, you can view each and every one of them!
>
> Many of the pre-computed databases from Dr. Ben Langmead [here](https://benlangmead.github.io/aws-indexes/k2) have a file called "inspect.txt" which has a line-by-line listing of each assembly used to create the database. If you do not have that, you can create it with `kraken-inspect --db /my/db/path --threads 12 > inspect.txt` or `kraken2-inspect --db /my/db/path --threads 12 > inspect.txt` depending on whether you're using Kraken or Kraken2 (respectively).

</details>

## Assembly

The cleaned reads are used to assemble contigs using SPAdes or SKESA `[Default: SPAdes]`. Contigs that have low compositional complexity are discarded. Please see the [contig filtering documentation](../modules/local/filter_contigs_biopython/README.md) for more information. Contigs from SPAdes also involves SNP and InDel correction/polishing to create the final genome assembly, which is done by using BWA, Pilon, and Samtools. Contigs from SKESA do not require this step.

> [!IMPORTANT]
> Outputs generated by SPAdes and SKESA cannot be compared even when using the same FastQ inputs.

> [!TIP]
> For many input FastQ files where you might want to decrease overall runtime, SKESA is useful to still obtain high SNP-level accuracy but at the expense of longer contigs. For tasks where contiguity (longer contigs) are important (e.g., gene neighborhood evaluations, operon detection), SPAdes is the more appropriate choice.

<details markdown="1">
<summary><strong>QC Steps</strong></summary>

- The contigs produced from an assembler software package must meet a minimum file size criteria `[Default: 1M]`. This is to prevent the analysis of a highly incomplete bacterial genome.

- The resulting contig file (after filtering out low coverage, short, and low compositional complexity contigs) must meet a minimum file size `[Default: 1M]`. This is to prevent the analysis of a highly incomplete bacterial genome.

- The cleaned paired-end reads are mapped onto the filtered assembly file in sequential steps (`[Default: 3]`), and the resulting binary paired-end alignment file must meet a minimum file size criteria `[Default: 6M]`. This is to prevent the analysis of an assembly file that has an unusually low read sequence amount.

- The assembly file goes through SNP and InDel corrections in sequential steps (`[Default: 3]`), and the resulting assembly file must meet a minimum file size criteria `[Default: 1M]`. This is to prevent further analysis of an unusually incomplete genome.

- The final error-corrected assembly file must meet a minimum file size criteria `[Default: 1M]`. This is to ensure that the final assembly file is not unexpectedly small or incomplete.

- If singletons (single-end reads) exist after read cleaning, they are mapped onto the assembly file and the resulting binary single-end alignment file must meet a minimum file size criteria `[Default: 1k]`. This is to ensure that read depth calculations can be performed on the single-end reads.

> [!TIP]
> Discarded contigs from filtering are stored in `Assembly/[assembler]/[sample]/[sample]-[assembler].discarded-contigs.fa.gz`. You can view the reason for each individual contig being discarded by `zcat discarded-contigs.fa.gz | grep '>'` and within the contig name there will be 1 or more reasons listed. For example "Failed=complexityFailed=lengthFailed=gc_content" had 3 independent reasons for being removed, whereas "Failed=length" was simply too short of a contig.

> [!TIP]
> Contig filtering statistics are stored in `Assembly/[assembler]/[sample]/[sample]-[assembler].filter-contigs-stats.txt`. There you'll find total contig counts and cumulative lengths for input, removed, and saved. Also, for coverage statistics there are minimum, average, maximum, 25% quartile, 50% quartile (median), and 75% quartile coverage values. All of these statistics are meant to help guide alternative contig filtering if you have an unusual assembly that requires non-default parameters.

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
  - `[sample]-SPAdes_contigs.fasta`: Assembled contigs in FastA format.
  - `[sample]-SPAdes_graph.gfa`: Assembly graph in gfa format.
  - `[sample]-SPAdes.discarded-contigs.fa.gz`: Post-assembly contigs that were filtered out (i.e., discarded).
  - `[sample]-SPAdes.filter-contigs-stats.txt`: Post-assembly contig filtering statistics.
  - `[sample]-SPAdes.InDels-corrected.cnt.txt`: Number of InDels corrected in each round of corrections.
  - `[sample]-SPAdes.log.gz`: SPAdes log file.
  - `[sample]-SPAdes_params.txt.gz`: Command used to perform the SPAdes analysis.
  - `[sample]-SPAdes_warnings.log`: Log file that lists warnings when forming a genome assembly.
  - `[sample]-SPAdes.SNPs-corrected.cnt.txt`: Number of SNPs corrected in each round of corrections.

</details>

### SKESA

Strategic K-mer Extension for Scrupulous Assemblies (SKESA) is a software that is based on DeBruijn graphs that forms a genome assembly from read sequences.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/SKESA/[sample]/`
  - `[sample]-SKESA_contigs.fasta`: Assembled contigs in FastA format.
  - `[sample]-SKESA.discarded-contigs.fa.gz`: Post-assembly contigs that were filtered out (i.e., discarded).
  - `[sample]-SKESA.filter-contigs-stats.txt`: Post-assembly contig filtering statistics.

</details>

## Assembly metrics and classification

### Quality assessment and coverage

QUAST is used to perform quality assessment on the assembly file to report metrics such as N50, cumulative length, longest contig length, and GC composition. Using bedtools, the BAM file created with the assembly is used to calculate genome size and the coverage of paired-end and single-end reads.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/QA/[sample]/`
  - `[sample]-[assembler].QuastSummary.tsv`: Assembly metrics such as N50, cumulative length, longest contig length, and GC composition.
  - `[sample]-[assembler].GenomeCoverage.tsv`: Genome coverage information.
  - `[sample]-[assembler].Clean_Reads-Bases.tsv`: Number of cleaned bases.

</details>

### Multilocus sequence typing (MLST)

The final assembly file is scanned against PubMLST typing schemes to determine the MLST for each sample. Unless a specific typing scheme is specified, the best typing scheme for each sample is automatically selected.

<details markdown="1">
<summary>Output files</summary>

- `Summaries/`
  - `Summary.MLST.tsv`: Summary of the MLST results for all samples.

</details>

<br />

<details markdown="1">
<summary>MLST output interpretation</summary>

| Symbol | Meaning                               | Length          | Identity       |
|--------|---------------------------------------|-----------------|----------------|
| `n`    | exact intact allele                   | 100%            | 100%           |
| `~n`   | novel full length allele similar to n | 100%            | &ge; `--minid` |
| `n?`   | partial match to known allele         | &ge; `--mincov` | &ge; `--minid` |
| `-`    | allele missing                        | &lt; `--mincov` | &lt; `--minid` |
| `n,m`  | multiple alleles                      |                 |                |

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

The GenBank file is parsed for 16S rRNA gene records (with BioPython). If there are no 16S rRNA gene records, Barrnap, with relaxed settings, is used to find partial 16S rRNA genes using the assembly file. BLAST+ (blastn) is then used to align these gene records to its database, where the best alignment is filtered out based on bit score. The default database is an NCBI currated set of 16S rRNA genes of species Type strains, but best matches are not always perfect.

> [!NOTE]
> Some assembled genomes do not contain identifiable 16S rRNA sequences and therefore 16S is not able to be classified. If the classification of 16S rRNA sequences is required, the sample must be re-sequenced.

> [!IMPORTANT]
> The 16S rRNA classification produced should **not** be used as a definitive classification as some taxa have 16S sequences that are extremely similar between different species.
>
> When the top bitscore match is < 100% identity or < 100% alignment, you should be extra cautious in what species you're reporting.

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

  - `16S.[sample]-[assembler].fa`: 16S rRNA gene sequence of the best BLAST alignment in FastA format.
  - `[assembler].16S_top_genus_RDP.tsv`: Summary of the best RDP match for each sample.
  - `[assembler].16S_top_species_BLAST.tsv`: Summary of the best BLAST alignment for each sample.

- `SSU/BLAST/`
  - `[sample]-[assembler].blast.tsv.gz`: Full, not yet bitscore sorted, BLASTn output for each 16S rRNA gene record in tab-separated value (TSV) format using the BLAST outfmt 6 standard with additional taxonomy fields

- `SSU/RDP/`
  - `[sample]-[assembler].rdp.tsv`: RDP classification output for each 16S rRNA gene record in tab-separated value (TSV) format

</details>

### Assembly taxonomic classification

GTDB-Tk is a taxonomic classification tool that uses the Genome Database Taxonomy [GTDB](https://gtdb.ecogenomic.org/). GTDB-Tk can be used as a quality check for the assembly file.

<details markdown="1">
<summary>Output files</summary>

- `Summaries/`
  - `Summary.GTDB-Tk.tsv`: Summary of the GTDB-Tk taxonomic classification for each sample.

</details>

## Summaries

Concatenation of output metrics for all samples.

> [!NOTE]
> The first column for **all** "Summary" files contains the sample name.
>
> Each Summary file is sorted based on sample names for easy cross-comparisons.

<details markdown="1">
<summary>Output files</summary>

- `Summaries/`
  - `Summary.16S_Genus_RDP.tsv`: RDP Classifier best matches from predicted 16S ribosomal RNA genes (to genus-level)
  - `Summary.16S_Species_BLAST.tsv`: Top bitscore BLASTn alignments for each 16S rRNA gene (to species-level)
  - `Summary.Adapter_QC_Trim_Reads.Metrics.tsv`: Sequence metrics after adapter clipping and quality trimming
  - `Summary.Annotation_Checksums.tsv`: Checksum (hash) values for annotated GenBank output files
  - `Summary.Assembly_Checksums.tsv`: Checksum (hash) values for final output assembly FastA output files
  - `Summary.Assembly_Depth.tsv`: Assembly depth of coverage mean and standard deviation values (units in "x")
  - `Summary.Assembly_Metrics.tsv`: Assembly metrics (e.g., N50, cumulative length, longest contig length, and GC composition)
  - `Summary.CheckM2.tsv`: Estimation percentages on completeness and contamination of each genome assembly
  - `Summary.Clean_and_Overlapped.tsv`: Counts and percentages from overlapping sister reads
  - `Summary.Clean_Reads_Aligned.tsv`: Paired, singleton, and total cleaned reads mapped onto the assembly statistics (e.g., mean and standard deviation depths, basepairs mapped, assembly size)
  - `Summary.Clean_Reads_Checksums.tsv`: Checksum (hash) values for final output cleaned reads FastQ output files
  - `Summary.Clean_Reads.Metrics.tsv`: Sequence metrics after all read cleaning steps
  - `Summary.Downsampled_Reads.Metrics.tsv`: Sequence metrics after subsampling the read set (if performed)
  - `Summary-[Fastp,Trimmomatic].Adapter_and_QC_Trim.tsv`: Number of discarded reads and singleton reads that remain after adapter clipping and quality trimming
  - `Summary.Input_Checksums.tsv`: Checksum (hash) values for FastQ input files
  - `Summary.Input_Reads.Metrics.tsv`: Sequence metrics on the initial user-provided input sequences
  - `Summary.Kraken2.tsv`: Counts and proportions of unclassified, top 3 genera, top 3 species with k-mer matches with the cleaned reads
  - `Summary.Kraken.tsv`: Counts and proportions of unclassified, top 3 genera, top 3 species with k-mer matches with the cleaned reads
  - `Summary.MLST.tsv`: MLST genotyping results
  - `Summary.PhiX_Removal.tsv`: Number of reads discarded and retained after PhiX k-mer match removal
  - `Summary.PhiX_Removed_Reads.Metrics.tsv`: Sequence metrics after PhiX was removed from the reads
  - `Summary.QC_File_Checks.tsv`: All QC file checks detailing if a sample passes or fails after each process
  - `Summary-Report.xlsx`: Microsoft Excel workbook where each file in the Summaries directory is added to a separate worksheet within the workbook

</details>

## Pipeline information

Information about the pipeline execution, output logs, error logs, and QC file checks for each sample are stored here.

> [!NOTE]
> Pipeline execution files have the date and time appended to the filename using the following shorthand notation: year (yyyy), month (MM), day (dd), hour (HH), minute (mm), second (ss).

<details markdown="1">
<summary>Pipeline information</summary>

- `pipeline_info/`
  - `software_versions.yml`: All software packages used and their version information
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
  - `[sample].Summary.Hostile.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria for after host removal using Hostile `[Default: 25M]`.
  - `[sample].SRA_Human_Scrubber_FastQ_File.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria for after host removal using SRA-Human-Scrubber `[Default: 25M]`.
  - `[sample].BBTools-Repair-removed_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria after repairing broken sister reads from SRA-Human-Scrubber `[Default: 25M]`.
  - `[sample].PhiX_Genome.tsv`: Details if the input PhiX reference genome meets the minimum file size criteria `[Default: 5k]`.
  - `[sample].PhiX-removed_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria after PhiX reads have been removed `[Default: 25M]`.
  - `[sample].Adapters_FastA.tsv`: Details if the input adapters reference file meets the minimum file size criteria `[Default: 10k]`.
  - `[sample].Adapter-removed_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria after adapter sequences have been removed `[Default: 25M]`.
  - `[sample].Non-overlapping_FastQ_Files.tsv`: Details if both reads (R1,R2) meet the minimum file size criteria after removing overlapping reads `[Default: 20M]`.
  - `[sample].Raw_Assembly_File.tsv`: Details if the genome assembly file produced by an assembler software package meets the minimum file size criteria `[Default: 1M]`.
  - `[sample].Filtered_Assembly_File.tsv`: Details if the genome assembly file meets the minimum file size criteria after low compositional complexity contigs are discarded `[Default: 1M]`.
  - `[sample].Binary_PE_Alignment_Map_File.tsv`: Details if the binary paired-end (PE) alignment file meets the minimum file size criteria after the cleaned paired-end reads are mapped onto the filtered genome assembly `[Default: 6M]`.
  - `[sample].Polished_Assembly_File.tsv`: Details if the genome assembly file meets the minimum file size criteria after SNP and InDel corrections are performed `[Default: 1M]`.
  - `[sample].Final_Corrected_Assembly_FastA_File.tsv`: Details if the final error-corrected genome assembly file meets the minimum file size criteria `[Default: 1M]`.
  - `[sample].Binary_SE_Alignment_Map_File.tsv`: Details if the single-end (SE) alignment file meets the minimum file size criteria after the cleaned singleton reads are mapped onto the final genome assembly file `[Default: 1k]`.
  - `[sample].Annotated_GenBank_File.tsv`: Details if the annotated GenBank file meets the minimum file size criteria `[Default: 3M]`.
  - `[sample].SSU_Extracted_File.tsv`: Details if the extracted 16S rRNA gene sequence file meets the minimum file size criteria `[Default: 500b]`.
  - `[sample]-[assembler].SSU_Renamed_File.tsv`: Details if the 16S rRNA gene sequence file meets the minimum file size criteria after sample identifiers are added to each sequence `[Default: 500b]`.
  - `[sample].16S_BLASTn_Output_File.tsv`: Details if the BLASTn output file meets the minimum file size criteria `[Default: 10b]`.
  - `[sample].Filtered_16S_BLASTn_File.tsv`: Details if the best BLASTn alignment sequence meets the minimum file size criteria `[Default: 10b]`.

</details>
