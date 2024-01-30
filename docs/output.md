# wf-paired-end-illumina-assembly: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Input quality control](#input-quality-control) - trimming and contaminant removal
- [Taxonomic classification of trimmed reads](#taxonomic-classification-of-trimmed-reads)
- [Assembly](#assembly) of trimmed reads
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

> *Note: `[sample]` is a unique identifier that is parsed from input FastQ filenames and excludes everything after [R1/R2].*

> *Note: `[assembler]` is the name of the assembler used to assemble contigs. [Default: SPAdes].*


## Input quality control

Input files must meet a size criteria to be processed within this pipeline. If this check passes, the input files go through host removal, downsampling, PhiX read removal, and adapter trimming.

### Initial FastQ file check

***QC Step***: Input files are checked to ensure that they meet a minimum file size to be processed within this pipeline. This is to prevent unusually small input sets from wasting compute time processing data that will not yield a usable bacterial genome assembly.

### Host read removal

Host read removal can be skipped or performed by Hostile and/or NCBI SRA-Human-Scrubber by specifying `--host_remove {both,hostile,sra-human-scrubber,skip}`. For SRA-Human-Scrubber, reads are repaired using BBTools to discard broken sister reads. Information about the number of reads discarded and retained are saved in the output directory.

***QC Step***: FastQ files after host removal are checked to ensure that they meet a minimum file size before continuing on to downstream processes.

<details markdown="1">
<summary>Output files</summary>

- `CleanedReads/Hostile`
  - `[sample].Summary.Hostile-Removal.tsv`: Summary of the number of reads discarded and retained from Hostile.

- `CleanedReads/SRA-Human-Scrubber`
  - `[sample].Summary.SRA-Human-Scrubber-Removal.tsv`: Summary of the number of reads discarded and retained from SRA-Human-Scrubber.

</details>

### Remove PhiX reads

PhiX reads are commonly used as a positive control for Illumina sequencing. PhiX reads must be removed to form a higher quality assembly. A default [PhiX reference file](../bin/PhiX_NC_001422.1.fasta) is used within this pipeline and can be found on [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1/). These PhiX sequences are removed using BBDuk.

***QC Step***: The PhiX reference genome that is used with BBDuk to remove sequence reads must be accessible and meet a minimum file size `[Default: 5k]`.

***QC Step***: FastQ files after PhiX removal are checked to ensure that they meet a minimum file size before continuing on to downstream processes `[Default: 25M]`. This is to halt analysis of a sample mostly containing PhiX reads rather than the sample DNA itself.

<details markdown="1">
<summary>Output files</summary>

- `CleanedReads/BBDUK`
  - `[sample].Summary.PhiX.tsv`: Number of reads discarded and retained from BBDuk.

- `Summaries/`
  - `Summary.PhiX.tab`: Summary of reads discarded and retained for each sample.

</details>

### Adapter clipping and quality trimming

Illumina instruments are able to detect and remove adapter sequences, but sometimes adapters can end up in the FastQ output due to sequencing errors. A default [adapters reference file](../bin/adapters_Nextera_NEB_TruSeq_NuGEN_ThruPLEX.fas) is used within this pipeline and is a concatenation of [Trimmomatic's adapters](https://github.com/usadellab/Trimmomatic/tree/main/adapters) and Illumina adapter sequences from [NCBI's UniVec database](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/). Trimmomatic also performs quality trimming via a sliding window where a sequence must have a specified mean Phred score of 30 for it to be retained. If the Phred score of the first or last basepair of the sequence is below 10, it is removed. Broken sister reads are retained for downstream processes.

***QC Step***: The adapters reference file that is used with Trimmomatic to remove sequence reads must be accessible and meet a minimum file size `[Default: 10k]`.

***QC Step***: FastQ files after removing adapter sequences are checked to ensure that they meet a minimum file size before continuing on to downstream processes `[Default: 25M]`. This is to halt analysis of a sample that is primarily contaminated with artificial Illumina sequences.

<details markdown="1">
<summary>Output files</summary>

- `CleanedReads/Trimmomatic`
  - `[sample].trimmomatic.tsv`: Summary of the number of reads discarded and retained from Trimmomatic.

</details>

### Merge overlapping sister reads

Overlapping content between sister reads that are at least 80% similar are collapsed into a singleton read.

***QC Step***: After singletons have been formed from overlapping read pairs, the remaining non-overlapping paired-end (R1,R2) read files must meet a minimum file size `[Default: 20M]`. This is to prevent the analysis of paired-end read files that are heavily misconstructed into small fragment sizes.

<details markdown="1">
<summary>Output files</summary>

- `CleanedReads`
  - `[sample]_single.fq.gz`: Final cleaned singleton reads.
  - `[sample]_R[1/2].paired.fq.gz`: Final cleaned paired reads.

- `CleanedReads/FLASH`
  - `[sample].overlap.tsv`: Number of reads that were overlapped into singleton reads.
  - `[sample].clean-reads.tsv`: Number of non-overlapping reads.

</details>

## Taxonomic classification of trimmed reads

> *Note: Taxonomic classification tools will be skipped if a database is not specified.*

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

The cleaned and trimmed reads are used to assemble contigs using SPAdes or SKESA `[Default: SPAdes]`. Contigs that have low compositional complexity are discarded. Contigs from SPAdes require polishing to create the final genome assembly, which is done by using BWA, Pilon, and Samtools.

***QC Step***: The contigs produced from an assembler software package must meet a minimum file size criteria `[Default: 1M]`. This is to prevent the analysis of a highly incomplete bacterial genome.

***QC Step***: The resulting contig file after low compositional complexity contigs are discarded must meet a minimum file size `[Default: 1M]`. This is to prevent the analysis of a highly incomplete bacterial genome.

***QC Step***: The cleaned paired-end reads are mapped onto the filtered assembly file in sequential steps (`[Default: 3]`), and the resulting binary paired-end alignment file must meet a minimum file size criteria `[Default: 25M]`. This is to prevent the analysis of an assembly file that has an unusally low read sequence amount.

***QC Step***: The assembly file goes through SNP and InDel corrections in sequential steps (`[Default: 3]`), and the resulting assembly file must meet a minimum file size criteria `[Default: 1M]`. This is to prevent further analysis of an unusually incomplete genome.

***QC Step***: The final error-corrected assembly file must meet a minimum file size criteria `[Default: 1M]`. This is to ensure that the final assembly file is not unexpectedly small or incomplete.

***QC Step***: If singletons (single-end reads) exist after read cleaning, they are mapped onto the assembly file and the resulting binary single-end alignment file must meet a minimum file size criteria `[Default: 1k]`. This is to ensure that read depth calculations can be performed on the single-end reads.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/`
  - `[sample]-[assembler]_contigs.fna`: Final genome assembly.

</details>

### SPAdes

SPAdes is a k-mer based software that forms a genome assembly from read sequences.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/SPAdes/[sample]/`
  - `[sample]-SPAdes.log.gz`: SPAdes log file.
  - `[sample]-SPAdes_graph.gfa`: Assembly graph in gfa format.
  - `[sample]-SPAdes_warnings.log`: Log file that lists warnings when forming a genome assembly.
  - `[sample]-SPAdes_params.txt.gz`: Command used to perform the SPAdes analysis.
  - `[sample]-SPAdes_contigs.fasta`: Assembled contigs in fasta format.
  - `[sample]-SPAdes.SNPs-corrected.cnt.txt`: Number of SNPs corrected in each round of corrections.
  - `[sample]-SPAdes.InDels-corrected.cnt.txt`: Number of InDels corrected in each round of corrections.

</details>

### SKESA

Strategic K-mer Extension for Scrupulous Assemblies (SKESA) is a software that is based on DeBruijn graphs that forms a genome assembly from read sequences.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/SKESA/[sample]/`
  - `[sample]-SKESA_contigs.fasta`: Assembled contigs in fasta format.

</details>

## Assembly information

### Quality assessment and coverage

QUAST is used to perform quality assessment on the assembly file to report metrics such as N50, cumulative length, longest contig length, and GC composition. Using bedtools, the BAM file created with the assembly is used to calculate genome size and the coverage of paired-end and single-end reads.

<details markdown="1">
<summary>Output files</summary>

- `Assembly/QA/[sample]/`
  - `[sample]-[assembler].QuastSummary.tsv`: Assembly metrics such as N50, cumulative length, longest contig length, and GC composition.
  - `[sample]-[assembler].GenomeCoverage.tsv`: Genome coverage information.
  - `[sample]-[assembler].CleanedReads-Bases.tsv`: Number of cleaned bases.

- `Summaries/`
  - `Summary.Assemblies.tab`: Assembly metrics such as N50, cumulative length, longest contig length, and GC composition for each sample.
  - `Summary.GenomeCoverage.tab`: Summary of the overall genome coverage for each sample.
  - `Summary.CleanedReads-Bases.tab`: Summary of the number of cleaned bases for each sample.
  - `Summary.CleanedReads-AlignmentStats.tab`: Summary of the genome size and coverages of the paired-end and single-end reads for each sample.


</details>

### Multilocus sequence typing (MLST)

The final assembly file is scanned against PubMLST typing schemes to determine the MLST for each sample.

<details markdown="1">
<summary>Output files</summary>

- `Summaries/`
  - `Summary.MLST.tab`: Summary of the MLST results for all samples.

</details>

### Genome annotation

The final assembly file is annotated to identify and label features using Prokka.

***QC Step***: The resulting annotated genbank file must meet a minimum file size criteria `[Default: 3M]`. This is to prevent further analysis of a highly incomplete annotation set.

<details markdown="1">
<summary>Output files</summary>

- `Annotation/Prokka`
  - `[sample]-[assembler].gbk`: Annotated genome in genbank file format.

</details>

### 16S rRNA

The genbank file is parsed for 16S rRNA gene records. If there are no 16S rRNA gene records, Barrnap is used to predict 16S rRNA genes using the assembly file. BLAST is then used to align these gene records to it's database, where the best alignment is filtered out based on bitscore.

***QC Step***: The extracted 16S rRNA gene sequence file from the genome assembly must meet a minimum file size criteria `[Default: 500b]`. This is to prevent the classification of a highly incomplete 16S rRNA gene sequence.

***QC Step***: The sample identifiers in the 16S rRNA gene sequence file get renamed and this resulting file must meet a minimum file size critiera `[Default: 500b]`. This is to prevent the classification of a highly incomplete 16S rRNA gene sequence.

***QC Step***: The 16S rRNA gene sequnce file is classified using BLASTn and the resulting output file must meet a minimum file size criteria `[Default: 10b]`. This is to ensure that the BLASTn output file contains at least one alignment with taxonomic information to be filtered and reported.

***QC Step***: The best 16S BLASTn alignment sequence filtered by bitscore is saved to a file in FastA format and must meet a minimum file size criteria `[Default: 10b]`. This is to ensure that the resulting file contains alignment information that can be parsed and reported.

<details markdown="1">
<summary>Output files</summary>

- `Summaries/`
  - `Summary.16S.tab`: Summary of the best BLAST alignment for each sample.

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


## Pipeline information

Information about the pipeline execution is stored here, along with output logs, error logs, and QC file checks for each sample.

<details markdown="1">
<summary>Pipeline information</summary>

- `pipeline_info/`
  - `software_versions.yml`: Summary of the software packages used in each process and their version information.
  - `pipeline_dag_YYYY-MM-DD_HH-MM-SS.html`: Direct acrylic graph (DAG) image of the workflow that gives a visual representation of how each process connects to each other.
  - `execution_trace_YYYY-MM-DD_HH-MM-SS.txt`: Text-based summary report detailing the work directory hash, runtime, CPU usage, memory usage, etc. for each process.
  - `execution_report_YYYY-MM-DD_HH-MM-SS.html`: Summary report of all processes, including processes that passed/failed, resource usage, etc.
  - `execution_timeline_YYYY-MM-DD_HH-MM-SS.html`: Summary report detailing the runtime and memory usage of each process.

</details>

<details markdown="1">
<summary>Process log information</summary>

- `pipeline_info/process_logs/`
  - `[sample].[process].command.out`: Output log file for each sample in a given process.
  - `[sample].[process].command.err`: Error log file for each sample in a given process.

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
  - `[sample].Annotated_GenBank_File.tsv`: Details if the annotated genbank file meets the minimum file size criteria `[Default: 3M]`.
  - `[sample].SSU_Extracted_File.tsv`: Details if the extracted 16S rRNA gene sequence file meets the minimum file size critieria `[Default: 500b]`.
  - `[sample]-[assembler].SSU_Renamed_File.tsv`: Details if the 16S rRNA gene sequence file meets the minimum file size criteria after sample identifiers are added to each sequence `[Default: 500b]`.
  - `[sample].16S_BLASTn_Output_File.tsv`: Details if the BLASTn output file meets the minimum file size criteria `[Default: 10b]`.
  - `[sample].Filtered_16S_BLASTn_File.tsv`: Details if the best BLASTn alignment sequence meets the minimum file size criteria `[Default: 10b]`.

</details>
