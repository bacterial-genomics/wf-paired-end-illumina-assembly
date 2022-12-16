#!/usr/bin/env nextflow


/*
==============================================================================
                              wf-paired-end-illumina-assembly                              
==============================================================================
usage: nextflow run main.nf [--help]
----------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     wf-paired-end-illumina-assembly v${version}
    =========================================
    Usage:
    nextflow run -profile <docker|singularity> main.nf --inpath <input directory> --outpath <directory for results>

    Run with test data:
    nextflow run main.nf -profile test,<docker|singularity>
    
    Input/output options:
      --inpath             Path to input data directory containing FastQ assemblies. Recognized extensions are:  fastq.gz, fq.gz.
      --outpath            The output directory where the results will be saved.
    Analysis options:
      --size               Specify file size that is used to verify minimum file sizes in all processes in BYTES, e.g. 1000.
      --kraken1_db         Specify path to database for Kraken1. Default database is Mini Kraken.
      --kraken2_db         Specify path to database for Kraken2. Default database is Mini Kraken.
      --blast_db           Specify path to 16S ribosomal database for BLAST. Default database is NCBI's 16S ribosomal database.
      --bigdata            Whether or not to use more compute resources. Options are true, false (default).
      --max_memory         Specify memory limit on your machine/infrastructure, e.g. '128.GB'. Useful to ensure workflow doesn't request too many resources.
      --max_time           Specify time limit for each process, e.g. '240.h'. Useful to ensure workflow doesn't request too many resources.
      --max_cpus           Specify CPU limit on your machine/infrastructure, e.g. 16. Useful to ensure workflow doesn't request too many resources.
    Profile options:
      -profile singularity Use Singularity images to run the workflow. Will pull and convert Docker images from Dockerhub if not locally available.
      -profile docker      Use Docker images to run the workflow. Will pull images from Dockerhub if not locally available.
      -profile conda       TODO: this is not implemented yet.
    Other options:
      -resume              Re-start a workflow using cached results. May not behave as expected with containerization profiles docker or singularity.
      -stub                Use example output files for any process with an uncommented stub block. For debugging/testing purposes.
      -name                Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

version = "1.0.0"
nextflow.enable.dsl=2

if (params.help) {
    helpMessage()
    exit 0
}

if (params.version){
    println "VERSION: $version"
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Validate inpath parameter
if (!params.inpath) {
    System.err.println "ERROR: parameter inpath must be specified"
    exit 1
}
File inpathFileObj = new File(params.inpath)
if (!inpathFileObj.exists()){
    System.err.println "ERROR: $params.inpath doesn't exist"
    exit 1
}

// Validate outpath parameter
File outpathFileObj = new File(params.outpath)
if (outpathFileObj.exists()){
    // Per the config file, outpath stores log & trace files so it is created before this point
    // Check that outpath only contains a trace file created this hour
    dayAndHour = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
    outFiles = outpathFileObj.list()
    if (!(outFiles[0] ==~ /trace.($dayAndHour).txt/ && outFiles.size() == 1)) {
        // If it contains an older trace file or other files, warn the user
        System.out.println "WARNING: $params.outpath already exists. Output files will be overwritten."
    }
} else {
    outpathFileObj.mkdirs()
}

// Set logpath parameter
File logpathFileObj = new File(params.logpath)
if (logpathFileObj.exists()){
    System.out.println "WARNING: $params.logpath already exists. Log files will be overwritten."
} else {
    logpathFileObj.mkdirs()
}

// Set optional params; will be overwritten if user specifies
params.kraken1_db = "Pre-loaded MiniKraken1"
params.kraken2_db = "Pre-loaded MiniKraken2"
params.blast_db = "Pre-loaded 16S rRNA"
params.size = "null"

// Print parameters used
log.info """
    =====================================
    wf-paired-end-illumina-assembly $version
    =====================================
    inpath:             ${params.inpath}
    outpath:            ${params.outpath}
    logpath:            ${params.logpath}
    workDir:            ${workflow.workDir}
    kraken1_db:         ${params.kraken1_db}
    kraken2_db:         ${params.kraken2_db}
    blast_db:           ${params.blast_db}
    =====================================
    """
    .stripIndent()

/*
========================================================================================
                 Import local custom modules and subworkflows                 
========================================================================================
*/

include { INFILE_HANDLING } from "./modules/local/infile_handling.nf"
include { REMOVE_PHIX } from "./modules/local/remove_phix.nf"
include { TRIMMOMATIC } from "./modules/local/trimmomatic.nf"
include { EXTRACT_SINGLETONS } from "./modules/local/extract_singleton.nf"
include { KRAKEN_ONE; KRAKEN_TWO; } from "./modules/local/kraken.nf"
include { SPADES } from "./modules/local/spades.nf"
include { FILTER_CONTIGS } from "./modules/local/filter_contigs.nf"
include { CLEAN_READS } from "./modules/local/clean_reads.nf"
include { CLEANED_COVERAGE } from "./modules/local/cleaned_coverage.nf"
include { MLST } from "./modules/local/mlst.nf"
include { ANNOTATE } from "./modules/local/annotate.nf"
include { EXTRACT_RECORDS } from "./modules/local/extract_records.nf"
include { BARRNAP } from "./modules/local/barrnap.nf"
include { BLAST } from "./modules/local/blast.nf"
include { FILTER_BLAST } from "./modules/local/filter_blast.nf"
include { QA } from "./modules/local/qa.nf"
include { GENOME_COVERAGE } from "./modules/local/genome_coverage.nf"

/*
========================================================================================
                   Import nf-core modules and subworkflows                    
========================================================================================
*/

// None

/*
========================================================================================
                            Run the main workflow                             
========================================================================================
*/

workflow {

    // SETUP: Define input, output
    // Double asterisk looks in specified directory and recursively all subdirectories
    input_ch = Channel.fromFilePairs(params.inpath+'/**R{1,2}*.{fastq,fq}.gz', checkIfExists: true)
    output_ch = Channel.fromPath(params.outpath)
    
    // SETUP: Define optional database inputs
    kraken1_db_ch = file(params.kraken1_db)
    kraken2_db_ch = file(params.kraken2_db)
    blast_db_ch = file(params.blast_db)

    // SETUP: Define empty channels to concatenate certain outputs
    ch_versions = Channel.empty()
    alnstats_summary_ch = Channel.empty()
    blast_summary_ch = Channel.empty()
    ssu_species_ch = Channel.empty()
    genome_cov_summary_ch = Channel.empty()
    mlst_summary_ch = Channel.empty()
    assembly_summary_ch = Channel.empty()
    cleaned_summary_ch = Channel.empty()

    // PROCESS: Read files from input directory, validate and stage input files
    INFILE_HANDLING (
        input_ch
    )

    // Collect version info
    ch_versions = ch_versions.mix(INFILE_HANDLING.out.versions)

    // PROCESS: Run bbduk to remove PhiX reads
    REMOVE_PHIX (
        INFILE_HANDLING.out.input
    )

    // Collect version info
    ch_versions = ch_versions.mix(REMOVE_PHIX.out.versions)

    // PROCESS: Run trimmomatic to clip adapters and do quality trimming
    TRIMMOMATIC (
        REMOVE_PHIX.out.nophix
    )

    // Collect version info
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)

    // PROCESS: Run flash to merge overlapping sister reads into singleton reads
    EXTRACT_SINGLETONS (
        INFILE_HANDLING.out.input.join(TRIMMOMATIC.out.trimmo)
    )

    // Collect version info
    ch_versions = ch_versions.mix(EXTRACT_SINGLETONS.out.versions)

    // PROCESS: Run kraken1 on paired reads
    KRAKEN_ONE (
        EXTRACT_SINGLETONS.out.gzip_reads,
        kraken1_db_ch
    )

    // Collect version info
    ch_versions = ch_versions.mix(KRAKEN_ONE.out.versions)

    // PROCESS: Run kraken2 on paired reads
    KRAKEN_TWO (
        EXTRACT_SINGLETONS.out.gzip_reads,
        kraken2_db_ch
    )

    // Collect version info
    ch_versions = ch_versions.mix(KRAKEN_TWO.out.versions)

    // PROCESS: Run SPAdes to assemble contigs
    SPADES (
        EXTRACT_SINGLETONS.out.gzip_reads
    )

    // Collect version info
    ch_versions = ch_versions.mix(SPADES.out.versions)

    // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
    FILTER_CONTIGS (
        EXTRACT_SINGLETONS.out.gzip_reads.join(SPADES.out.contigs)
    )

    // Collect version info
    ch_versions = ch_versions.mix(FILTER_CONTIGS.out.versions)

    // PROCESS: Use BWA/Samtools/Pilon to correct contigs with PE reads
    CLEAN_READS (
        EXTRACT_SINGLETONS.out.gzip_reads.join(FILTER_CONTIGS.out.uncorrected_contigs)
    )

    // Collect version info
    ch_versions = ch_versions.mix(CLEAN_READS.out.versions)

    // PROCESS: Run Bedtools to calculate coverage
    CLEANED_COVERAGE (
        CLEAN_READS.out.bam
    )

    // Collect all Summary Stats and concatenate into one file
    alnstats_summary_ch = alnstats_summary_ch.mix(CLEANED_COVERAGE.out.summary_alnstats)
    alnstats_summary_ch.collectFile(name: 'Summary.Illumina.CleanedReads-AlnStats.tab', storeDir: "${params.outpath}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(CLEANED_COVERAGE.out.versions)

    // PROCESS: Run MLST to find MLST for each assembly
    MLST (
        CLEAN_READS.out.bam.join(CLEAN_READS.out.base_fna)
    )

    // Collect all MLST Summaries and concatenate into one file
    mlst_summary_ch = mlst_summary_ch.mix(MLST.out.summary_mlst)
    mlst_summary_ch.collectFile(name: 'Summary.MLST.tab', storeDir: "${params.outpath}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(MLST.out.versions)

    // PROCESS: Run Prokka to annotate reads
    ANNOTATE (
        CLEAN_READS.out.bam.join(CLEAN_READS.out.base_fna)
    )

    // Collect version info
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)

    // PROCESS: Extract records from annotation file
    EXTRACT_RECORDS (
        ANNOTATE.out.annotation.join(CLEAN_READS.out.base_fna)
    )

    // Collect version info
    ch_versions = ch_versions.mix(EXTRACT_RECORDS.out.versions)

    // PROCESS: Run Barrnap to predict ribosomal RNA genes
    BARRNAP (
        ANNOTATE.out.annotation.join(CLEAN_READS.out.base_fna).join(EXTRACT_RECORDS.out.extracted_rna)
    )

    // Collect version info
    ch_versions = ch_versions.mix(BARRNAP.out.versions)

    // PROCESS: Run Blast on predicted ribosomal RNA genes
    BLAST (
        BARRNAP.out.extracted_base.join(CLEAN_READS.out.base_fna),
        blast_db_ch
    )

    // Collect version info
    ch_versions = ch_versions.mix(BLAST.out.versions)

    // PROCESS: Filter Blast output for best score
    FILTER_BLAST (
        BLAST.out.blast_tsv.join(CLEAN_READS.out.base_fna)
    )

    // Collect all BLAST Summaries and concatenate into one file
    blast_summary_ch = blast_summary_ch.mix(FILTER_BLAST.out.blast_summary)
    blast_summary_ch.collectFile(name: 'Summary.16S.tab', storeDir: "${params.outpath}/qa")

    // Collect all BLAST Top Species Summaries and concatenate into one file
    ssu_species_ch = ssu_species_ch.mix(FILTER_BLAST.out.ssu_species)
    ssu_species_ch.collectFile(name: '16S-top-species.tsv', storeDir: "${params.outpath}/ssu")

    // Collect version info
    ch_versions = ch_versions.mix(FILTER_BLAST.out.versions)

    // PROCESS: Run QUAST for quality assessment 
    QA (
        EXTRACT_SINGLETONS.out.gzip_reads.join(CLEAN_READS.out.base_fna)
    )

    // Collect all Assembly Summaries and concatenate into one file
    assembly_summary_ch = assembly_summary_ch.mix(QA.out.summary_assemblies)
    assembly_summary_ch.collectFile(name: 'Summary.Assemblies.tab', keepHeader: true, storeDir: "${params.outpath}/qa")

    // Collect all Cleaned Read/Base Summaries and concatenate into one file
    cleaned_summary_ch = cleaned_summary_ch.mix(QA.out.summary_reads)
    cleaned_summary_ch.collectFile(name: 'Summary.Illumina.CleanedReads-Bases.tab', storeDir: "${params.outpath}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(QA.out.versions)

    // PROCESS: Calculate genome coverage
    GENOME_COVERAGE (
        QA.out.qa_summaries.join(CLEANED_COVERAGE.out.summary_stats)
    )

    // Collect all Genome Coverage Summaries and concatenate into one file
    genome_cov_summary_ch = genome_cov_summary_ch.mix(GENOME_COVERAGE.out.genome_coverage)
    genome_cov_summary_ch.collectFile(name: 'Summary.Illumina.GenomeCoverage.tab', storeDir: "${params.outpath}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(GENOME_COVERAGE.out.versions)
    
    // PATTERN: Collate method version information
    ch_versions.collectFile(name: 'software_versions.yml', storeDir: params.logpath)
}

/*
========================================================================================
                        Completion e-mail and summary                         
========================================================================================
*/

workflow.onComplete {
    log.info """
                |=====================================
                |Pipeline Execution Summary
                |=====================================
                |Workflow Version : ${version}
                |Nextflow Version : ${nextflow.version}
                |Command Line     : ${workflow.commandLine}
                |Resumed          : ${workflow.resume}
                |Completed At     : ${workflow.complete}
                |Duration         : ${workflow.duration}
                |Success          : ${workflow.success}
                |Exit Code        : ${workflow.exitStatus}
                |Launch Dir       : ${workflow.launchDir}
                |=====================================
             """.stripMargin()
}

workflow.onError {
    def err_msg = """
                     |=====================================
                     |Error summary
                     |=====================================
                     |Completed at : ${workflow.complete}
                     |exit status  : ${workflow.exitStatus}
                     |workDir      : ${workflow.workDir}
                     |Error Report :
                     |${workflow.errorReport ?: '-'}
                     |=====================================
                  """.stripMargin()
    log.info err_msg
}


