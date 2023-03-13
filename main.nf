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
      --inpath             Path containing FastQ read sets. Recognized extensions are .fastq, fastq.gz, .fq, fq.gz
      --outpath            Path where results will be saved
    Analysis options:
      --kraken1_db         Path containing Kraken1 database files (*.kdb, *.idk). Default: Mini Kraken.
      --kraken2_db         Path containing Kraken2 database files (*.k2d). Default: Mini Kraken.
      --blast_db           Path to 16S ribosomal database for BLAST. Default: NCBI's 16S ribosomal database.
      --max_memory         Memory (RAM) limit, e.g. '128.GB'
      --max_time           Time (duration) limit for each process, e.g. '240.h'
      --max_cpus           Number of CPUs to use, e.g. 16
    Profile options:
      -profile singularity Use Singularity images to run the workflow. Will pull and convert Docker images from Dockerhub if not locally available.
      -profile docker      Use Docker images to run the workflow. Will pull images from Dockerhub if not locally available.
      -profile conda       TODO: this is not implemented yet.
    Other options:
      -resume              Re-start a workflow using cached results. May not behave as expected with containerization profiles docker or singularity.
      -name                Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

version = "1.1.0"
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

// Validate input parameters
if (!params.inpath) {
    System.err.println "ERROR: inpath parameter must be specified"
    exit 1
}

// Check inpath parameter
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

include { INFILE_HANDLING_UNIX } from "./modules/local/infile_handling_unix/main.nf"
include { REMOVE_PHIX_BBDUK } from "./modules/local/remove_phix_bbduk/main.nf"
include { TRIM_READS_TRIMMOMATIC } from "./modules/local/trim_reads_trimmomatic/main.nf"
//include { TRIM_READS_FASTP } from "./modules/local/trim_reads_fastp/main.nf"
include { OVERLAP_PAIRED_READS_FLASH } from "./modules/local/overlap_paired_reads_flash/main.nf"
//include { OVERLAP_PAIRED_READS_PEAR } from "./modules/local/overlap_paired_reads_pear/main.nf"
include { READ_CLASSIFY_KRAKEN_ONE; READ_CLASSIFY_KRAKEN_TWO; } from "./modules/local/read_classify_kraken/main.nf"
//include { READ_CLASSIFY_CENTRIFUGE } from "./modules/local/read_classify_centrifuge/main.nf"
include { ASSEMBLE_SPADES } from "./modules/local/assemble_spades/main.nf"
//include { ASSEMBLE_SKESA } from "./modules/local/assemble_skesa/main.nf"
include { FILTER_CONTIGS_BIOPYTHON } from "./modules/local/filter_contigs_biopython/main.nf"
include { POLISH_ASSEMBLY_BWA_PILON } from "./modules/local/polish_assembly_bwa_pilon/main.nf"
//include { POLISH_ASSEMBLY_UNICYCLER } from "./modules/local/polish_assembly_unicycler/main.nf"
include { EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS } from "./modules/local/extract_read_alignment_depths_bedtools/main.nf"
include { MLST_MLST } from "./modules/local/mlst_mlst/main.nf"
//include { MLST_SRST2 } from "./modules/local/mlst_srst2/main.nf"
include { ANNOTATE_PROKKA } from "./modules/local/annotate_prokka/main.nf"
//include { ANNOTATE_BACTA } from "./modules/local/annotate_bacta/main.nf"
include { EXTRACT_16S_BIOPYTHON } from "./modules/local/extract_16S_biopython/main.nf"
include { EXTRACT_16S_BARRNAP } from "./modules/local/extract_16S_barrnap/main.nf"
//include { 16S_EXTRACT_RNAMMER } from "./modules/local/16S_extract_rnammer/main.nf"
include { ALIGN_16S_BLAST } from "./modules/local/align_16S_blast/main.nf"
include { BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON } from "./modules/local/best_16S_blastn_bitscore_taxon_python/main.nf"
include { QA_ASSEMBLY_QUAST } from "./modules/local/qa_assembly_quast/main.nf"
//include { QA_ASSEMBLY_BUSCO } from "./modules/local/qa_assembly_busco/main.nf"
//include { QA_ASSEMBLY_CAT } from "./modules/local/qa_assembly_cat/main.nf"
//include { QA_ASSEMBLY_CHECKM2 } from "./modules/local/qa_assembly_checkm2/main.nf"
include { CALCULATE_COVERAGE_UNIX } from "./modules/local/calculate_coverage_unix/main.nf"

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
    // SETUP: Define input
    // Double asterisk looks in specified directory and recursively all subdirectories
    input_ch = Channel.fromFilePairs(params.inpath+'/**_{,R}{1,2}*{fastq,fq}{,.gz}', maxDepth: 2, checkIfExists: true)
    
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
    qc_filecheck_ch = Channel.empty()

    // PROCESS: Read files from input directory, validate and stage input files
    INFILE_HANDLING_UNIX (
        input_ch
    )

    // Collect version info
    ch_versions = ch_versions.mix(INFILE_HANDLING_UNIX.out.versions)

    // PROCESS: Run bbduk to remove PhiX reads
    REMOVE_PHIX_BBDUK (
        INFILE_HANDLING_UNIX.out.input
    )

    // Collect version info
    ch_versions = ch_versions.mix(REMOVE_PHIX_BBDUK.out.versions)

    // PROCESS: Run trimmomatic to clip adapters and do quality trimming
    TRIM_READS_TRIMMOMATIC (
        REMOVE_PHIX_BBDUK.out.phix_removed
    )

    // Collect version info
    ch_versions = ch_versions.mix(TRIM_READS_TRIMMOMATIC.out.versions)

    // PROCESS: Run flash to merge overlapping sister reads into singleton reads
    OVERLAP_PAIRED_READS_FLASH (
        INFILE_HANDLING_UNIX.out.input.join(TRIM_READS_TRIMMOMATIC.out.trimmo)
    )

    // Collect version info
    ch_versions = ch_versions.mix(OVERLAP_PAIRED_READS_FLASH.out.versions)

    // PROCESS: Run kraken1 on paired cleaned reads
    READ_CLASSIFY_KRAKEN_ONE (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads,
        kraken1_db_ch
    )

    // Collect version info
    ch_versions = ch_versions.mix(READ_CLASSIFY_KRAKEN_ONE.out.versions)

    // PROCESS: Run kraken2 on paired cleaned reads
    READ_CLASSIFY_KRAKEN_TWO (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads,
        kraken2_db_ch
    )

    // Collect version info
    ch_versions = ch_versions.mix(READ_CLASSIFY_KRAKEN_TWO.out.versions)

    // PROCESS: Run SPAdes to assemble contigs with cleaned paired reads and cleaned singletons
    ASSEMBLE_SPADES (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads
    )

    // Collect version info
    ch_versions = ch_versions.mix(ASSEMBLE_SPADES.out.versions)

    // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
    FILTER_CONTIGS_BIOPYTHON (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads.join(ASSEMBLE_SPADES.out.contigs)
    )

    // Collect version info
    ch_versions = ch_versions.mix(FILTER_CONTIGS_BIOPYTHON.out.versions)

    // PROCESS: Use BWA/Samtools/Pilon to correct contigs with cleaned PE reads
    POLISH_ASSEMBLY_BWA_PILON (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads.join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
    )

    // Collect version info
    ch_versions = ch_versions.mix(POLISH_ASSEMBLY_BWA_PILON.out.versions)

    // PROCESS: Run Bedtools to extract coverage from the pre-computed BAM alignment file
    EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS (
        POLISH_ASSEMBLY_BWA_PILON.out.bam
    )

    // Collect all Summary Stats and concatenate into one file
    alnstats_summary_ch = alnstats_summary_ch.mix(EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.summary_alnstats)
    alnstats_summary_ch.collectFile(name: 'Summary.Illumina.CleanedReads-AlnStats.tab', storeDir: "${params.outpath}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.versions)

    // PROCESS: Run MLST to find MLST for each polished assembly
    MLST_MLST (
        POLISH_ASSEMBLY_BWA_PILON.out.bam.join(POLISH_ASSEMBLY_BWA_PILON.out.base_fna)
    )

    // Collect all MLST Summaries and concatenate into one file
    mlst_summary_ch = mlst_summary_ch.mix(MLST_MLST.out.summary_mlst)
    mlst_summary_ch.collectFile(name: 'Summary.MLST.tab', storeDir: "${params.outpath}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(MLST_MLST.out.versions)

    // PROCESS: Run Prokka on the polished assembly to annotate the contigs
    ANNOTATE_PROKKA (
        POLISH_ASSEMBLY_BWA_PILON.out.bam.join(POLISH_ASSEMBLY_BWA_PILON.out.base_fna)
    )

    // Collect version info
    ch_versions = ch_versions.mix(ANNOTATE_PROKKA.out.versions)

    // PROCESS: Attempt to extract 16S rRNA gene records from annotation file
    EXTRACT_16S_BIOPYTHON (
        ANNOTATE_PROKKA.out.annotation.join(POLISH_ASSEMBLY_BWA_PILON.out.base_fna)
    )

    // Collect version info
    ch_versions = ch_versions.mix(EXTRACT_16S_BIOPYTHON.out.versions)

    // PROCESS: Extract 16S rRNA gene sequences with Barrnap if missing from 16S_EXTRACT_BIOPYTHON
    EXTRACT_16S_BARRNAP (
        ANNOTATE_PROKKA.out.annotation.join(POLISH_ASSEMBLY_BWA_PILON.out.base_fna).join(EXTRACT_16S_BIOPYTHON.out.extracted_rna)
    )

    // Collect version info
    ch_versions = ch_versions.mix(EXTRACT_16S_BARRNAP.out.versions)

    // PROCESS: Run Blast on predicted 16S ribosomal RNA genes
    ALIGN_16S_BLAST (
        EXTRACT_16S_BARRNAP.out.extracted_base.join(POLISH_ASSEMBLY_BWA_PILON.out.base_fna),
        blast_db_ch
    )

    // Collect version info
    ch_versions = ch_versions.mix(ALIGN_16S_BLAST.out.versions)

    // PROCESS: Filter Blast output for best alignment, based on bitscore
    BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON (
        ALIGN_16S_BLAST.out.blast_tsv.join(POLISH_ASSEMBLY_BWA_PILON.out.base_fna)
    )

    // Collect all BLAST Summaries and concatenate into one file
    blast_summary_ch = blast_summary_ch.mix(BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.blast_summary)
    blast_summary_ch.collectFile(name: 'Summary.16S.tab', storeDir: "${params.outpath}/qa")

    // Collect all BLAST Top Species Summaries and concatenate into one file
    ssu_species_ch = ssu_species_ch.mix(BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.ssu_species)
    ssu_species_ch.collectFile(name: '16S-top-species.tsv', storeDir: "${params.outpath}/ssu")

    // Collect version info
    ch_versions = ch_versions.mix(BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.versions)

    // PROCESS: Run QUAST on the polished assembly for quality assessment and
    //  report the number of cleaned basepairs used to form the assembly
    QA_ASSEMBLY_QUAST (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads.join(POLISH_ASSEMBLY_BWA_PILON.out.base_fna)
    )

    // Collect all Assembly Summaries and concatenate into one file
    assembly_summary_ch = assembly_summary_ch.mix(QA_ASSEMBLY_QUAST.out.summary_assemblies)
    assembly_summary_ch.collectFile(name: 'Summary.Assemblies.tab', keepHeader: true, storeDir: "${params.outpath}/qa")

    // Collect all Cleaned Read/Base Summaries and concatenate into one file
    cleaned_summary_ch = cleaned_summary_ch.mix(QA_ASSEMBLY_QUAST.out.summary_reads)
    cleaned_summary_ch.collectFile(name: 'Summary.Illumina.CleanedReads-Bases.tab', storeDir: "${params.outpath}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(QA_ASSEMBLY_QUAST.out.versions)

    // PROCESS: Calculate genome assembly depth of coverage
    CALCULATE_COVERAGE_UNIX (
        QA_ASSEMBLY_QUAST.out.qa_summaries.join(EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.summary_stats)
    )

    // Collect all Genome Coverage Summaries and concatenate into one file
    genome_cov_summary_ch = genome_cov_summary_ch.mix(CALCULATE_COVERAGE_UNIX.out.genome_coverage)
    genome_cov_summary_ch.collectFile(name: 'Summary.Illumina.GenomeCoverage.tab', storeDir: "${params.outpath}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(CALCULATE_COVERAGE_UNIX.out.versions)
    
    // PATTERN: Collate method version information
    ch_versions.collectFile(name: 'software_versions.yml', storeDir: params.logpath)

    // Collect all QC File Checks and concatenate into one file
    qc_filecheck_ch = qc_filecheck_ch.concat(
        INFILE_HANDLING_UNIX.out.qc_filecheck,
        REMOVE_PHIX_BBDUK.out.qc_phix_genome_filecheck,
        REMOVE_PHIX_BBDUK.out.qc_phix_removed_filecheck,
        TRIM_READS_TRIMMOMATIC.out.qc_adapters_filecheck,
        TRIM_READS_TRIMMOMATIC.out.qc_removed_adapters_filecheck,
        OVERLAP_PAIRED_READS_FLASH.out.qc_filecheck,
        ASSEMBLE_SPADES.out.qc_filecheck,
        POLISH_ASSEMBLY_BWA_PILON.out.filtered_asm_filecheck,
        POLISH_ASSEMBLY_BWA_PILON.out.pe_alignment_filecheck,
        POLISH_ASSEMBLY_BWA_PILON.out.polished_asm_filecheck,
        POLISH_ASSEMBLY_BWA_PILON.out.corrected_asm_filecheck,
        POLISH_ASSEMBLY_BWA_PILON.out.se_alignment_filecheck,
        ANNOTATE_PROKKA.out.qc_filecheck,
        EXTRACT_16S_BARRNAP.out.qc_ssu_extracted_filecheck,
        EXTRACT_16S_BARRNAP.out.qc_ssu_renamed_filecheck,
        ALIGN_16S_BLAST.out.qc_filecheck,
        BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.qc_filecheck
    )
    qc_filecheck_ch.collectFile(name: 'Summary.QC_File_Checks.tab', storeDir: "${params.outpath}/qa", sort: {it.getSimpleName()})
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
                     |Error Summary
                     |=====================================
                     |Completed at     : ${workflow.complete}
                     |Exit Code        : ${workflow.exitStatus}
                     |Nextflow workDir : ${workflow.workDir}
                     |Error Report     :
                     |${workflow.errorReport ?: '-'}
                     |=====================================
                  """.stripMargin()
    log.info err_msg
}