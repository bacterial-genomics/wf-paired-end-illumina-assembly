/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAssembly.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.kraken1_db, params.kraken2_db, params.blast_db, params.gtdb_db, params.busco_db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet or directory not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// CONFIGS: Import configs for this workflow
//

// BUSCO config
if(params.busco_config){
    ch_busco_config_file = Channel.fromPath( "${params.busco_config}" )
} else {
    ch_busco_config_file = []
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { INFILE_HANDLING_UNIX                    } from "../modules/local/infile_handling_unix/main"
include { VALIDATE_FASTQ_SEQFU                    } from "../modules/local/validate_fastq_seqfu/main.nf"

include { CALCULATE_METRICS_FASTQ_SEQKIT as CALC_STATS_INPUT_FQ_SEQKIT } from "../modules/local/calculate_metrics_fastq_seqkit/main"
include { CALCULATE_METRICS_FASTQ_SEQKIT as CALC_STATS_NOPHIX_FQ_SEQKIT } from "../modules/local/calculate_metrics_fastq_seqkit/main"
include { CALCULATE_METRICS_FASTQ_SEQKIT as CALC_STATS_TRIM_FQ_SEQKIT } from "../modules/local/calculate_metrics_fastq_seqkit/main"
include { CALCULATE_METRICS_FASTQ_SEQKIT as CALC_STATS_CLEANEDREADS_FQ_SEQKIT } from "../modules/local/calculate_metrics_fastq_seqkit/main"
include { CALCULATE_METRICS_FASTQ_SEQTK           } from "../modules/local/calculate_metrics_fastq_seqtk/main"

include { REMOVE_PHIX_BBDUK                       } from "../modules/local/remove_phix_bbduk/main"
include { TRIM_READS_TRIMMOMATIC                  } from "../modules/local/trim_reads_trimmomatic/main"
include { TRIM_READS_FASTP                        } from "../modules/local/trim_reads_fastp/main"
include { OVERLAP_PAIRED_READS_FLASH              } from "../modules/local/overlap_paired_reads_flash/main"
//include { OVERLAP_PAIRED_READS_PEAR             } from "../modules/local/overlap_paired_reads_pear/main"

include { KRAKEN1_DB_PREPARATION_UNIX             } from "../modules/local/kraken1_db_preparation_unix/main"
include { READ_CLASSIFY_KRAKEN_ONE                } from "../modules/local/read_classify_kraken/main"
include { KRAKEN2_DB_PREPARATION_UNIX             } from "../modules/local/kraken2_db_preparation_unix/main"
include { READ_CLASSIFY_KRAKEN_TWO                } from "../modules/local/read_classify_kraken2/main"
// include { READ_CLASSIFY_CENTRIFUGE                } from "../modules/local/read_classify_centrifuge/main"
// include { READ_CLASSIFY_KRAKENUNIQ                } from "../modules/local/read_classify_krakenuniq/main"
// include { READ_CLASSIFY_METAPHLAN                 } from "../modules/local/read_classify_metaphlan/main"

include { EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS  } from "../modules/local/extract_read_alignment_depths_bedtools/main"

include { MLST_MLST                               } from "../modules/local/mlst_mlst/main"
//include { MLST_SRST2                            } from "../modules/local/mlst_srst2/main"

include { ANNOTATE_PROKKA                         } from "../modules/local/annotate_prokka/main"
//include { ANNOTATE_BAKTA                        } from "../modules/local/annotate_bakta/main"

include { EXTRACT_16S_BIOPYTHON                   } from "../modules/local/extract_16S_biopython/main"
include { EXTRACT_16S_BARRNAP                     } from "../modules/local/extract_16S_barrnap/main"
//include { 16S_EXTRACT_RNAMMER                   } from "../modules/local/16S_extract_rnammer/main"
include { BLAST_DB_PREPARATION_UNIX               } from "../modules/local/blast_db_preparation_unix/main"
include { ALIGN_16S_BLAST                         } from "../modules/local/align_16S_blast/main"
include { BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON   } from "../modules/local/best_16S_blastn_bitscore_taxon_python/main"
include { CLASSIFY_16S_RDP                        } from "../modules/local/classify_16S_rdp/main"
include { SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON     } from "../modules/local/split_multifasta_assembly_biopython/main"

include { CONVERT_TSV_TO_EXCEL_PYTHON             } from "../modules/local/convert_tsv_to_excel_python/main"
include { CREATE_EXCEL_RUN_SUMMARY_PYTHON         } from "../modules/local/create_excel_run_summary_python/main"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                             } from "../subworkflows/local/input_check"
include { HOST_REMOVAL                            } from "../subworkflows/local/host_removal"
include { DOWNSAMPLE                              } from "../subworkflows/local/downsampling"
include { ASSEMBLE_CONTIGS                        } from "../subworkflows/local/assemble_contigs"
include { ASSEMBLY_ASSESSMENT                     } from "../subworkflows/local/assembly_assessment"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Assembler input
if ( toLower(params.assembler) == "skesa" ) {
    var_assembler_name = "SKESA"
} else {
    var_assembler_name = "SPAdes"
}

// NCBI's SRA Human Scrubber
if (params.sra_scrubber_db) {
    ch_sra_scrubber_db_file = file(params.sra_scrubber_db, checkIfExists: true)
} else {
    ch_sra_scrubber_db_file = []
}

// PhiX Reference
if (params.phix_reference) {
    ch_phix_reference = Channel
                            .fromPath(params.phix_reference, checkIfExists: true)
                            .map{
                                file ->
                                    if ( file.extension in ['fasta', 'fas', 'fa', 'fna'] ) {
                                        [ file ]
                                    } else {
                                        error("PhiX reference file not in supported extension: .fasta, .fas, .fa, .fna")
                                    }
                            }
                            .collect()
} else {
    error("Path to PhiX reference not specified. Please supply a PhiX reference file in FastA format via `--phix_reference` parameter.")
}

// Adapter Reference
if (params.adapter_reference) {
    ch_adapter_reference = Channel
                            .fromPath(params.adapter_reference, checkIfExists: true)
                            .map{
                                file ->
                                    if ( file.extension in ['fasta', 'fas', 'fa', 'fna'] ) {
                                        [ file ]
                                    } else {
                                        error("Adapter reference file not in supported extension: .fasta, .fas, .fa, .fna")
                                    }
                            }
                            .collect()
} else {
    ch_adapter_reference = []
}

// CAT
if (params.cat_db) {
    ch_cat_db_file = file(params.cat_db, checkIfExists: true)
} else {
    ch_cat_db_file = []
}

// CheckM2
if (params.checkm2_db) {
    ch_checkm2_db_file = file(params.checkm2_db, checkIfExists: true)
} else {
    ch_checkm2_db_file = []
}

// GTDB
if (params.gtdb_db) {
    ch_gtdbtk_db_file = file(params.gtdb_db, checkIfExists: true)
} else {
    ch_gtdbtk_db_file = []
}

// Mash database for GTDB-Tk
if (params.mash_db) {
    ch_mash_db_file = file(params.mash_db)
} else {
    ch_mash_db_file = []
}

// BUSCO
if (params.busco_db) {
    ch_busco_db_file = file(params.busco_db, checkIfExists: true)
} else {
    ch_busco_db_file = []
}

// kraken
if (params.kraken1_db) {
    ch_kraken1_db_file = file(params.kraken1_db, checkIfExists: true)
} else {
    ch_kraken1_db_file = []
}

// kraken2
if (params.kraken2_db) {
    ch_kraken2_db_file = file(params.kraken2_db, checkIfExists: true)
} else {
    ch_kraken2_db_file = []
}

// NCBI BLAST
if (params.blast_db) {
    ch_blast_db_file = file(params.blast_db, checkIfExists: true)
} else {
    ch_blast_db_file = Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert params.assembler to lowercase
def toLower(it) {
    it.toString().toLowerCase()
}

// Check QC filechecks for a failure
def qcfilecheck(process, qcfile, inputfile) {
    qcfile.map{ meta, file -> [ meta, [file] ] }
            .join(inputfile)
            .map{ meta, qc, input ->
                data = []
                qc.flatten().each{ data += it.readLines() }

                if ( data.any{ it.contains('FAIL') } ) {
                    line = data.last().split('\t')
                    log.warn("${line[1]} QC check failed during process ${process} for sample ${line.first()}")
                } else {
                    [ meta, input ]
                }
            }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY {

    // SETUP: Define empty channel
    ch_versions             = Channel.empty()
    ch_qc_filecheck         = Channel.empty()
    ch_output_summary_files = Channel.empty()

    /*
    ================================================================================
                            Preprocessing and Cleaning FastQ files
    ================================================================================
    */

    // SUBWORKFLOW: Check input for samplesheet or pull inputs from directory
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Check input files meet size criteria
    INFILE_HANDLING_UNIX (
        INPUT_CHECK.out.raw_reads
    )
    // INFILE_HANDLING_UNIX
    //     .view { file -> println "DEBUG: From INFILE_HANDLING_UNIX, emitting file: ${file}" }
    // INFILE_HANDLING_UNIX
    //     .view { item -> println "DEBUG: From INFILE_HANDLING_UNIX, channel item: ${item}" }
    ch_versions        = ch_versions.mix(INFILE_HANDLING_UNIX.out.versions)
    ch_qc_filecheck    = ch_qc_filecheck.concat(INFILE_HANDLING_UNIX.out.qc_filecheck)
    ch_infile_handling = qcfilecheck(
                            "INFILE_HANDLING_UNIX",
                            INFILE_HANDLING_UNIX.out.qc_filecheck,
                            INFILE_HANDLING_UNIX.out.input
                          )

    ch_infile_handling = ch_infile_handling
                            .map{
                                meta, file ->
                                    meta['assembler'] = "${var_assembler_name}"
                                    [ meta, file]
                            }

    ch_infile_checksum = INFILE_HANDLING_UNIX.out.checksums
                            .collectFile(
                                name:       "Summary.Input_Checksums.tsv",
                                keepHeader: true,
                                sort:       { file -> file.text },
                                storeDir:   "${params.outdir}/Summaries"
                            )
                            .view { collectedFiles -> println "DEBUG: From INFILE_HANDLING_UNIX.out.checksums, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_infile_checksum)

    // PROCESS: Calculate input FastQ metrics for each sample with SeqKit
    CALC_STATS_INPUT_FQ_SEQKIT (
        INPUT_CHECK.out.raw_reads,
        "Input_Reads"
    )
    ch_versions = ch_versions.mix(CALC_STATS_INPUT_FQ_SEQKIT.out.versions)

    // Collect raw input read/base summaries and concatenate into one file
    ch_input_reads_metrics_summary = CALC_STATS_INPUT_FQ_SEQKIT.out.output
                                        .collectFile(
                                            name:       "Summary.Input_Reads.Metrics.tsv",
                                            keepHeader: true,
                                            sort:       { file -> file.text },
                                            storeDir:   "${params.outdir}/Summaries"
                                        )
    ch_output_summary_files = ch_output_summary_files.mix(ch_input_reads_metrics_summary)

    // Check input files are valid FastQ format
    VALIDATE_FASTQ_SEQFU (
        INPUT_CHECK.out.raw_reads
    )
    ch_versions              = ch_versions.mix(VALIDATE_FASTQ_SEQFU.out.versions)
    ch_qc_filecheck          = ch_qc_filecheck.concat(VALIDATE_FASTQ_SEQFU.out.qc_filecheck)
    ch_validate_fastq_format = qcfilecheck(
                                  "VALIDATE_FASTQ_SEQFU",
                                  VALIDATE_FASTQ_SEQFU.out.qc_filecheck,
                                  VALIDATE_FASTQ_SEQFU.out.input
                                )

    // SUBWORKFLOW: Remove host from FastQ files
    HOST_REMOVAL (
        ch_infile_handling,
        ch_sra_scrubber_db_file
    )
    // HOST_REMOVAL
    //     .view { file -> println "DEBUG: From HOST_REMOVAL, emitting file: ${file}" }
    // HOST_REMOVAL
    //     .view { item -> println "DEBUG: From HOST_REMOVAL, channel item: ${item}" }
    ch_versions             = ch_versions.mix(HOST_REMOVAL.out.versions)
    ch_output_summary_files = ch_output_summary_files.mix(HOST_REMOVAL.out.output_summary_files)

    // SUBWORKFLOW: Downsample FastQ files
    DOWNSAMPLE (
        HOST_REMOVAL.out.host_removed_reads
    )
    ch_versions             = ch_versions.mix(DOWNSAMPLE.out.versions)
    ch_output_summary_files = ch_output_summary_files.mix(DOWNSAMPLE.out.output_summary_files)

    // PROCESS: Run bbduk to remove PhiX reads
    REMOVE_PHIX_BBDUK (
        DOWNSAMPLE.out.reads,
        ch_phix_reference
    )
    // REMOVE_PHIX_BBDUK
    //     .view { file -> println "DEBUG: From REMOVE_PHIX_BBDUK, emitting file: ${file}" }
    // REMOVE_PHIX_BBDUK
    //     .view { item -> println "DEBUG: From REMOVE_PHIX_BBDUK, channel item: ${item}" }
    ch_versions     = ch_versions.mix(REMOVE_PHIX_BBDUK.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(REMOVE_PHIX_BBDUK.out.qc_filecheck)
    ch_removed_phix = qcfilecheck(
                        "REMOVE_PHIX_BBDUK",
                        REMOVE_PHIX_BBDUK.out.qc_filecheck,
                        REMOVE_PHIX_BBDUK.out.phix_removed_reads
                      )

    // Collect PhiX removal summaries and concatenate into one file
    ch_phix_removal_summary = REMOVE_PHIX_BBDUK.out.summary
                                .collectFile(
                                    name:       "Summary.PhiX_Removal.tsv",
                                    keepHeader: true,
                                    sort:       { file -> file.text },
                                    storeDir:   "${params.outdir}/Summaries"
                                )
                                .view { collectedFiles -> println "DEBUG: From REMOVE_PHIX_BBDUK.out.summary, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_phix_removal_summary)

    // PROCESS: Calculate PhiX-free FastQ metrics for each sample with SeqKit
    CALC_STATS_NOPHIX_FQ_SEQKIT (
        REMOVE_PHIX_BBDUK.out.phix_removed_reads,
        "PhiX_Removed_Reads"
    )
    ch_versions = ch_versions.mix(CALC_STATS_NOPHIX_FQ_SEQKIT.out.versions)

    // Collect PhiX-free read/base summaries and concatenate into one file
    ch_nophix_reads_metrics_summary = CALC_STATS_NOPHIX_FQ_SEQKIT.out.output
                                        .collectFile(
                                            name:       "Summary.PhiX_Removed_Reads.Metrics.tsv",
                                            keepHeader: true,
                                            sort:       { file -> file.text },
                                            storeDir:   "${params.outdir}/Summaries"
                                        )
    ch_output_summary_files = ch_output_summary_files.mix(ch_nophix_reads_metrics_summary)

    if ( toLower(params.trim_reads_tool) == "trimmomatic" ) {
        // PROCESS: Run trimmomatic to clip adapters and do quality trimming
        TRIM_READS_TRIMMOMATIC (
            ch_removed_phix,
            ch_adapter_reference
        )
        // TRIM_READS_TRIMMOMATIC
        //     .view { file -> println "DEBUG: From TRIM_READS_TRIMMOMATIC, emitting file: ${file}" }
        // TRIM_READS_TRIMMOMATIC
        //     .view { item -> println "DEBUG: From TRIM_READS_TRIMMOMATIC, channel item: ${item}" }
        ch_versions     = ch_versions.mix(TRIM_READS_TRIMMOMATIC.out.versions)
        ch_qc_filecheck = ch_qc_filecheck.concat(TRIM_READS_TRIMMOMATIC.out.qc_filecheck)
        ch_trim_reads   = qcfilecheck(
                            "TRIM_READS_TRIMMOMATIC",
                            TRIM_READS_TRIMMOMATIC.out.qc_filecheck,
                            TRIM_READS_TRIMMOMATIC.out.fastq_adapters_removed
                            )

        // Collect read trimming summaries and concatenate into one file
        ch_trimmomatic_summary = TRIM_READS_TRIMMOMATIC.out.summary
                                    .collectFile(
                                        name:       "Summary.Adapter_and_QC_Trim.tsv",
                                        keepHeader: true,
                                        sort:       { file -> file.text },
                                        storeDir:   "${params.outdir}/Summaries"
                                    )
                                    .view { collectedFiles -> println "DEBUG: From TRIM_READS_TRIMMOMATIC.out.summary, collected files: ${collectedFiles}" }
        ch_output_summary_files = ch_output_summary_files.mix(ch_trimmomatic_summary)

        // PROCESS: Calculate adapter-and-quality-trimmed FastQ metrics for each sample with SeqKit
        CALC_STATS_TRIM_FQ_SEQKIT (
            TRIM_READS_TRIMMOMATIC.out.fastq_adapters_removed,
            "Adapter_QC_Trim_Reads"
        )
        ch_versions = ch_versions.mix(CALC_STATS_TRIM_FQ_SEQKIT.out.versions)

        // Collect adapter-and-quality-trimmed read/base summaries and concatenate into one file
        ch_trim_reads_metrics_summary = CALC_STATS_TRIM_FQ_SEQKIT.out.output
                                            .collectFile(
                                                name:       "Summary.Adapter_QC_Trim_Reads.Metrics.tsv",
                                                keepHeader: true,
                                                sort:       { file -> file.text },
                                                storeDir:   "${params.outdir}/Summaries"
                                            )
        ch_output_summary_files = ch_output_summary_files.mix(ch_trim_reads_metrics_summary)

    } else if ( toLower(params.trim_reads_tool) == "fastp" ) {
        // Do not use 'adapters_Nextera_NEB_TruSeq_NuGEN_ThruPLEX.fas' for fastp
        ch_adapter_reference = ch_adapter_reference.map{ it[0].getSimpleName() == "adapters_Nextera_NEB_TruSeq_NuGEN_ThruPLEX.fas" }
                                ? []
                                : ch_adapter_reference

        TRIM_READS_FASTP (
            ch_removed_phix,
            ch_adapter_reference
        )
        // TRIM_READS_FASTP
        //     .view { file -> println "DEBUG: From TRIM_READS_FASTP, emitting file: ${file}" }
        // TRIM_READS_FASTP
        //     .view { item -> println "DEBUG: From TRIM_READS_FASTP, channel item: ${item}" }
        ch_versions     = ch_versions.mix(TRIM_READS_FASTP.out.versions)
        ch_qc_filecheck = ch_qc_filecheck.concat(TRIM_READS_FASTP.out.qc_filecheck)
        ch_trim_reads   = qcfilecheck(
                            "TRIM_READS_FASTP",
                            TRIM_READS_FASTP.out.qc_filecheck,
                            TRIM_READS_FASTP.out.fastq_adapters_removed
                          )

        ch_fastp_summary = TRIM_READS_FASTP.out.summary
                                .collectFile(
                                    name:       "Summary.Adapter_and_QC_Trim.tsv",
                                    keepHeader: true,
                                    sort:       { file -> file.text },
                                    storeDir:   "${params.outdir}/Summaries"
                                )
                                .view { collectedFiles -> println "DEBUG: From TRIM_READS_FASTP.out.summary, collected files: ${collectedFiles}" }
        ch_output_summary_files = ch_output_summary_files.mix(ch_fastp_summary)

        // PROCESS: Calculate adapter-and-quality-trimmed FastQ metrics for each sample with SeqKit
        CALC_STATS_TRIM_FQ_SEQKIT (
            TRIM_READS_FASTP.out.fastq_adapters_removed,
            "Adapter_QC_Trim_Reads"
        )
        ch_versions = ch_versions.mix(CALC_STATS_TRIM_FQ_SEQKIT.out.versions)

        // Collect adapter-and-quality-trimmed read/base summaries and concatenate into one file
        ch_trim_reads_metrics_summary = CALC_STATS_TRIM_FQ_SEQKIT.out.output
                                            .collectFile(
                                                name:       "Summary.Adapter_QC_Trim_Reads.Metrics.tsv",
                                                keepHeader: true,
                                                sort:       { file -> file.text },
                                                storeDir:   "${params.outdir}/Summaries"
                                            )
        ch_output_summary_files = ch_output_summary_files.mix(ch_trim_reads_metrics_summary)
    }

    // PROCESS: Run flash to merge overlapping sister reads into singleton reads
    OVERLAP_PAIRED_READS_FLASH (
        ch_trim_reads
    )
    // OVERLAP_PAIRED_READS_FLASH
    //     .view { file -> println "DEBUG: From OVERLAP_PAIRED_READS_FLASH, emitting file: ${file}" }
    // OVERLAP_PAIRED_READS_FLASH
    //     .view { item -> println "DEBUG: From OVERLAP_PAIRED_READS_FLASH, channel item: ${item}" }
    ch_versions      = ch_versions.mix(OVERLAP_PAIRED_READS_FLASH.out.versions)
    ch_qc_filecheck  = ch_qc_filecheck.concat(OVERLAP_PAIRED_READS_FLASH.out.qc_filecheck)
    ch_overlap_flash = qcfilecheck(
                            "OVERLAP_PAIRED_READS_FLASH",
                            OVERLAP_PAIRED_READS_FLASH.out.qc_filecheck,
                            OVERLAP_PAIRED_READS_FLASH.out.cleaned_fastq_files
                        )

    ch_cleanedreads_checksum = OVERLAP_PAIRED_READS_FLASH.out.checksums
                            .collectFile(
                                name:       "Summary.Clean_Reads_Checksums.tsv",
                                keepHeader: true,
                                sort:       { file -> file.text },
                                storeDir:   "${params.outdir}/Summaries"
                            )
                            .view { collectedFiles -> println "DEBUG: From OVERLAP_PAIRED_READS_FLASH.out.checksums, collected files: ${collectedFiles}" }
    ch_output_summary_files  = ch_output_summary_files.mix(ch_cleanedreads_checksum)

    // Collect singleton read summaries and concatenate into one file
    ch_overlap_summary = OVERLAP_PAIRED_READS_FLASH.out.summary
                                .collectFile(
                                    name:       "Summary.Clean_and_Overlapped.tsv",
                                    keepHeader: true,
                                    sort:       { file -> file.text },
                                    storeDir:   "${params.outdir}/Summaries"
                                )
                                .view { collectedFiles -> println "DEBUG: From OVERLAP_PAIRED_READS_FLASH.out.summary, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_overlap_summary)

    /*
    ================================================================================
                                    Clean_Reads Assessment
    ================================================================================
    */
    // PROCESS: Calculate cleaned FastQ metrics for each sample with SeqKit
    CALC_STATS_CLEANEDREADS_FQ_SEQKIT (
        OVERLAP_PAIRED_READS_FLASH.out.cleaned_fastq_files,
        "Clean_Reads"
    )
    ch_versions = ch_versions.mix(CALC_STATS_CLEANEDREADS_FQ_SEQKIT.out.versions)

    // Collect cleaned read/base summaries and concatenate into one file
    ch_cleaned_reads_metrics_summary = CALC_STATS_CLEANEDREADS_FQ_SEQKIT.out.output
                                        .collectFile(
                                            name:       "Summary.Clean_Reads.Metrics.tsv",
                                            keepHeader: true,
                                            sort:       { file -> file.text },
                                            storeDir:   "${params.outdir}/Summaries"
                                        )
    ch_output_summary_files = ch_output_summary_files.mix(ch_cleaned_reads_metrics_summary)

    // NOTE: This is just a backup working solution; use SeqKit only for now due to extra metrics provided
    // // PROCESS: Calculate cleaned FastQ metrics for each sample with SeqKit
    // CALCULATE_METRICS_FASTQ_SEQTK (
    //     OVERLAP_PAIRED_READS_FLASH.out.cleaned_fastq_files
    // )
    // ch_versions = ch_versions.mix(CALCULATE_METRICS_FASTQ_SEQTK.out.versions)
    // // Collect cleaned read/base summaries and concatenate into one file
    // ch_cleaned_reads_metrics_summary = CALCULATE_METRICS_FASTQ_SEQTK.out.output
    //                                     .collectFile(
    //                                         name:       "Summary.Cleaned_Reads.Metrics.tsv",
    //                                         keepHeader: true,
    //                                         sort:       { file -> file.text },
    //                                         storeDir:   "${params.outdir}/Summaries"
    //                                     )
    //                                     .view { collectedFiles -> println "DEBUG: From ch_cleaned_reads_metrics_summary, collected files: ${collectedFiles}" }
    // ch_output_summary_files = ch_output_summary_files.mix(ch_cleaned_reads_metrics_summary)

    // KAT K-mer plot

    // seqkit stats for TSV

    /*
    ================================================================================
                                    Taxonomic information
    ================================================================================
    */

    // Prepare kraken1 database for use
    if ( ch_kraken1_db_file ) {
        if ( ch_kraken1_db_file.extension in ['gz', 'tgz'] ) {
            // Add meta information
            ch_kraken1_db = Channel.of(ch_kraken1_db_file)
                                .map{
                                    db ->
                                        def meta = [:]
                                        meta['id'] = db.getSimpleName()
                                        [ meta, db ]
                                }
            // Expects to be .tar.gz!
            KRAKEN1_DB_PREPARATION_UNIX (
                ch_kraken1_db
            )
            ch_versions       = ch_versions.mix(KRAKEN1_DB_PREPARATION_UNIX.out.versions)
            ch_db_for_kraken1 = KRAKEN1_DB_PREPARATION_UNIX.out.db.collect()

        } else if ( ch_kraken1_db_file.isDirectory() ) {
            ch_db_for_kraken1 = Channel
                                    .fromPath(
                                        "${ch_kraken1_db_file}/{database}.{idx,kdb}",
                                        checkIfExists: true
                                    )
                                    .combine(
                                        Channel
                                            .fromPath(
                                                "${ch_kraken1_db_file}/taxonomy/{names,nodes}.dmp",
                                                checkIfExists: true
                                            )
                                    )
                                    .collect()
                                    .map{
                                        file ->
                                            if (file.size() >= 4) {
                                                [ file[0].getParent() ]
                                            } else {
                                                error("Kraken requires 'database.{idx,kdb}' and 'taxonomy/{names,nodes}.dmp' files!")
                                            }
                                    }
                                    .collect()

        } else {
            log.error("Unsupported object given to --kraken1_db, database must be supplied as either a directory or a .tar.gz file!")
            ch_db_for_kraken1 = Channel.empty()
        }

    } else {
        log.warn("Kraken could not be performed - database not specified using --kraken1_db!")
        ch_db_for_kraken1 = Channel.empty()
    }

    // PROCESS: Run kraken1 on paired cleaned reads
    READ_CLASSIFY_KRAKEN_ONE (
        ch_overlap_flash,
        ch_db_for_kraken1
    )
    // READ_CLASSIFY_KRAKEN_ONE
    //     .view { file -> println "DEBUG: From READ_CLASSIFY_KRAKEN_ONE, emitting file: ${file}" }
    // READ_CLASSIFY_KRAKEN_ONE
    //     .view { item -> println "DEBUG: From READ_CLASSIFY_KRAKEN_ONE, channel item: ${item}" }
    ch_versions = ch_versions.mix(READ_CLASSIFY_KRAKEN_ONE.out.versions)

    // Collect kraken summaries and concatenate into one file
    ch_kraken_one_summary = READ_CLASSIFY_KRAKEN_ONE.out.summary
                                .collectFile(
                                    name:       "Summary.Kraken.tsv",
                                    keepHeader: true,
                                    sort:       { file -> file.text },
                                    storeDir:   "${params.outdir}/Summaries"
                                )
                                .view { collectedFiles -> println "DEBUG: From READ_CLASSIFY_KRAKEN_ONE.out.summary, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_kraken_one_summary)

    // Prepare kraken2 database for use
    if ( ch_kraken2_db_file ) {
        if ( ch_kraken2_db_file.extension in ['gz', 'tgz'] ) {
            // Add meta information
            ch_kraken2_db = Channel.of(ch_kraken2_db_file)
                                .map{
                                    db ->
                                        def meta = [:]
                                        meta['id'] = db.getSimpleName()
                                        [ meta, db ]
                                }
            // Expects to be .tar.gz!
            KRAKEN2_DB_PREPARATION_UNIX (
                ch_kraken2_db
            )
            ch_versions       = ch_versions.mix(KRAKEN2_DB_PREPARATION_UNIX.out.versions)
            ch_db_for_kraken2 = KRAKEN2_DB_PREPARATION_UNIX.out.db.collect()

        } else if ( ch_kraken2_db_file.isDirectory() ) {
            ch_db_for_kraken2 = Channel
                                    .fromPath( "${ch_kraken2_db_file}/*.k2d" )
                                    .collect()
                                    .map{
                                        file ->
                                            if (file.size() >= 3) {
                                                [ file[0].getParent() ]
                                            } else {
                                                error("Kraken2 requires '{hash,opts,taxo}.k2d' files!")
                                            }
                                    }
                                    .collect()

        } else {
            log.error("Unsupported object given to --kraken2_db, database must be supplied as either a directory or a .tar.gz file!")
            ch_db_for_kraken2 = Channel.empty()
        }

    } else {
        log.warn("Kraken2 could not be performed - database not specified using --kraken2_db!")
        ch_db_for_kraken2 = Channel.empty()
    }

    // PROCESS: Run kraken2 on paired cleaned reads
    READ_CLASSIFY_KRAKEN_TWO (
        ch_overlap_flash,
        ch_db_for_kraken2
    )
    // READ_CLASSIFY_KRAKEN_TWO
    //     .view { file -> println "DEBUG: From READ_CLASSIFY_KRAKEN_TWO, emitting file: ${file}" }
    // READ_CLASSIFY_KRAKEN_TWO
    //     .view { item -> println "DEBUG: From READ_CLASSIFY_KRAKEN_TWO, channel item: ${item}" }
    ch_versions = ch_versions.mix(READ_CLASSIFY_KRAKEN_TWO.out.versions)

    // Collect kraken2 summaries and concatenate into one file
    ch_kraken_two_summary = READ_CLASSIFY_KRAKEN_TWO.out.summary
                                .collectFile(
                                    name:       "Summary.Kraken2.tsv",
                                    keepHeader: true,
                                    sort:       { file -> file.text },
                                    storeDir:   "${params.outdir}/Summaries"
                                )
                                .view { collectedFiles -> println "DEBUG: From READ_CLASSIFY_KRAKEN_TWO.out.summary, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_kraken_two_summary)


    /*
    ================================================================================
                                    Assembly
    ================================================================================
    */

    ASSEMBLE_CONTIGS (
        ch_overlap_flash,
        var_assembler_name
    )
    // ASSEMBLE_CONTIGS
    //     .view { file -> println "DEBUG: From ASSEMBLE_CONTIGS, emitting file: ${file}" }
    // ASSEMBLE_CONTIGS
    //     .view { item -> println "DEBUG: From ASSEMBLE_CONTIGS, channel item: ${item}" }
    ch_versions     = ch_versions.mix(ASSEMBLE_CONTIGS.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(ASSEMBLE_CONTIGS.out.qc_filecheck)

    ch_assembly_checksum = ASSEMBLE_CONTIGS.out.checksums
                            .collectFile(
                                name:       "Summary.Assembly_Checksums.tsv",
                                keepHeader: true,
                                sort:       { file -> file.text },
                                storeDir:   "${params.outdir}/Summaries"
                            )
                            .view { collectedFiles -> println "DEBUG: From OVERLAP_PAIRED_READS_FLASH.out.checksums, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_assembly_checksum)

    /*
    ================================================================================
                            Assembly Information
    ================================================================================
    */

    // PROCESS: Run Bedtools to extract coverage from the pre-computed BAM alignment file
    EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS (
        ASSEMBLE_CONTIGS.out.bam_files
    )
    // EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS
    //     .view { file -> println "DEBUG: From EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS, emitting file: ${file}" }
    // EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS
    //     .view { item -> println "DEBUG: From EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS, channel item: ${item}" }
    ch_versions = ch_versions.mix(EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.versions)

    // Collect alignment summary stats and concatenate into one file
    ch_alignment_stats_summary = EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.summary
                                    .map{ meta, file -> file }
                                    .collectFile(
                                        name:     "Summary.Clean_Reads_Aligned.tsv",
                                        keepHeader: true,
                                        sort:      { file -> file.text },
                                        storeDir: "${params.outdir}/Summaries"
                                    )
                                .view { collectedFiles -> println "DEBUG: From EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.summary, collected files: ${collectedFiles}" }
    ch_output_summary_files    = ch_output_summary_files.mix(ch_alignment_stats_summary)

    // PROCESS: Run MLST to find MLST for each polished assembly
    MLST_MLST (
        ASSEMBLE_CONTIGS.out.assembly_file
    )
    // MLST_MLST
    //     .view { file -> println "DEBUG: From MLST_MLST, emitting file: ${file}" }
    // MLST_MLST
    //     .view { item -> println "DEBUG: From MLST_MLST, channel item: ${item}" }
    ch_versions = ch_versions.mix(MLST_MLST.out.versions)

    // Collect MLST Summaries and concatenate into one file
    ch_mlst_summary = MLST_MLST.out.summary
                        .collectFile(
                            name:     "Summary.MLST.tsv",
                            keepHeader: true,
                            sort:      { file -> file.text },
                            storeDir: "${params.outdir}/Summaries"
                        )
                        .view { collectedFiles -> println "DEBUG: From MLST_MLST.out.summary, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_mlst_summary)

    // PROCESS: Annotate the polished assembly using Prokka
    ANNOTATE_PROKKA (
        ASSEMBLE_CONTIGS.out.assembly_file
    )
    // ANNOTATE_PROKKA
    //     .view { file -> println "DEBUG: From ANNOTATE_PROKKA, emitting file: ${file}" }
    // ANNOTATE_PROKKA
    //     .view { item -> println "DEBUG: From ANNOTATE_PROKKA, channel item: ${item}" }
    ch_versions = ch_versions.mix(ANNOTATE_PROKKA.out.versions)
    ch_genbank  = qcfilecheck(
                    "ANNOTATE_PROKKA",
                    ANNOTATE_PROKKA.out.qc_filecheck,
                    ANNOTATE_PROKKA.out.prokka_genbank_file
                )

    ch_annotation_checksum = ANNOTATE_PROKKA.out.checksums
                            .collectFile(
                                name:       "Summary.Annotation_Checksums.tsv",
                                keepHeader: true,
                                sort:       { file -> file.text },
                                storeDir:   "${params.outdir}/Summaries"
                            )
                            .view { collectedFiles -> println "DEBUG: From ANNOTATE_PROKKA.out.checksums, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_annotation_checksum)


    /*
    ================================================================================
                            Evaluate 16S
    ================================================================================
    */

    // PROCESS: Attempt to extract 16S rRNA gene records from prokka_genbank_file file
    EXTRACT_16S_BIOPYTHON (
        ch_genbank.join(ASSEMBLE_CONTIGS.out.assembly_file)
    )
    // EXTRACT_16S_BIOPYTHON
    //     .view { file -> println "DEBUG: From EXTRACT_16S_BIOPYTHON, emitting file: ${file}" }
    // EXTRACT_16S_BIOPYTHON
    //     .view { item -> println "DEBUG: From EXTRACT_16S_BIOPYTHON, channel item: ${item}" }
    ch_versions = ch_versions.mix(EXTRACT_16S_BIOPYTHON.out.versions)

    // PROCESS: Extract 16S rRNA gene sequences with Barrnap if missing from 16S_EXTRACT_BIOPYTHON
    EXTRACT_16S_BARRNAP (
        ASSEMBLE_CONTIGS.out.assembly_file
            .join(EXTRACT_16S_BIOPYTHON.out.extracted_rna)
    )
    // EXTRACT_16S_BARRNAP
    //     .view { file -> println "DEBUG: From EXTRACT_16S_BARRNAP, emitting file: ${file}" }
    // EXTRACT_16S_BARRNAP
    //     .view { item -> println "DEBUG: From EXTRACT_16S_BARRNAP, channel item: ${item}" }
    ch_versions      = ch_versions.mix(EXTRACT_16S_BARRNAP.out.versions)
    ch_qc_filecheck  = ch_qc_filecheck.concat(EXTRACT_16S_BARRNAP.out.qc_filecheck)
    ch_extracted_rna = qcfilecheck(
                            "EXTRACT_16S_BARRNAP",
                            EXTRACT_16S_BARRNAP.out.qc_filecheck,
                            EXTRACT_16S_BARRNAP.out.extracted_rna
                        )

    // Prepare BLAST database for use
    if ( ch_blast_db_file ) {
        if ( ch_blast_db_file.extension in ['gz', 'tgz'] ) {
            // Add meta information
            ch_blast_db = Channel.of(ch_blast_db_file)
                            .map{
                                db ->
                                    def meta = [:]
                                    meta['id'] = db.getSimpleName()
                                    [ meta, db ]
                            }

            // Expects to be .tar.gz!
            BLAST_DB_PREPARATION_UNIX (
                ch_blast_db
            )
            ch_versions     = ch_versions.mix(BLAST_DB_PREPARATION_UNIX.out.versions)
            ch_db_for_blast = BLAST_DB_PREPARATION_UNIX.out.db.collect()

        } else if ( ch_blast_db_file.isDirectory() ) {
            ch_db_for_blast = Channel
                                    .fromPath( "${ch_blast_db_file}/{16S_ribosomal_RNA.n*,taxdb.b*}" )
                                    .collect()
                                    .map{
                                        file ->
                                            if (file.size() >= 3) {
                                                [ file[0].getSimpleName(), file ]
                                            } else {
                                                error("16S_ribosomal_RNA BLAST database requires at least '16S_ribosomal_RNA.{nin,nsq,nhr}' files.")
                                            }
                                    }
                                    .collect()

        } else {
            error("Unsupported object given to --blast_db, database must be supplied as either a directory or a .tar.gz file!")
        }

    } else {
        error("Missing 16S ribosomal RNA database! Database must be supplied to `--blast_db` as either a directory or a .tar.gz file!")
    }

    // PROCESS: Run Blast on predicted 16S ribosomal RNA genes
    ALIGN_16S_BLAST (
        ch_extracted_rna,
        ch_db_for_blast
    )
    // ALIGN_16S_BLAST
    //     .view { file -> println "DEBUG: From ALIGN_16S_BLAST, emitting file: ${file}" }
    // ALIGN_16S_BLAST
    //     .view { item -> println "DEBUG: From ALIGN_16S_BLAST, channel item: ${item}" }
    ch_versions     = ch_versions.mix(ALIGN_16S_BLAST.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(ALIGN_16S_BLAST.out.qc_filecheck)
    ch_blast_output = qcfilecheck(
                        "ALIGN_16S_BLAST",
                        ALIGN_16S_BLAST.out.qc_filecheck,
                        ALIGN_16S_BLAST.out.blast_output
                    )

    // PROCESS: Run RDP Classifier on predicted 16S ribosomal RNA genes
    CLASSIFY_16S_RDP (
        EXTRACT_16S_BARRNAP.out.extracted_rna
    )
    // CLASSIFY_16S_RDP
    //     .view { file -> println "DEBUG: From CLASSIFY_16S_RDP, emitting file: ${file}" }
    // CLASSIFY_16S_RDP
    //     .view { item -> println "DEBUG: From CLASSIFY_16S_RDP, channel item: ${item}" }
    ch_versions = ch_versions.mix(CLASSIFY_16S_RDP.out.versions)

    ch_rdp_summary = qcfilecheck(
                        "CLASSIFY_16S_RDP",
                        CLASSIFY_16S_RDP.out.qc_filecheck,
                        CLASSIFY_16S_RDP.out.rdp_tsv
                    )

    // PROCESS: Concatenate RDP summaries
    ch_rdp_summary = ch_rdp_summary
                        .map{ meta, file ->       // Map to only include the files
                            if (file.exists() && file.size() > 0) {
                                return file
                            } else {
                                error "File does not exist or empty size: ${file}"
                            }
                        }
                        .view { file -> println "DEBUG: From ch_rdp_summary, File to be collected: ${file}" }
                        .collectFile(
                            name:       "${var_assembler_name}.16S_top_genus_RDP.tsv",
                            keepHeader: true,
                            sort:       { file -> file.text },
                            storeDir:   "${params.outdir}/SSU"
                        )
                        .collectFile(
                            name:       "Summary.16S_Genus_RDP.tsv",
                            keepHeader: true,
                            sort:       { file -> file.text },
                            storeDir:   "${params.outdir}/Summaries"
                        )
                        .view { collectedFiles -> println "DEBUG: From ch_rdp_summary, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_rdp_summary)

    // PROCESS: Filter Blast output for best alignment, based on bitscore
    BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON (
        ch_blast_output
    )
    BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.summary
        .view { file -> println "DEBUG: From BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.summary, emitting file: ${file}" }
    BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.summary
        .view { item -> println "DEBUG: From BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.summary, channel item: ${item}" }
    ch_versions  = ch_versions.mix(BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.versions)
    ch_top_blast = qcfilecheck(
                        "BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON",
                        BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.qc_filecheck,
                        BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.summary
                    )

    // Collect top BLASTn species and concatenate into one file
    ch_top_blast = ch_top_blast
                        .map{ meta, file -> file
                          if (file.exists() && file.size() > 0) {
                              return file
                          } else {
                              error "File does not exist or empty size: ${file}"
                          }
                        }
                        .collectFile(
                            name:       "${var_assembler_name}.16S_top_species_BLAST.tsv",
                            keepHeader: true,
                            sort:       { file -> file.text },
                            storeDir:   "${params.outdir}/SSU"
                        )
                        .collectFile(
                            name:       "Summary.16S_Species_BLAST.tsv",
                            keepHeader: true,
                            sort:       { file -> file.text },
                            storeDir:   "${params.outdir}/Summaries"
                        )
                        .view { collectedFiles -> println "DEBUG: From ch_top_blast, collected files: ${collectedFiles}" }
    ch_output_summary_files = ch_output_summary_files.mix(ch_top_blast)

    /*
    ================================================================================
                          Perform assessment on final assembly file
    ================================================================================
    */

    ASSEMBLY_ASSESSMENT (
        ASSEMBLE_CONTIGS.out.assembly_file,
        EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.summary,
        ch_busco_config_file,
        ch_busco_db_file,
        ch_mash_db_file,
        ch_gtdbtk_db_file,
        ch_checkm2_db_file,
        ch_cat_db_file
    )
    ch_versions             = ch_versions.mix(ASSEMBLY_ASSESSMENT.out.versions)
    ch_qc_filecheck         = ch_qc_filecheck.concat(ASSEMBLY_ASSESSMENT.out.qc_filecheck)
    ch_output_summary_files = ch_output_summary_files.mix(ASSEMBLY_ASSESSMENT.out.summary_files)

    /*
    ================================================================================
                        Collect QC information
    ================================================================================
    */

    // Collect QC file checks and concatenate into one file
    ch_qc_filecheck = ch_qc_filecheck
                          .map{ meta, file -> file }
                          .collect()
                          .flatten()
                          .collectFile(
                              name:       "Summary.QC_File_Checks.tsv",
                              keepHeader: true,
                              sort:       { file -> file.text },
                              storeDir:   "${params.outdir}/Summaries",
                          )

    // Execute when the ch_qc_filecheck channel completes
    ch_qc_filecheck.subscribe { files ->
        int failCount = 0
        int passCount = 0
        int fileCount = 0

        // Loop through the collected files and count passes/fails
        new File("${params.outdir}/pipeline_info/qc_file_checks/").eachFileMatch(~'.*\\.tsv') { file ->
            fileCount++  // Increment the file counter for each file
            file.eachLine { line ->
                if (line.contains('FAIL')) {
                    System.err.println('ERROR: ' + line + ' in file: ' + file.name)
                    System.out.println('ERROR: ' + line + ' in file: ' + file.name)
                    failCount++
                } else if (line.contains('PASS')) {
                    passCount++
                }
            }
        }

        // Print the total number of files found
        System.out.println('INFO: found ' + fileCount + ' QC test results')

        // Print the total number of fails and passes
        System.err.println('INFO: found ' + failCount + ' QC failure steps across all files')
        System.out.println('INFO: found ' + passCount + ' QC passed steps, out of total: ' + fileCount)

        // Extra warning if there are any failures
        if (failCount > 0) {
            System.err.println('ERROR: found ' + failCount + ' QC files that failed!')
            System.out.println('ERROR: found ' + failCount + ' QC files that failed!')
        }
    }

    ch_output_summary_files = ch_output_summary_files.mix(ch_qc_filecheck)

    /*
    ================================================================================
                        Convert TSV outputs to Excel XLSX
    ================================================================================
    */

    if (params.create_excel_outputs) {
        tab_colors_file = file("${projectDir}/modules/local/create_excel_run_summary_python/resources/xlsx_tab_color_key.txt")

        // Collect summary files into the list_of_files variable
        list_of_files = ch_output_summary_files
            .view { item ->
                println "DEBUG: From ch_output_summary_files, Received item: ${item.getClass().getName()} - ${item}"
            }
            .filter { item ->
                // Extract the file if it's a tuple or list
                def file = (item instanceof List) ? item[1] : item  // Assuming file is the second element in the tuple
                def fileName = file.getName()  // Convert Path to String
                if (fileName.startsWith("Summary") && file.size() > 0) {
                    println "DEBUG: From ch_output_summary_files, Valid summary file found: ${fileName} (Size: ${file.size()} bytes)"
                    return true
                } else {
                    if (file.size() == 0) {
                        println "DEBUG: From ch_output_summary_files, Skipping empty summary file: ${fileName}"
                    } else {
                        println "DEBUG: From ch_output_summary_files, Skipping non-summary file: ${fileName}"
                    }
                    return false
                }
            }
            .collect()
            // Debugging print to show all files collected
            .view { files -> println "DEBUG: Files passed to CREATE_EXCEL_RUN_SUMMARY_PYTHON: ${files}" }

        // Pass both variables to the process
        def workflow_version = workflow.manifest.version ?: 'dev'
        CREATE_EXCEL_RUN_SUMMARY_PYTHON(list_of_files, tab_colors_file, workflow_version)
        ch_versions = ch_versions.mix(CREATE_EXCEL_RUN_SUMMARY_PYTHON.out.versions)
    }

    /*
    ================================================================================
                        Collect version information
    ================================================================================
    */

    // PATTERN: Collate method for version information
    ch_versions
        .unique()
        .collectFile(
            name:     "software_versions.yml",
            storeDir: params.logpath
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
