//
// Perform downsampling on input FastQ files
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { ESTIMATE_GENOME_SIZE_KMC           } from "../../modules/local/estimate_genome_size_kmc/main"
include { COUNT_TOTAL_BP_INPUT_READS_SEQTK   } from "../../modules/local/count_total_bp_input_reads_seqtk/main"
include { COUNT_TOTAL_BP_INPUT_READS_SEQKIT  } from "../../modules/local/count_total_bp_input_reads_seqkit/main"
include { ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX } from "../../modules/local/estimate_original_input_depth_unix/main"
include { SUBSAMPLE_READS_TO_DEPTH_SEQTK     } from "../../modules/local/subsample_reads_to_depth_seqtk/main"
include { SUBSAMPLE_READS_TO_DEPTH_SEQKIT    } from "../../modules/local/subsample_reads_to_depth_seqkit/main"
include { CALCULATE_METRICS_FASTQ_SEQTK  as CALC_STATS_DOWNSAMPLE_FQ_SEQTK  } from "../../modules/local/calculate_metrics_fastq_seqtk/main"
include { CALCULATE_METRICS_FASTQ_SEQKIT as CALC_STATS_DOWNSAMPLE_FQ_SEQKIT } from "../../modules/local/calculate_metrics_fastq_seqkit/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert params.assembler to lowercase
def toLower(it) {
    it.toString().toLowerCase()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DOWNSAMPLE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOWNSAMPLE {

    take:
    ch_raw_reads    // channel: [ val(meta), [reads] ]

    main:
    ch_versions             = Channel.empty()
    ch_output_summary_files = Channel.empty()

    // Handle too much raw data, subsample the input FastQ files
    // Only calculate genome size if genome_size is unknown or a depth given to subsample
    if (!params.genome_size && params.depth > 0) {
        // Estimate the genome size from the input R1 FastQ file
        ESTIMATE_GENOME_SIZE_KMC (
            ch_raw_reads
        )

        ch_versions = ch_versions.mix(ESTIMATE_GENOME_SIZE_KMC.out.versions)

    } else {
        if (params.depth <= 0) {
            log.info("Depth is set to <= 0x. No subsampling to perform and therefore no genome size estimation required.")
        } else {
            log.info("Using the user-input genome size of ${params.genome_size}bp")
        }
        // pass the genome_size val onto the next depth channel
        // Skip genome size estimation based on user input, but
        //  still consider downsampling with specified genome_size input value
    }

    // Only if specified depth is less than wanted depth, subsample infiles
    if (params.depth > 0) {
        log.info("Estimating if the input exceeds ${params.depth}x")

        // Subsample with seqtk
        if ( toLower(params.subsample_tool) == "seqtk" ) {
            
            // Use the genome size to figure out the expected depth
            COUNT_TOTAL_BP_INPUT_READS_SEQTK (
                ch_raw_reads
            )

            ch_versions = ch_versions.mix(COUNT_TOTAL_BP_INPUT_READS_SEQTK.out.versions)

            ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX (
                COUNT_TOTAL_BP_INPUT_READS_SEQTK.out.input_total_bp
                        .join(ESTIMATE_GENOME_SIZE_KMC.out.genome_size)
            )

            ch_versions = ch_versions.mix(ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX.out.versions)

            SUBSAMPLE_READS_TO_DEPTH_SEQTK (
                ch_raw_reads.join(ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX.out.fraction_of_reads_to_use)
            )

            ch_versions = ch_versions.mix(SUBSAMPLE_READS_TO_DEPTH_SEQTK.out.versions)

            // Collect subsampled reads
            ch_downsampled_reads = SUBSAMPLE_READS_TO_DEPTH_SEQTK.out.reads

            // PROCESS: Calculate downsampled FastQ metrics for each sample with Seqtk
            CALC_STATS_DOWNSAMPLE_FQ_SEQTK (
                SUBSAMPLE_READS_TO_DEPTH_SEQTK.out.reads,
                "Downsampled_Reads"
            )

            ch_versions = ch_versions.mix(CALC_STATS_DOWNSAMPLE_FQ_SEQTK.out.versions)

            // Collect cleaned read/base summaries and concatenate into one file
            ch_downsampled_reads_metrics_summary = CALC_STATS_DOWNSAMPLE_FQ_SEQTK.out.output
                                                        .collectFile(
                                                            name:       "Summary.Downsampled_Reads.Metrics.tsv",
                                                            keepHeader: true,
                                                            sort:       { file -> file.text },
                                                            storeDir:   "${params.outdir}/Summaries"
                                                        )

            ch_output_summary_files = ch_output_summary_files.mix(ch_downsampled_reads_metrics_summary)

        // Subsample with SeqKit
        } else if ( toLower(params.subsample_tool) == "seqkit" ) {

            // // Use the genome size to figure out the expected depth
            // CALC_STATS_INPUT_FQ_SEQKIT (
            //     ch_raw_reads,
            //     "Input_for_Subsampling_Reads"
            // )
            // ch_versions = ch_versions.mix(CALC_STATS_INPUT_FQ_SEQKIT.out.versions)

            COUNT_TOTAL_BP_INPUT_READS_SEQKIT (
                ch_raw_reads
            )

            ch_versions = ch_versions.mix(COUNT_TOTAL_BP_INPUT_READS_SEQKIT.out.versions)

            ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX (
                COUNT_TOTAL_BP_INPUT_READS_SEQKIT.out.input_total_bp
                        .join(ESTIMATE_GENOME_SIZE_KMC.out.genome_size)
            )

            ch_versions = ch_versions.mix(ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX.out.versions)

            // Subsample with seqkit
            SUBSAMPLE_READS_TO_DEPTH_SEQKIT (
                ch_raw_reads.join(ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX.out.fraction_of_reads_to_use)
            )

            ch_versions = ch_versions.mix(SUBSAMPLE_READS_TO_DEPTH_SEQKIT.out.versions)

            // Collect subsampled reads
            ch_downsampled_reads = SUBSAMPLE_READS_TO_DEPTH_SEQKIT.out.reads

            // PROCESS: Calculate downsampled FastQ metrics for each sample with SeqKit
            CALC_STATS_DOWNSAMPLE_FQ_SEQKIT (
                SUBSAMPLE_READS_TO_DEPTH_SEQKIT.out.reads,
                "Downsampled_Reads"
            )

            ch_versions = ch_versions.mix(CALC_STATS_DOWNSAMPLE_FQ_SEQKIT.out.versions)

            // Collect cleaned read/base summaries and concatenate into one file
            ch_downsampled_reads_metrics_summary = CALC_STATS_DOWNSAMPLE_FQ_SEQKIT.out.output
                                                        .collectFile(
                                                            name:       "Summary.Downsampled_Reads.Metrics.tsv",
                                                            keepHeader: true,
                                                            sort:       { file -> file.text },
                                                            storeDir:   "${params.outdir}/Summaries"
                                                        )

            ch_output_summary_files = ch_output_summary_files.mix(ch_downsampled_reads_metrics_summary)
        }
    } else {
        // Skip subsampling and pass raw reads to PhiX removal
        // Collect raw reads
        ch_downsampled_reads = ch_raw_reads
    }

    emit:
    reads                = ch_downsampled_reads    // channel: [ val(meta), [reads] ]
    versions             = ch_versions
    output_summary_files = ch_output_summary_files
}
