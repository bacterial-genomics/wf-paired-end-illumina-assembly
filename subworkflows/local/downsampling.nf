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
include { ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX } from "../../modules/local/estimate_original_input_depth_unix/main"
include { SUBSAMPLE_READS_TO_DEPTH_SEQTK     } from "../../modules/local/subsample_reads_to_depth_seqtk/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DOWNSAMPLE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOWNSAMPLE {

    take:
    ch_raw_reads    // channel: [ val(meta), [reads]. [qc_filechecks] ]

    main:
    ch_versions = Channel.empty()

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
            println "Depth is set to ${params.genome_size}x. No subsampling to perform and therefore no genome size estimation required."
        } else {
            println "Using the user-input genome size of ${params.genome_size}bp"
        }
        // pass the genome_size val onto the next depth channel
        // Skip genome size estimation based on user input, but
        //  still consider downsampling with specified genome_size input value
    }

    if (params.depth > 0) {
        println "Estimating if the input exceeds ${params.depth}x"

        // Use the genome size to figure out the expected depth
        COUNT_TOTAL_BP_INPUT_READS_SEQTK (
            ch_raw_reads
        )

        ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX (
            COUNT_TOTAL_BP_INPUT_READS_SEQTK.out.total_bp
                    .join(ESTIMATE_GENOME_SIZE_KMC.out.genome_size)
        )

        // Only if specified depth is less than wanted depth, subsample infiles
        SUBSAMPLE_READS_TO_DEPTH_SEQTK (
            ch_raw_reads.join(ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX.out.fraction_of_reads_to_use)
        )

        // Collect subsampled reads
        ch_downsampled_reads = SUBSAMPLE_READS_TO_DEPTH_SEQTK.out.reads

        // Collect version info
        ch_versions = ch_versions
                        .mix(COUNT_TOTAL_BP_INPUT_READS_SEQTK.out.versions)
                        .mix(ESTIMATE_ORIGINAL_INPUT_DEPTH_UNIX.out.versions)
                        .mix(SUBSAMPLE_READS_TO_DEPTH_SEQTK.out.versions)
    } else {
        // Skip subsampling and pass raw reads to PhiX removal
        // Collect raw reads
        ch_downsampled_reads = ch_raw_reads
    }

    emit:
    reads     = ch_downsampled_reads    // channel: [ val(meta), [reads]. [qc_filechecks] ]
    versions  = ch_versions
}
