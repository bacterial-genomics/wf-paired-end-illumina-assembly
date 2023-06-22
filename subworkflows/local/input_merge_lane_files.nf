//
// Merge FastQ files with multiple lanes
// Adapted from https://github.com/nf-core/mag
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { CONVERT_SAMPLESHEET_PYTHON } from "../../modules/local/convert_samplesheet_python/main"
include { MERGE_LANE_FILES_PYTHON } from "../../modules/local/merge_lane_files_python/main"

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN INPUT_MERGE_LANE_FILES WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_MERGE_LANE_FILES {

    take:
    input

    main:
    ch_versions = Channel.empty()

    if (params.merge_lanes) {
        // If merged_lanes parameter is used, merge
        // multiple lanes based on 'sample' column
        MERGE_LANE_FILES_PYTHON(input)
        merged_output = MERGE_LANE_FILES_PYTHON.out.lanes_merged_samplesheet
        ch_versions = ch_versions.mix(MERGE_LANE_FILES_PYTHON.out.versions)
    } else {
        // If merged_lanes parameter is not used,
        // output the input channel
        merged_output = input
    }

    emit:
    output = merged_output
    versions = ch_versions
}
