//
// Check input samplesheet and get read channels
// Adapted from https://github.com/nf-core/mag
//

include { CONVERT_SAMPLESHEET_PYTHON } from "../../modules/local/convert_samplesheet_python/main"
include { MERGE_LANE_FILES_PYTHON } from "../../modules/local/merge_lane_files_python/main"

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Subworkflow of INPUT_CHECK
// Can't run a process within a groovy function
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

// Main workflow
workflow INPUT_CHECK {
    take:
    ch_input

    main:
    ch_versions = Channel.empty()

    if (hasExtension(ch_input, "csv")) {
        // Merge lanes if merged_lanes parameter is used
        INPUT_MERGE_LANE_FILES(ch_input)

        // Extracts read files from samplesheet CSV and distribute into channels
        ch_input_rows = INPUT_MERGE_LANE_FILES.out.output
            .splitCsv(header: true)
            .map { row ->
                    if (row.size() == 3) {
                        def id = row.sample
                        def sr1 = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : false
                        def sr2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : false
                        // Check if given combination is valid
                        if (!id) exit 1, "Invalid input samplesheet: sample can not be empty."
                        if (!sr1) exit 1, "Invalid input samplesheet: fastq_1 can not be empty."
                        if (!sr2) exit 1, "Invalid input samplesheet: fastq_2 can not be empty."
                        return [ id, sr1, sr2 ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 3."
                    }
                }
        // Separate reads
        ch_raw_reads = ch_input_rows
            .map { id, sr1, sr2 ->
                        def meta = [:]
                        meta.id  = id
                        return [ meta, [ sr1, sr2 ] ]
                }

        // Collect version info
        ch_versions = ch_versions.mix(INPUT_MERGE_LANE_FILES.out.versions)
    } else if (hasExtension(ch_input, "tsv")) {
        // Merge lanes if merged_lanes parameter is used
        INPUT_MERGE_LANE_FILES(ch_input)

        // Extracts read files from samplesheet TSV and distribute into channels
        ch_input_rows = INPUT_MERGE_LANE_FILES.out.output
            .splitCsv(header: true, sep:'\t')
            .map { row ->
                    if (row.size() == 3) {
                        def id = row.sample
                        def sr1 = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : false
                        def sr2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : false
                        // Check if given combination is valid
                        if (!id) exit 1, "Invalid input samplesheet: sample can not be empty."
                        if (!sr1) exit 1, "Invalid input samplesheet: fastq_1 can not be empty."
                        if (!sr2) exit 1, "Invalid input samplesheet: fastq_2 can not be empty."
                        return [ id, sr1, sr2 ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 3."
                    }
                }
        // Separate reads
        ch_raw_reads = ch_input_rows
            .map { id, sr1, sr2 ->
                        def meta = [:]
                        meta.id  = id
                        return [ meta, [ sr1, sr2 ] ]
                }

        // Collect version info
        ch_versions = ch_versions.mix(INPUT_MERGE_LANE_FILES.out.versions)
    } else if (hasExtension(ch_input, "xlsx") || hasExtension(ch_input, "xls") || hasExtension(ch_input, "ods")) {
        // Convert samplesehet to TSV format
        CONVERT_SAMPLESHEET_PYTHON(ch_input)

        // Merge lanes if merged_lanes parameter is used
        INPUT_MERGE_LANE_FILES(CONVERT_SAMPLESHEET_PYTHON.out.converted_samplesheet)

        // Extracts read files from TSV samplesheet created
        // in above process and distribute into channels
        ch_input_rows = INPUT_MERGE_LANE_FILES.out.output
            .splitCsv(header: true, sep:'\t')
            .map { row ->
                    if (row.size() == 3) {
                        def id = row.sample
                        def sr1 = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : false
                        def sr2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : false
                        // Check if given combination is valid
                        if (!id) exit 1, "Invalid input samplesheet: sample can not be empty."
                        if (!sr1) exit 1, "Invalid input samplesheet: fastq_1 can not be empty."
                        if (!sr2) exit 1, "Invalid input samplesheet: fastq_2 can not be empty."
                        return [ id, sr1, sr2 ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 3."
                    }
                }
        // Separate reads
        ch_raw_reads = ch_input_rows
            .map { id, sr1, sr2 ->
                        def meta = [:]
                        meta.id  = id
                        return [ meta, [ sr1, sr2 ] ]
                }

        // Collect version info
        ch_versions = ch_versions
                        .mix(INPUT_MERGE_LANE_FILES.out.versions)
                        .mix(CONVERT_SAMPLESHEET_PYTHON.out.versions)
    } else {
        // Read from FilePairs if no samplesheet is given
        ch_raw_reads = Channel
            .fromFilePairs(ch_input+'/**_{,R}{1,2}*{fastq,fq}{,.gz}', maxDepth: 2, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${input}\nNB: Path needs to be enclosed in quotes!" }
            .map { row ->
                        def meta = [:]
                        meta.id  = row[0]
                        return [ meta, row[1] ]
                }
        ch_input_rows = Channel.empty()
    }

    // Ensure sample IDs are unique
    ch_input_rows
        .map { id, sr1, sr2 -> id }
        .toList()
        .map { ids -> if( ids.size() != ids.unique().size() ) {exit 1, "ERROR: input samplesheet contains duplicated sample IDs!" } }

    emit:
    raw_reads = ch_raw_reads
    versions = ch_versions
}