//
// Remove background host reads via `--host_remove {both,hostile,sra-human-scrubber}`
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { REMOVE_HOST_HOSTILE            } from "../../modules/local/remove_host_hostile/main"
include { REMOVE_HOST_SRA_HUMAN_SCRUBBER } from "../../modules/local/remove_host_sra_human_scrubber/main"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert params.assembler to lowercase
def toLower(it) {
    it.toString().toLowerCase()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN HOST_REMOVAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HOST_REMOVAL {

    take:
    ch_infile_handling  // channel: [ val(meta), [raw_fastq_files (R1, R2)] ]

    main:
    ch_versions      = Channel.empty()

    // Update meta to include meta.assembler
    if (toLower(params.host_remove) == "sra-human-scrubber") {
        // sra-human-scrubber removal tool
        // PROCESS: Run sra-human-scrubber to remove background host DNA read sequences
        REMOVE_HOST_SRA_HUMAN_SCRUBBER (
            ch_infile_handling
        )
        ch_versions = ch_versions.mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.versions)
        ch_host_removed_reads = REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.sra_human_scrubber_removed
    } else if (toLower(params.host_remove) == "hostile") {
        // hostile removal tool
        // PROCESS: Run hostile to remove background host DNA read sequences
        REMOVE_HOST_HOSTILE (
            ch_infile_handling
        )
        ch_versions = ch_versions.mix(REMOVE_HOST_HOSTILE.out.versions)
        ch_host_removed_reads = REMOVE_HOST_HOSTILE.out.hostile_removed
    } else if (toLower(params.host_remove) == "both") {
        // sra-human-scrubber removal tool then hostile removal tool 
        // PROCESS: Run sra-human-scrubber first and then hostile to remove background host DNA read sequences
        REMOVE_HOST_SRA_HUMAN_SCRUBBER (
            ch_infile_handling
        )
        ch_versions = ch_versions.mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.versions)

        REMOVE_HOST_HOSTILE (
            REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.sra_human_scrubber_removed
        )
        ch_versions = ch_versions.mix(REMOVE_HOST_HOSTILE.out.versions)
        ch_host_removed_reads = REMOVE_HOST_HOSTILE.out.hostile_removed
    } else {
        // Defaulting to no host removal at all
        ch_host_removed_reads = ch_infile_handling
    }

    emit:
    host_removed_reads = ch_host_removed_reads // channel: [ val(meta), [R1 FastQ file, R2 FastQ file] ]
    versions           = ch_versions
}
