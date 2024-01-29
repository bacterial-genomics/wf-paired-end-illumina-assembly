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
include { REMOVE_HOST_HOSTILE                } from "../../modules/local/remove_host_hostile/main"
include { REMOVE_HOST_SRA_HUMAN_SCRUBBER     } from "../../modules/local/remove_host_sra_human_scrubber/main"
include { UPDATE_DB_SRA_HUMAN_SCRUBBER       } from "../../modules/local/update_db_sra_human_scrubber/main"
include { REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR } from "../../modules/local/remove_broken_pairs_bbtools_repair/main"

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
    ch_infile_handling       // channel: [ val(meta), [raw_fastq_files (R1, R2)] ]
    ch_sra_scrubber_db_file  // channel: [ sra-human-scrubber database file ]

    main:
    ch_versions      = Channel.empty()
    ch_qc_filechecks = Channel.empty()

    // Database handling
    if ( !ch_sra_scrubber_db_file.isEmpty() ) {
        if (ch_sra_scrubber_db_file.extension in ['gz', 'tgz']) {
            // Expects to be .tar.gz!
            PREPARE_DB_SRA_HUMAN_SCRUBBER (
                ch_sra_scrubber_db_file
            )
            ch_versions = ch_versions.mix(PREPARE_DB_SRA_HUMAN_SCRUBBER.out.versions)
            ch_db_for_sra_human_scrubber = PREPARE_DB_SRA_HUMAN_SCRUBBER.out.db

        } else {
            ch_db_for_sra_human_scrubber = Channel.fromPath(ch_sra_scrubber_db_file)
        }
    } else {
        ch_db_for_sra_human_scrubber = []
    }

    if ( params.update_sra_human_scrubber_db ) {
        UPDATE_DB_SRA_HUMAN_SCRUBBER (
            Channel.of("Update_SRA_Human_Scrubber_DB")
        )
        ch_versions = ch_versions.mix(UPDATE_DB_SRA_HUMAN_SCRUBBER.out.versions)
        ch_db_for_sra_human_scrubber = UPDATE_DB_SRA_HUMAN_SCRUBBER.out.db
    }

    // Perform host removal
    if ( toLower(params.host_remove) == "sra-human-scrubber" && ch_db_for_sra_human_scrubber ) {
        // sra-human-scrubber removal tool
        // PROCESS: Run sra-human-scrubber to remove background host DNA read sequences
        REMOVE_HOST_SRA_HUMAN_SCRUBBER (
            ch_infile_handling,
            ch_db_for_sra_human_scrubber
        )
        ch_versions      = ch_versions.mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.versions)
        ch_qc_filechecks = ch_qc_filechecks.mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.qc_filecheck)

        // sra-human-scrubber non-default "-x" removes reads instead of masks
        //   with N, so it's essential to discard broken pairs or else the
        //   assembly polishing step with bwa mem and samtools mapping fails
        //   with an error for the such as:
        //    [mem_sam_pe] paired reads have different names: "SRR16343585.13690", "SRR16343585.13689"
        //    [mem_sam_pe] paired reads have different names: "SRR16343585.13682", "SRR16343585.13681"
        // PROCESS: run BBTools' repair.sh to discard broken sister reads
        //   (singletons) to aggressively remove host sequence reads
        REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR (
            REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.sra_human_scrubber_removed
        )
        ch_versions           = ch_versions.mix(REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.versions)
        ch_qc_filechecks      = ch_qc_filechecks.mix(REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.qc_filecheck)
        ch_host_removed_reads = REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.fastq_removed_broken_pairs

    } else if ( toLower(params.host_remove) == "hostile" ) {
        // hostile removal tool
        // PROCESS: Run hostile to remove background host DNA read sequences
        REMOVE_HOST_HOSTILE (
            ch_infile_handling
        )
        ch_versions           = ch_versions.mix(REMOVE_HOST_HOSTILE.out.versions)
        ch_qc_filechecks      = ch_qc_filechecks.mix(REMOVE_HOST_HOSTILE.out.qc_filecheck)
        ch_host_removed_reads = REMOVE_HOST_HOSTILE.out.hostile_removed

    } else if ( toLower(params.host_remove) == "both" && ch_db_for_sra_human_scrubber ) {
        // sra-human-scrubber removal tool then hostile removal tool
        // PROCESS: Run sra-human-scrubber first and then hostile to remove background host DNA read sequences
        REMOVE_HOST_SRA_HUMAN_SCRUBBER (
            ch_infile_handling,
            ch_db_for_sra_human_scrubber
        )
        ch_versions      = ch_versions.mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.versions)
        ch_qc_filechecks = ch_qc_filechecks.mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.qc_filecheck)

        // sra-human-scrubber non-default "-x" removes reads instead of masks
        //   with N, so it's essential to discard broken pairs or else the
        //   assembly polishing step with bwa mem and samtools mapping fails
        //   with an error for the such as:
        //    [mem_sam_pe] paired reads have different names: "SRR16343585.13690", "SRR16343585.13689"
        //    [mem_sam_pe] paired reads have different names: "SRR16343585.13682", "SRR16343585.13681"
        // PROCESS: run BBTools' repair.sh to discard broken sister reads
        //   (singletons) to aggressively remove host sequence reads
        REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR (
            REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.sra_human_scrubber_removed
        )
        ch_versions      = ch_versions.mix(REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.versions)
        ch_qc_filechecks = ch_qc_filechecks.mix(REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.qc_filecheck)

        // hostile removal tool
        // PROCESS: Run hostile to remove background host DNA read sequences
        REMOVE_HOST_HOSTILE (
            REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.fastq_removed_broken_pairs
        )
        ch_versions           = ch_versions.mix(REMOVE_HOST_HOSTILE.out.versions)
        ch_qc_filechecks      = ch_qc_filechecks.mix(REMOVE_HOST_HOSTILE.out.qc_filecheck)
        ch_host_removed_reads = REMOVE_HOST_HOSTILE.out.hostile_removed

    } else if ( toLower(params.host_remove) == "skip" ) {
        // User-specified skip host removal
        log.warn("Host removal user-specified to skip.")
        ch_host_removed_reads = ch_infile_handling

    } else {
        // Defaulting to no host removal
        log.warn("Host removal is not being performed.")
        ch_host_removed_reads = ch_infile_handling
    }

    emit:
    host_removed_reads = ch_host_removed_reads // channel: [ val(meta), [host_removed_fastq_files (R1, R2)] ]
    qc_filecheck       = ch_qc_filechecks
    versions           = ch_versions
}
