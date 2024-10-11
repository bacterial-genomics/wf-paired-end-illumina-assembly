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
include { CALCULATE_METRICS_FASTQ_SEQKIT as CALC_STATS_HOSTILE_REMOVED_FQ_SEQKIT } from "../../modules/local/calculate_metrics_fastq_seqkit/main"
include { CALCULATE_METRICS_FASTQ_SEQKIT as CALC_STATS_REPAIR_REMOVED_FQ_SEQKIT  } from "../../modules/local/calculate_metrics_fastq_seqkit/main"
include { CALCULATE_METRICS_FASTQ_SEQKIT as CALC_STATS_SCRUB_REMOVED_FQ_SEQKIT   } from "../../modules/local/calculate_metrics_fastq_seqkit/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FUNCTIONS
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
    RUN HOST_REMOVAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HOST_REMOVAL {

    take:
    ch_infile_handling       // channel: [ val(meta), [raw_fastq_files (R1, R2)] ]
    ch_sra_scrubber_db_file  // channel: [ sra-human-scrubber database file ]

    main:
    ch_versions             = Channel.empty()
    ch_qc_filechecks        = Channel.empty()
    ch_output_summary_files = Channel.empty()

    // Database handling
    if ( !ch_sra_scrubber_db_file.isEmpty() ) {
        if (ch_sra_scrubber_db_file.extension in ['gz', 'tgz']) {
            // Expects to be .tar.gz!
            PREPARE_DB_SRA_HUMAN_SCRUBBER (
                ch_sra_scrubber_db_file
            )
            ch_versions = ch_versions.mix(PREPARE_DB_SRA_HUMAN_SCRUBBER.out.versions)
            ch_db_for_sra_human_scrubber = PREPARE_DB_SRA_HUMAN_SCRUBBER.out.db.collect()

        } else {
            ch_db_for_sra_human_scrubber = Channel.fromPath(ch_sra_scrubber_db_file).collect()
        }
    } else {
        ch_db_for_sra_human_scrubber = []
    }

    if ( params.update_sra_human_scrubber_db ) {
        UPDATE_DB_SRA_HUMAN_SCRUBBER (
            Channel.of("Update_SRA_Human_Scrubber_DB")
        )
        ch_versions = ch_versions.mix(UPDATE_DB_SRA_HUMAN_SCRUBBER.out.versions)
        ch_db_for_sra_human_scrubber = UPDATE_DB_SRA_HUMAN_SCRUBBER.out.db.collect()
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

        // Collect output files
        ch_cleaned_fastq = qcfilecheck(
                                "REMOVE_HOST_SRA_HUMAN_SCRUBBER",
                                REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.qc_filecheck,
                                REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.host_removed_reads
                            )

        // Collect removal summaries and concatenate into one file
        ch_scrub_summary = REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.summary
                                    .collectFile(
                                        name:       "Summary.SRA_Human_Scrubbed.tsv",
                                        keepHeader: true,
                                        sort:       { file -> file.text },
                                        storeDir:   "${params.outdir}/Summaries"
                                    )
                                    .view { collectedFiles -> println "DEBUG: From REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.summary, collected files: ${collectedFiles}" }

        ch_output_summary_files = ch_output_summary_files.mix(ch_scrub_summary)

        // PROCESS: Calculate human removed FastQ metrics for each sample with SeqKit
        CALC_STATS_SCRUB_REMOVED_FQ_SEQKIT (
            REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.host_removed_reads,
            "SRA_Scrubbed_Reads"
        )

        ch_versions = ch_versions.mix(CALC_STATS_SCRUB_REMOVED_FQ_SEQKIT.out.versions)

        // Collect scrubbed counts into one file
        ch_scrubbed_reads_metrics = CALC_STATS_SCRUB_REMOVED_FQ_SEQKIT.out.output
                                                    .collectFile(
                                                        name:       "Summary.SRA_Scrubbed_Reads.Metrics.tsv",
                                                        keepHeader: true,
                                                        sort:       { file -> file.text },
                                                        storeDir:   "${params.outdir}/Summaries"
                                                    )

        ch_output_summary_files = ch_output_summary_files.mix(ch_scrubbed_reads_metrics)

        // sra-human-scrubber non-default "-x" removes reads instead of masks
        //   with N, so it's essential to discard broken pairs or else the
        //   assembly polishing step with bwa mem and samtools mapping fails
        //   with an error, such as:
        //    [mem_sam_pe] paired reads have different names: "SRR16343585.13690", "SRR16343585.13689"
        //    [mem_sam_pe] paired reads have different names: "SRR16343585.13682", "SRR16343585.13681"
        // PROCESS: run BBTools' repair.sh to discard broken sister reads
        //   (singletons) to aggressively remove host sequence reads
        REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR (
            ch_cleaned_fastq
        )
        ch_versions           = ch_versions.mix(REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.versions)

        // Collect output files
        ch_host_removed_reads = qcfilecheck(
                                    "REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR",
                                    REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.qc_filecheck,
                                    REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.repaired_reads
                                )

        // Collect removal summaries and concatenate into one file
        ch_repair_summary = REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.summary
                                    .collectFile(
                                        name:       "Summary.BBTools_Repair_Removal.tsv",
                                        keepHeader: true,
                                        sort:       { file -> file.text },
                                        storeDir:   "${params.outdir}/Summaries"
                                    )
                                    .view { collectedFiles -> println "DEBUG: From REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.summary, collected files: ${collectedFiles}" }

        ch_output_summary_files = ch_output_summary_files.mix(ch_repair_summary)

        // PROCESS: Calculate human removed FastQ metrics for each sample with SeqKit
        CALC_STATS_REPAIR_REMOVED_FQ_SEQKIT (
            REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.repaired_reads,
            "BBTools_Repaired_Reads"
        )

        ch_versions = ch_versions.mix(CALC_STATS_REPAIR_REMOVED_FQ_SEQKIT.out.versions)

        // Collect counts from repair.sh removed broken pairs into one file
        ch_repaired_reads_metrics = CALC_STATS_REPAIR_REMOVED_FQ_SEQKIT.out.output
                                                    .collectFile(
                                                        name:       "Summary.BBTools_Repaired_Reads.Metrics.tsv",
                                                        keepHeader: true,
                                                        sort:       { file -> file.text },
                                                        storeDir:   "${params.outdir}/Summaries"
                                                    )

        ch_output_summary_files = ch_output_summary_files.mix(ch_repaired_reads_metrics)

        // Collect QC File Checks
        ch_qc_filechecks      = ch_qc_filechecks
                                    .mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.qc_filecheck)
                                    .mix(REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.qc_filecheck)

    } else if ( toLower(params.host_remove) == "hostile" ) {
        // hostile removal tool

        // PROCESS: Run hostile to remove background host DNA read sequences
        REMOVE_HOST_HOSTILE (
            ch_infile_handling
        )
        ch_versions           = ch_versions.mix(REMOVE_HOST_HOSTILE.out.versions)

        // Collect output files
        ch_host_removed_reads = qcfilecheck(
                                    "REMOVE_HOST_HOSTILE",
                                    REMOVE_HOST_HOSTILE.out.qc_filecheck,
                                    REMOVE_HOST_HOSTILE.out.host_removed_reads
                                )

        // Collect QC File Checks
        ch_qc_filechecks      = ch_qc_filechecks.mix(REMOVE_HOST_HOSTILE.out.qc_filecheck)

        // Collect removal summaries and concatenate into one file
        ch_hostile_summary = REMOVE_HOST_HOSTILE.out.summary
                                    .collectFile(
                                        name:       "Summary.BBTools_Repair_Removal.tsv",
                                        keepHeader: true,
                                        sort:       { file -> file.text },
                                        storeDir:   "${params.outdir}/Summaries"
                                    )
                                    .view { collectedFiles -> println "DEBUG: From REMOVE_HOST_HOSTILE.out.summary, collected files: ${collectedFiles}" }

        ch_output_summary_files = ch_output_summary_files.mix(ch_hostile_summary)

        // PROCESS: Calculate human removed FastQ metrics for each sample with SeqKit
        CALC_STATS_HOSTILE_REMOVED_FQ_SEQKIT (
            REMOVE_HOST_HOSTILE.out.host_removed_reads,
            "Hostile_Removed_Reads"
        )

        ch_versions = ch_versions.mix(CALC_STATS_HOSTILE_REMOVED_FQ_SEQKIT.out.versions)

        // Collect hostile removal counts into one file
        ch_hostile_removed_reads_metrics = CALC_STATS_HOSTILE_REMOVED_FQ_SEQKIT.out.output
                                                    .collectFile(
                                                        name:       "Summary.Hostile_Removed_Reads.Metrics.tsv",
                                                        keepHeader: true,
                                                        sort:       { file -> file.text },
                                                        storeDir:   "${params.outdir}/Summaries"
                                                    )

        ch_output_summary_files = ch_output_summary_files.mix(ch_hostile_removed_reads_metrics)

    } else if ( toLower(params.host_remove) == "both" && ch_db_for_sra_human_scrubber ) {
        // sra-human-scrubber removal tool (+ repair step), then hostile removal tool

        // PROCESS: Run sra-human-scrubber first and then hostile to remove background host DNA read sequences
        REMOVE_HOST_SRA_HUMAN_SCRUBBER (
            ch_infile_handling,
            ch_db_for_sra_human_scrubber
        )
        ch_versions = ch_versions.mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.versions)

        // Collect output files
        ch_sra_host_removal = qcfilecheck(
                                "REMOVE_HOST_SRA_HUMAN_SCRUBBER",
                                REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.qc_filecheck,
                                REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.host_removed_reads
                              )

        // Collect removal summaries and concatenate into one file
        ch_scrub_summary = REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.summary
                                    .collectFile(
                                        name:       "Summary.SRA_Human_Scrubbed.tsv",
                                        keepHeader: true,
                                        sort:       { file -> file.text },
                                        storeDir:   "${params.outdir}/Summaries"
                                    )
                                    .view { collectedFiles -> println "DEBUG: From REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.summary, collected files: ${collectedFiles}" }

        ch_output_summary_files = ch_output_summary_files.mix(ch_scrub_summary)

        // PROCESS: Calculate human removed FastQ metrics for each sample with SeqKit
        CALC_STATS_SCRUB_REMOVED_FQ_SEQKIT (
            REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.host_removed_reads,
            "SRA_Scrubbed_Reads"
        )

        ch_versions = ch_versions.mix(CALC_STATS_SCRUB_REMOVED_FQ_SEQKIT.out.versions)

        // Collect scrubbed counts into one file
        ch_scrubbed_reads_metrics = CALC_STATS_SCRUB_REMOVED_FQ_SEQKIT.out.output
                                                    .collectFile(
                                                        name:       "Summary.SRA_Scrubbed_Reads.Metrics.tsv",
                                                        keepHeader: true,
                                                        sort:       { file -> file.text },
                                                        storeDir:   "${params.outdir}/Summaries"
                                                    )

        ch_output_summary_files = ch_output_summary_files.mix(ch_scrubbed_reads_metrics)

        // sra-human-scrubber non-default "-x" removes reads instead of masks
        //   with N, so it's essential to discard broken pairs or else the
        //   assembly polishing step with bwa mem and samtools mapping fails
        //   with an error for the such as:
        //    [mem_sam_pe] paired reads have different names: "SRR16343585.13690", "SRR16343585.13689"
        //    [mem_sam_pe] paired reads have different names: "SRR16343585.13682", "SRR16343585.13681"
        // PROCESS: run BBTools' repair.sh to discard broken sister reads
        //   (singletons) to aggressively remove host sequence reads
        REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR (
            ch_sra_host_removal
        )
        ch_versions      = ch_versions.mix(REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.versions)

        // Collect output files
        ch_repair_sra    = qcfilecheck(
                                "REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR",
                                REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.qc_filecheck,
                                REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.repaired_reads
                            )

        // Collect removal summaries and concatenate into one file
        ch_repair_summary = REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.summary
                                    .collectFile(
                                        name:       "Summary.BBTools_Repair_Removal.tsv",
                                        keepHeader: true,
                                        sort:       { file -> file.text },
                                        storeDir:   "${params.outdir}/Summaries"
                                    )
                                    .view { collectedFiles -> println "DEBUG: From REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.summary, collected files: ${collectedFiles}" }

        ch_output_summary_files = ch_output_summary_files.mix(ch_repair_summary)

        // PROCESS: Calculate human removed FastQ metrics for each sample with SeqKit
        CALC_STATS_REPAIR_REMOVED_FQ_SEQKIT (
            REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.repaired_reads,
            "BBTools_Repaired_Reads"
        )

        ch_versions = ch_versions.mix(CALC_STATS_REPAIR_REMOVED_FQ_SEQKIT.out.versions)

        // Collect counts from repair.sh removed broken pairs into one file
        ch_repaired_reads_metrics = CALC_STATS_REPAIR_REMOVED_FQ_SEQKIT.out.output
                                                    .collectFile(
                                                        name:       "Summary.BBTools_Repaired_Reads.Metrics.tsv",
                                                        keepHeader: true,
                                                        sort:       { file -> file.text },
                                                        storeDir:   "${params.outdir}/Summaries"
                                                    )

        ch_output_summary_files = ch_output_summary_files.mix(ch_repaired_reads_metrics)

        // hostile removal tool
        // PROCESS: Run hostile to remove background host DNA read sequences
        REMOVE_HOST_HOSTILE (
            ch_repair_sra
        )
        ch_versions = ch_versions.mix(REMOVE_HOST_HOSTILE.out.versions)

        // Collect output files
        ch_host_removed_reads = qcfilecheck(
                                    "REMOVE_HOST_HOSTILE",
                                    REMOVE_HOST_HOSTILE.out.qc_filecheck,
                                    REMOVE_HOST_HOSTILE.out.host_removed_reads
                                )

        // Collect QC File Checks
        ch_qc_filechecks      = ch_qc_filechecks
                                    .mix(REMOVE_HOST_SRA_HUMAN_SCRUBBER.out.qc_filecheck)
                                    .mix(REMOVE_BROKEN_PAIRS_BBTOOLS_REPAIR.out.qc_filecheck)
                                    .mix(REMOVE_HOST_HOSTILE.out.qc_filecheck)

        // Collect removal summaries and concatenate into one file
        ch_hostile_summary = REMOVE_HOST_HOSTILE.out.summary
                                    .collectFile(
                                        name:       "Summary.Hostile_Human_Removed.tsv",
                                        keepHeader: true,
                                        sort:       { file -> file.text },
                                        storeDir:   "${params.outdir}/Summaries"
                                    )
                                    .view { collectedFiles -> println "DEBUG: From REMOVE_HOST_HOSTILE.out.summary, collected files: ${collectedFiles}" }

        ch_output_summary_files = ch_output_summary_files.mix(ch_hostile_summary)
        // PROCESS: Calculate human removed FastQ metrics for each sample with SeqKit
        CALC_STATS_HOSTILE_REMOVED_FQ_SEQKIT (
            REMOVE_HOST_HOSTILE.out.host_removed_reads,
            "Hostile_Removed_Reads"
        )

        ch_versions = ch_versions.mix(CALC_STATS_HOSTILE_REMOVED_FQ_SEQKIT.out.versions)

        // Collect hostile removal counts into one file
        ch_hostile_removed_reads_metrics = CALC_STATS_HOSTILE_REMOVED_FQ_SEQKIT.out.output
                                                    .collectFile(
                                                        name:       "Summary.Hostile_Removed_Reads.Metrics.tsv",
                                                        keepHeader: true,
                                                        sort:       { file -> file.text },
                                                        storeDir:   "${params.outdir}/Summaries"
                                                    )

        ch_output_summary_files = ch_output_summary_files.mix(ch_hostile_removed_reads_metrics)

    } else if ( toLower(params.host_remove) == "skip" ) {
        // User-specified skip host removal
        log.warn("User specified to skip host removal!")
        ch_host_removed_reads = ch_infile_handling

    } else {
        // Defaulting to no host removal
        log.warn("Host removal is not being performed.")
        ch_host_removed_reads = ch_infile_handling
    }

    emit:
    host_removed_reads   = ch_host_removed_reads // channel: [ val(meta), [host_removed_fastq_files (R1, R2)] ]
    qc_filecheck         = ch_qc_filechecks
    versions             = ch_versions
    output_summary_files = ch_output_summary_files
}
