//
// Assess assembly file using QUAST, GTDB-Tk, CAT, CheckM2, and BUSCO.
// Parameters for each of these have to be used independent of each other.
//

import java.nio.file.Files
import java.nio.file.Paths

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { CAT_DB_PREPARATION_UNIX     } from "../../modules/local/cat_db_preparation_unix/main"
include { DOWNLOAD_CAT_DB_UNIX        } from "../../modules/local/download_cat_db_unix/main"
include { BUSCO_DB_PREPARATION_UNIX   } from "../../modules/local/busco_db_preparation_unix/main"
include { SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON } from "../../modules/local/split_multifasta_assembly_biopython/main"
include { GTDBTK_DB_PREPARATION_UNIX  } from "../../modules/local/gtdbtk_db_preparation_unix/main"
include { CHECKM2_DB_PREPARATION_UNIX } from "../../modules/local/checkm2_db_preparation_unix/main"

include { QA_ASSEMBLY_QUAST           } from "../../modules/local/qa_assembly_quast/main"
include { CALCULATE_COVERAGE_UNIX     } from "../../modules/local/calculate_coverage_unix/main"
include { CLASSIFY_CONTIGS_CAT        } from "../../modules/local/classify_contigs_cat/main"
include { ASSESS_ASSEMBLY_CHECKM2     } from "../../modules/local/assess_assembly_checkm2/main"

//
// MODULES: nf-core modules
//
include { BUSCO as QA_ASSEMBLY_BUSCO              } from "../../modules/nf-core/busco/main"
include { GTDBTK_CLASSIFYWF as QA_ASSEMBLY_GTDBTK } from "../../modules/nf-core/gtdbtk/classifywf/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    RUN ASSEMBLE_CONTIGS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY_ASSESSMENT {

    take:
    ch_assembly_file           // channel: [ val(meta), [ contigs.fasta ] ]
    ch_read_alignment_stats    // channel: [ val(meta), [ Clean_Reads-AlnStats.tsv ] ]
    ch_busco_config_file       // channel: busco_config.ini
    ch_busco_db_file           // channel: database
    ch_mash_db_file            // channel: database
    ch_gtdbtk_db_file          // channel: database
    ch_checkm2_db_file         // channel: database
    ch_cat_db_file             // channel: database

    main:
    ch_versions             = Channel.empty()
    ch_qc_filecheck         = Channel.empty()
    ch_output_summary_files = Channel.empty()

    /*
    ================================================================================
                        QUAST: Genome assembly evaluation tool
    ================================================================================
    */

    // PROCESS: Run QUAST on the polished assembly for quality assessment and
    //  report the number of cleaned basepairs used to form the assembly
    QA_ASSEMBLY_QUAST (
        ch_assembly_file
    )
    ch_versions = ch_versions.mix(QA_ASSEMBLY_QUAST.out.versions)

    // Collect assembly summaries and concatenate into one file
    ch_assembly_summary = QA_ASSEMBLY_QUAST.out.summary_assemblies
                              .collectFile(
                                  name:       "Summary.Assembly_Metrics.tsv",
                                  keepHeader: true,
                                  sort:       { file -> file.text },
                                  storeDir:   "${params.outdir}/Summaries"
                              )

    ch_output_summary_files = ch_output_summary_files.mix(ch_assembly_summary)

    /*
    ================================================================================
                        Calculate coverage of assembly
    ================================================================================
    */

    // PROCESS: Calculate genome assembly depth of coverage
    CALCULATE_COVERAGE_UNIX (
        QA_ASSEMBLY_QUAST.out.qa_summaries
            .join(ch_assembly_file)
            .join(ch_read_alignment_stats)
    )
    ch_versions = ch_versions.mix(CALCULATE_COVERAGE_UNIX.out.versions)

    // Collect genome coverage summaries and concatenate into one file
    ch_genome_cov_summary = CALCULATE_COVERAGE_UNIX.out.summary
                                .collectFile(
                                    name:       "Summary.Assembly_Depth.tsv",
                                    keepHeader: true,
                                    sort:       { file -> file.text },
                                    storeDir:   "${params.outdir}/Summaries"
                                )

    ch_output_summary_files = ch_output_summary_files.mix(ch_genome_cov_summary)

    /*
    ================================================================================
                        GTDB-Tk: taxonomic classification using a GTDB reference
    ================================================================================
    */

    // PROCESS: Classify assembly FastA file using GTDB-Tk
    if ( !params.skip_gtdbtk && ch_gtdbtk_db_file ) {
        if ( ch_gtdbtk_db_file.extension in ['gz', 'tgz'] ) {
            // Add meta information
            ch_gtdb_db = Channel.of(ch_gtdbtk_db_file)
                            .map{
                                db ->
                                    def meta = [:]
                                    meta['id'] = db.getSimpleName()
                                    [ meta, db ]
                            }

            // Expects to be .tar.gz!
            GTDBTK_DB_PREPARATION_UNIX (
                ch_gtdb_db
            )
            ch_versions      = ch_versions.mix(GTDBTK_DB_PREPARATION_UNIX.out.versions)
            ch_db_for_gtdbtk = GTDBTK_DB_PREPARATION_UNIX.out.db.collect()

        } else if ( ch_gtdbtk_db_file.isDirectory() ) {
            ch_db_for_gtdbtk = Channel
                                .fromPath( "${ch_gtdbtk_db_file}/*", type: 'dir', maxDepth: 1 )
                                .collect()
                                .map{
                                    [ it[0].getSimpleName(), it ]
                                }
                                .collect()

        } else {
            error("Unsupported object given to --gtdb_db, database must be supplied as either a directory or a .tar.gz file!")
        }
    } else {
        ch_db_for_gtdbtk = Channel.empty()
    }

    // PROCESS: Perform GTDB-Tk on assembly FastA file
    QA_ASSEMBLY_GTDBTK (
        ch_assembly_file,  // tuple val(meta)   , path("bins/*")
        ch_db_for_gtdbtk,  // tuple val(db_name), path("database/*")
        '/scratch',        // val use_pplacer_scratch_dir NOTE: current nf-core module doesn't even use this path!
        []                 // path mash_db
    )
    ch_versions             = ch_versions.mix(QA_ASSEMBLY_GTDBTK.out.versions)
    ch_output_summary_files = ch_output_summary_files.mix(QA_ASSEMBLY_GTDBTK.out.summary.map{ meta, file -> file })

    // Collect GTDB-Tk summaries and concatenate into one file
    ch_gtdbtk_summary = QA_ASSEMBLY_GTDBTK.out.summary
                            .map{ meta, file -> file }  // Map to only include the files
                            .collectFile(
                                name:       "Summary.Assemblies_Classified.tsv",
                                keepHeader: true,
                                sort:       { file -> file.text },
                                storeDir:   "${params.outdir}/Summaries"
                            )

    ch_output_summary_files = ch_output_summary_files.mix(ch_gtdbtk_summary)

    /*
    ================================================================================
                        CAT: Contig Annotation Tool
    ================================================================================
    */

    if ( ch_cat_db_file ) {
        if ( ch_cat_db_file.extension in ['gz', 'tgz'] ) {
            // Add meta information
            ch_cat_db = Channel.of(ch_cat_db_file)
                            .map{
                                db ->
                                    def meta = [:]
                                    meta['id'] = db.getSimpleName()
                                    [ meta, db ]
                            }

            // Expects to be .tar.gz!
            CAT_DB_PREPARATION_UNIX (
                ch_cat_db
            )
            ch_versions   = ch_versions.mix(CAT_DB_PREPARATION_UNIX.out.versions)
            ch_db_for_cat = CAT_DB_PREPARATION_UNIX.out.db.collect()

        } else if ( ch_cat_db_file.isDirectory() ) {
            ch_db_for_cat = Channel
                                .fromPath( "${ch_cat_db_file}/{db,tax}", type: 'dir', maxDepth: 1 )
                                .collect()

        } else {
            error("Unsupported object given to --cat_db, database must be supplied as either a directory or a .tar.gz file!")
        }

    } else if ( !ch_cat_db_file && params.download_cat_db ) {
        DOWNLOAD_CAT_DB_UNIX (
            Channel.of("Download CAT DB")
                    .map {
                        def meta = [:]
                        meta['id'] = it
                        [ meta ]
                    }
        )
        ch_versions    = ch_versions.mix(DOWNLOAD_CAT_DB_UNIX.out.versions)
        ch_cat_db_file = DOWNLOAD_CAT_DB_UNIX.out.db.collect()

        // Expects to be .tar.gz!
        CAT_DB_PREPARATION_UNIX (
            ch_cat_db_file
        )
        ch_versions   = ch_versions.mix(CAT_DB_PREPARATION_UNIX.out.versions)
        ch_db_for_cat = CAT_DB_PREPARATION_UNIX.out.db

    } else {
            ch_db_for_cat = Channel.empty()
    }

    CLASSIFY_CONTIGS_CAT (
        ch_assembly_file,
        ch_db_for_cat
    )
    ch_versions     = ch_versions.mix(CLASSIFY_CONTIGS_CAT.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.mix(CLASSIFY_CONTIGS_CAT.out.qc_filecheck)
    ch_cat_output   = qcfilecheck(
                        "CLASSIFY_CONTIGS_CAT",
                        CLASSIFY_CONTIGS_CAT.out.qc_filecheck,
                        CLASSIFY_CONTIGS_CAT.out.summary
                    )

    ch_output_summary_files = ch_output_summary_files.mix(ch_cat_output.map{ meta, file -> file })

    // Collect CAT summaries and concatenate into one file
    ch_classified_contigs_cat_summary = CLASSIFY_CONTIGS_CAT.out.output
                                            .collectFile(
                                                name:       "Summary.Contigs_Classified.tsv",
                                                keepHeader: true,
                                                sort:       { file -> file.text },
                                                storeDir:   "${params.outdir}/Summaries"
                                            )

    ch_output_summary_files = ch_output_summary_files.mix(ch_classified_contigs_cat_summary)

    /*
    ================================================================================
                        CheckM2: Check for completeness and contamination
    ================================================================================
    */

    if ( ch_checkm2_db_file ) {
        // Add meta information
        ch_checkm2_db = Channel.of(ch_checkm2_db_file)
                        .map{
                            db ->
                                def meta = [:]
                                meta['id'] = db.getSimpleName()
                                [ meta, db ]
                        }

        if ( ch_checkm2_db_file.extension in ['gz', 'tgz'] ) {
            // Expects to be .tar.gz!
            CHECKM2_DB_PREPARATION_UNIX (
                ch_checkm2_db
            )
            ch_versions       = ch_versions.mix(CHECKM2_DB_PREPARATION_UNIX.out.versions)
            ch_db_for_checkm2 = CHECKM2_DB_PREPARATION_UNIX.out.db.collect()

        } else if ( ch_checkm2_db_file.isDirectory() ) {
            // Parse directory for `.dmnd` file
            ch_db_for_checkm2 = Channel
                                .fromPath( "${ch_checkm2_db_file}/*.dmnd", maxDepth: 1 )
                                .collect()
        } else {
            error("Unsupported object given to --checkm2_db, database must be supplied as either a directory or a .tar.gz file!")
        }

    } else {
            ch_db_for_checkm2 = Channel.empty()
    }

    ASSESS_ASSEMBLY_CHECKM2 (
        ch_assembly_file,
        ch_db_for_checkm2
    )
    ch_versions        = ch_versions.mix(ASSESS_ASSEMBLY_CHECKM2.out.versions)
    ch_qc_filecheck    = ch_qc_filecheck.mix(ASSESS_ASSEMBLY_CHECKM2.out.qc_filecheck)
    ch_checkm2_output  = qcfilecheck(
                            "ASSESS_ASSEMBLY_CHECKM2",
                            ASSESS_ASSEMBLY_CHECKM2.out.qc_filecheck,
                            ASSESS_ASSEMBLY_CHECKM2.out.summary
                        )

    // Concatenate CheckM2 summaries
    ch_checkm2_output = ASSESS_ASSEMBLY_CHECKM2.out.summary
                            .map{ meta, file -> file }  // Map to only include the files
                            .collectFile(
                                name:       "Summary.Assembly_Completeness.tsv",
                                keepHeader: true,
                                sort:       { file -> file.text },
                                storeDir:   "${params.outdir}/Summaries"
    )

    ch_output_summary_files = ch_output_summary_files.mix(ch_checkm2_output)

    /*
    ================================================================================
                        BUSCO: predict genes on contigs
    ================================================================================
    */

    // PROCESS: Classify contigs with BUSCO
    if ( !params.skip_busco && ch_busco_db_file ) {
        if ( ch_busco_db_file.extension in ['gz', 'tgz'] ) {
            // Add meta information
            ch_busco_db = Channel.of(ch_busco_db_file)
                            .map{
                                db ->
                                    def meta = [:]
                                    meta['id'] = db.getSimpleName()
                                    [ meta, db ]
                            }

            // Expects to be tar.gz!
            BUSCO_DB_PREPARATION_UNIX (
                ch_busco_db
            )
            ch_versions     = ch_versions.mix(BUSCO_DB_PREPARATION_UNIX.out.versions)
            ch_db_for_busco = BUSCO_DB_PREPARATION_UNIX.out.db.collect()

        } else if ( ch_busco_db_file.isDirectory() ) {
            // Expects directory in <database>/lineages/<lineage>_odb10 format!
            ch_db_for_busco = Channel
                                .fromPath(ch_busco_db_file)
                                .map{
                                    db ->
                                        if ( db.getSimpleName().contains('odb10') ) {
                                            if ( db.getParent().getSimpleName() == "lineages" ) {
                                                db.getParent().getParent()
                                            } else {
                                                error("Unsupported object given to --busco_db, database directory must be in format `<database>/lineages/<lineage>_odb10`!")
                                            }
                                        } else {
                                            db
                                        }
                                }
                                .collect()
            // Debug prints
            println "DEBUG: Type of ch_db_for_busco: ${ch_db_for_busco.getClass()}"
            println "DEBUG: Contents of ch_db_for_busco: ${ch_db_for_busco.inspect()}"
        } else {
            error("Unsupported object given to --busco_db, database must be supplied as either a directory or a .tar.gz file!")
        }

        // NOTE: this defaults to "[auto]" not "auto" and gives error, no output
        // "ERROR: [auto]_odb10 is not a valid option for 'lineages'" without [0] for channel item in QA_ASSEMBLY_BUSCO
        ch_lineage_for_busco_db = Channel
                                    .of(ch_busco_db_file)
                                    .map{
                                        db ->
                                            db = db.getSimpleName()
                                            db.contains('odb10') ? db : 'auto'
                                    }
                                    .collect()

        // PROCESS: Split assembly FastA file into individual contig files
        SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON (
            ch_assembly_file
        )
        ch_versions = ch_versions.mix(SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON.out.versions)

        // PROCESS: Perform BUSCO analysis on contigs
            // SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON.out.split_multifasta_assembly_dir, // tuple val(meta), path('tmp_input/*')
            // ch_lineage_for_busco_db, // val lineage ; Required: lineage to check against, "auto" enables --auto-lineage instead
        QA_ASSEMBLY_BUSCO (
            ch_assembly_file,        // tuple val(meta), path('tmp_input/*')
            'genome',                // val mode ; Required: One of genome, proteins, or transcriptome
            'auto',                  // val lineage ; Required: lineage to check against, "auto" enables --auto-lineage instead
            ch_db_for_busco,         // path busco_lineages_path ; Recommended: path to busco lineages - downloads if not set
            ch_busco_config_file     // path config_file ; Optional: busco configuration file
        )
        ch_versions             = ch_versions.mix(QA_ASSEMBLY_BUSCO.out.versions)
        ch_output_summary_files = ch_output_summary_files.mix(QA_ASSEMBLY_BUSCO.out.batch_summary.map{ meta, file -> file })

        // Collect BUSCO summaries and concatenate into one file
        println "DEBUG: ch_busco_summary QA_ASSEMBLY_BUSCO.out.batch_summary = ${QA_ASSEMBLY_BUSCO.out.batch_summary}"
        ch_busco_summary = QA_ASSEMBLY_BUSCO.out.batch_summary
            .map { meta, file ->
                println "DEBUG: ch_busco_summary meta.id = ${meta.id} (${meta.getClass()})"
                println "DEBUG: ch_busco_summary file = ${file} (${file.getClass()})"
                return [meta.id, file] // meta.id and file path (or string representing file path)
            }
            .collectFile(
                name:       "Summary.BUSCO_Completeness.tsv",
                keepHeader: true,
                sort: {
                    // Adding debug to print all incoming elements
                    println "DEBUG: ch_busco_summary Sorting pair: ${it} (${it.getClass()})"
                    def filePath = (it[1] instanceof String) ? Paths.get(it[1]) : it[1]
                    println "DEBUG: ch_busco_summary Resolved filePath = ${filePath}"

                    // Check if the file exists before reading
                    if (Files.exists(filePath)) {
                        Files.readString(filePath)
                    } else {
                        println "ERROR: ch_busco_summary File does not exist at path: ${filePath}"
                        return ""
                    }
                },
                storeDir: "${params.outdir}/Summaries"
            ) { pair ->
                def sample_id = pair[0]
                def path = pair[1]

                println "DEBUG: ch_busco_summary Processing sample_id = ${sample_id}"
                println "DEBUG: ch_busco_summary Processing path = ${path} (${path.getClass()})"

                def filePath = (path instanceof String) ? Paths.get(path) : path

                if (!Files.exists(filePath)) {
                    println "ERROR: ch_busco_summary File does not exist for sample ${sample_id}: ${filePath}"
                    return ""
                }

                def lines = Files.readAllLines(filePath) // Read each line as a list of strings
                lines[0] = "Sample_name\t" + lines[0].replaceAll(/\s/, '_')

                def modifiedRows = lines.collect { line ->
                    def columns = line.split('\t').toList()
                    columns.remove(1) // Remove the second column
                    return columns
                }

                def idx = 0
                def finalRows = modifiedRows.collect { columns ->
                    def result
                    if (idx == 0) {
                        result = "Sample_name\t" + columns.join('\t')
                    } else {
                        result = "$sample_id\t" + columns.join('\t')
                    }
                    idx++
                    return result
                }

                return finalRows.join('\n') + '\n'
            }

        // ch_output_summary_files = ch_output_summary_files.mix(ch_busco_summary)
    }

    emit:
    versions      = ch_versions
    qc_filecheck  = ch_qc_filecheck
    summary_files = ch_output_summary_files
}
