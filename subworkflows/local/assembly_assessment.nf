//
// Assess assembly file using QUAST, GTDB-Tk, CAT, CheckM2, and BUSCO.
// Parameters for each of these have to be used independent of each other.
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { QA_ASSEMBLY_QUAST                       } from "../modules/local/qa_assembly_quast/main"
include { CLASSIFY_CONTIGS_CAT                    } from "../modules/local/classify_contigs_cat/main"
include { CLASSIFY_ASSEMBLY_CHECKM2               } from "../modules/local/qa_assembly_checkm2/main"
include { BUSCO_DB_PREPARATION_UNIX               } from "../modules/local/busco_db_preparation_unix/main"
include { GTDBTK_DB_PREPARATION_UNIX              } from "../modules/local/gtdbtk_db_preparation_unix/main"

//
// MODULES: nf-core modules
//
include { BUSCO as QA_ASSEMBLY_BUSCO              } from "../modules/nf-core/busco/main"
include { GTDBTK_CLASSIFYWF as QA_ASSEMBLY_GTDBTK } from "../modules/nf-core/gtdbtk/classifywf/main"

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
    ch_cleaned_fastq_files     // channel: [ val(meta), [ cleaned_fastq_files (R1, R2, single) ] ]
    ch_read_alignment_stats    // channel: [ val(meta), [ CleanedReads-AlnStats.tsv ] ]
    ch_busco_config_file       // channel: [ busco_config.ini ]
    ch_busco_db_file           // channel: [ database ]
    ch_gtdbtk_db_file          // channel: [ database ]
    ch_checkm2_db_file         // channel: [ database ]
    ch_cat_taxonomy_db_file    // channel: [ database ]
    ch_cat_alignment_db_file   // channel: [ database ]

    main:
    ch_versions = Channel.empty()

    /*
    ================================================================================
                        QUAST: Genome assembly evaluation tool
    ================================================================================
    */

    // PROCESS: Run QUAST on the polished assembly for quality assessment and
    //  report the number of cleaned basepairs used to form the assembly
    QA_ASSEMBLY_QUAST (
        ch_cleaned_fastq_files.join(ch_assembly_file)
    )
    ch_versions = ch_versions.mix(QA_ASSEMBLY_QUAST.out.versions)

    // Collect assembly summaries and concatenate into one file
    ch_assembly_summary = Channel.empty()
    ch_assembly_summary = ch_assembly_summary
                            .mix(QA_ASSEMBLY_QUAST.out.summary_assemblies)
                            .collectFile(
                                name:       "Summary.Assemblies.tsv",
                                keepHeader: true,
                                storeDir:   "${params.outdir}/Summaries"
                            )

    ch_output_summary_files = ch_output_summary_files.mix(ch_assembly_summary)

    // Collect cleaned read/base summaries and concatenate into one file
    ch_cleaned_summary = Channel.empty()
    ch_cleaned_summary = ch_cleaned_summary
                            .mix(QA_ASSEMBLY_QUAST.out.summary_reads)
                            .collectFile(
                                name:     "Summary.CleanedReads-Bases.tsv",
                                keepHeader: true,
                                storeDir: "${params.outdir}/Summaries"
                            )

    ch_output_summary_files = ch_output_summary_files.mix(ch_cleaned_summary)

    // PROCESS: Calculate genome assembly depth of coverage
    CALCULATE_COVERAGE_UNIX (
        QA_ASSEMBLY_QUAST.out.qa_summaries
                          .join(ch_read_alignment_stats)
    )
    ch_versions = ch_versions.mix(CALCULATE_COVERAGE_UNIX.out.versions)

    // Collect genome coverage summaries and concatenate into one file
    ch_genome_cov_summary = Channel.empty()
    ch_genome_cov_summary = ch_genome_cov_summary
                                .mix(CALCULATE_COVERAGE_UNIX.out.summary)
                                .collectFile(
                                    name:     "Summary.GenomeCoverage.tsv",
                                    keepHeader: true,
                                    storeDir: "${params.outdir}/Summaries"
                                )

    ch_output_summary_files = ch_output_summary_files.mix(ch_genome_cov_summary)

    /*
    ================================================================================
                        GTDB-Tk: taxonomic classification using a GTDB reference
    ================================================================================
    */

    // PROCESS: Classify assembly FastA file using GTDB-Tk
    if (!params.skip_gtdbtk && params.gtdb_db) {
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
            ch_db_for_gtdbtk = GTDBTK_DB_PREPARATION_UNIX.out.db

        } else if ( ch_gtdbtk_db_file.isDirectory() ) {
            ch_db_for_gtdbtk = Channel
                                .fromPath( "${ch_gtdbtk_db_file}/*", type: 'dir', maxDepth: 1 )
                                .collect()
                                .map{
                                    [ it[0].getSimpleName(), it ]
                                }

        } else {
            error("Unsupported object given to --gtdb_db, database must be supplied as either a directory or a .tar.gz file!")
        }

        // PROCESS: Perform GTDBTk on assembly FastA file
        QA_ASSEMBLY_GTDBTK (
            ch_assembly_file,
            ch_db_for_gtdbtk,
            ch_mash_db_file
        )
        ch_versions             = ch_versions.mix(QA_ASSEMBLY_GTDBTK.out.versions)
        ch_output_summary_files = ch_output_summary_files.mix(QA_ASSEMBLY_GTDBTK.out.summary.map{ meta, file -> file })
    }

    /*
    ================================================================================
                        CAT: Contig Annotation Tool
    ================================================================================
    */

    // TODO: DB prep

    CLASSIFY_CONTIGS_CAT (
        ch_assembly_file,
        ch_cat_alignment_db_file,
        ch_cat_taxonomy_db_file
    )

    /*
    ================================================================================
                        CheckM2: Check for completeness and contamination
    ================================================================================
    */

    // TODO: DB prep

    CLASSIFY_ASSEMBLY_CHECKM2 (
        ch_assembly_file,
        ch_checkm2_db_file
    )

    /*
    ================================================================================
                        BUSCO: predict genes on contigs
    ================================================================================
    */

    // PROCESS: Classify contigs with BUSCO
    if (!params.skip_busco && params.busco_db) {
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
            BUSCO_DB_PREPARATION_UNIX(
                ch_busco_db
            )
            ch_versions     = ch_versions.mix(BUSCO_DB_PREPARATION_UNIX.out.versions)
            ch_db_for_busco = BUSCO_DB_PREPARATION_UNIX.out.db

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
        } else {
            error("Unsupported object given to --busco_db, database must be supplied as either a directory or a .tar.gz file!")
        }

        ch_lineage_for_busco_db = Channel
                                    .of(ch_busco_db_file)
                                    .map{
                                        db ->
                                            db = db.getSimpleName()
                                            db.contains('odb10') ? db : 'auto'
                                    }

        // PROCESS: Split assembly FastA file into individual contig files
        SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON (
            ch_assembly_file
        )
        ch_versions = ch_versions.mix(SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON.out.versions)

        // PROCESS: Perform BUSCO analysis on contigs
        QA_ASSEMBLY_BUSCO (
            SPLIT_MULTIFASTA_ASSEMBLY_BIOPYTHON.out.split_multifasta_assembly_dir,
            ch_lineage_for_busco_db,
            ch_db_for_busco,
            ch_busco_config_file
        )
        ch_versions             = ch_versions.mix(QA_ASSEMBLY_BUSCO.out.versions)
        ch_output_summary_files = ch_output_summary_files.mix(QA_ASSEMBLY_BUSCO.out.batch_summary.map{ meta, file -> file })
    }

    emit:
    versions      = ch_versions
}
