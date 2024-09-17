//
// Assemble contigs using tool specified via `--assembler {spades,skesa}`
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { ASSEMBLE_CONTIGS_SKESA    } from "../../modules/local/assemble_contigs_skesa/main"
include { ASSEMBLE_CONTIGS_SPADES   } from "../../modules/local/assemble_contigs_spades/main"

include { FILTER_CONTIGS_BIOPYTHON  } from "../../modules/local/filter_contigs_biopython/main"

include { MAP_CONTIGS_BWA           } from "../../modules/local/map_contigs_bwa/main"
include { POLISH_ASSEMBLY_BWA_PILON } from "../../modules/local/polish_assembly_bwa_pilon/main"
// include { POLISH_ASSEMBLY_UNICYCLER } from "../../modules/local/polish_assembly_unicycler/main"

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

workflow ASSEMBLE_CONTIGS {

    take:
    ch_cleaned_reads    // channel: [ val(meta), [cleaned_fastq_files (R1, R2, single)] ]
    var_assembler_name  // var (str): assembler_name

    main:
    ch_versions       = Channel.empty()
    ch_checksums_file = Channel.empty()
    ch_qc_filechecks  = Channel.empty()

    // Update meta to include meta.assembler
    if ( var_assembler_name == "SKESA" ) {
        // SKESA assembler
        // PROCESS: Run SKESA to assemble contigs with cleaned paired reads and cleaned singletons
        //          which skips post-assembly mapping for SNP and InDel corrections too for speed.
        ASSEMBLE_CONTIGS_SKESA (
            ch_cleaned_reads
        )
        ch_versions = ch_versions.mix(ASSEMBLE_CONTIGS_SKESA.out.versions)
        ch_contigs  = qcfilecheck(
                        "ASSEMBLE_CONTIGS_SKESA",
                        ASSEMBLE_CONTIGS_SKESA.out.qc_filecheck,
                        ASSEMBLE_CONTIGS_SKESA.out.contigs
                      )

        // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
        FILTER_CONTIGS_BIOPYTHON (
            ch_contigs
        )
        ch_versions = ch_versions.mix(FILTER_CONTIGS_BIOPYTHON.out.versions)

        // PROCESS: Create BAM file for depth of coverage calculations
        MAP_CONTIGS_BWA (
            ch_cleaned_reads.join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
        )
        ch_versions = ch_versions.mix(MAP_CONTIGS_BWA.out.versions)

        // Collect output files
        ch_bam_files      = qcfilecheck(
                                "MAP_CONTIGS_BWA",
                                MAP_CONTIGS_BWA.out.qc_filecheck,
                                MAP_CONTIGS_BWA.out.bam
                            )
        ch_assembly_file  = qcfilecheck(
                                "MAP_CONTIGS_BWA",
                                MAP_CONTIGS_BWA.out.qc_filecheck,
                                MAP_CONTIGS_BWA.out.assembly
                            )

        // Collect QC File Checks
        ch_checksums_file = ch_checksums_file.mix(MAP_CONTIGS_BWA.out.checksums)
        ch_qc_filechecks  = ch_qc_filechecks
                                .mix(ASSEMBLE_CONTIGS_SKESA.out.qc_filecheck)
                                .mix(MAP_CONTIGS_BWA.out.qc_filecheck)

    } else {
        // Defaulting to SPAdes assembler
        // PROCESS: Run SPAdes to assemble contigs with cleaned paired reads and cleaned singletons
        ASSEMBLE_CONTIGS_SPADES (
            ch_cleaned_reads
        )
        ch_versions = ch_versions.mix(ASSEMBLE_CONTIGS_SPADES.out.versions)
        ch_contigs = qcfilecheck(
                        "ASSEMBLE_CONTIGS_SPADES",
                        ASSEMBLE_CONTIGS_SPADES.out.qc_filecheck,
                        ASSEMBLE_CONTIGS_SPADES.out.contigs
                    )

        // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
        FILTER_CONTIGS_BIOPYTHON (
            ch_contigs
        )
        ch_versions = ch_versions.mix(FILTER_CONTIGS_BIOPYTHON.out.versions)

        // PROCESS: Use BWA/Samtools/Pilon to SNP and InDel correct contigs with cleaned PE reads
        // NOTE:  The "path(cleaned_fastq_files)" is already input to this POLISH channel, but
        //        currently just the meta.id is used for the readset. Should be using the
        //        path(cleaned_fastq_files) items though.
        POLISH_ASSEMBLY_BWA_PILON (
            ch_cleaned_reads.join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
        )
        ch_versions = ch_versions.mix(POLISH_ASSEMBLY_BWA_PILON.out.versions)

        // Collect output files
        ch_bam_files      = qcfilecheck(
                                "POLISH_ASSEMBLY_BWA_PILON",
                                POLISH_ASSEMBLY_BWA_PILON.out.qc_filecheck,
                                POLISH_ASSEMBLY_BWA_PILON.out.bam
                            )
        ch_assembly_file  = qcfilecheck(
                                "POLISH_ASSEMBLY_BWA_PILON",
                                POLISH_ASSEMBLY_BWA_PILON.out.qc_filecheck,
                                POLISH_ASSEMBLY_BWA_PILON.out.assembly
                            )

        // Collect QC File Checks
        ch_checksums_file = ch_checksums_file.mix(POLISH_ASSEMBLY_BWA_PILON.out.checksums)
        ch_qc_filechecks  = ch_qc_filechecks
                                .mix(ASSEMBLE_CONTIGS_SPADES.out.qc_filecheck)
                                .mix(POLISH_ASSEMBLY_BWA_PILON.out.qc_filecheck)
    }

    emit:
    bam_files     = ch_bam_files            // channel: [ val(meta), [{paired,single}.bam] ]
    assembly_file = ch_assembly_file        // channel: [ val(meta), [assembly.fna] ]
    qc_filecheck  = ch_qc_filechecks
    checksums     = ch_checksums_file       // channel: [ val(meta), [assembly.fna] ]
    versions      = ch_versions
}
