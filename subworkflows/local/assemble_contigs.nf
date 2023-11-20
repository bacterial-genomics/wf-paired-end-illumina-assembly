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
def qcfilecheck(qcfile, inputfile) {
    qcfile.map{ meta, file -> [ meta, [file] ] }
            .join(inputfile)
            .map{ meta, qc, input ->
                data = []
                qc.flatten().each{ data += it.readLines() }

                if ( data.any{ it.contains('FAIL') } ) {
                    log.warn("QC check failed for sample: ${data.last().split('\t').first()}")
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
    ch_versions      = Channel.empty()
    ch_qc_filechecks = Channel.empty()

    // Update meta to include meta.assembler
    if ( var_assembler_name == "SKESA" ) {
        // SKESA assembler
        // PROCESS: Run SKESA to assemble contigs with cleaned paired reads and cleaned singletons
        ASSEMBLE_CONTIGS_SKESA (
            ch_cleaned_reads
        )
        ch_versions = ch_versions.mix(ASSEMBLE_CONTIGS_SKESA.out.versions)
        ch_contigs = qcfilecheck(ASSEMBLE_CONTIGS_SKESA.out.qc_filecheck, ASSEMBLE_CONTIGS_SKESA.out.contigs)

        // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
        FILTER_CONTIGS_BIOPYTHON (
            ch_contigs
        )
        ch_versions = ch_versions.mix(FILTER_CONTIGS_BIOPYTHON.out.versions)

        // PROCESS: Use BWA/Samtools/Pilon to correct contigs with cleaned PE reads
        MAP_CONTIGS_BWA (
            ch_cleaned_reads.join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
        )
        ch_versions = ch_versions.mix(MAP_CONTIGS_BWA.out.versions)

        // Collect output files
        ch_bam_files      = qcfilecheck(MAP_CONTIGS_BWA.out.qc_filecheck, MAP_CONTIGS_BWA.out.bam)
        ch_assembly_file  = qcfilecheck(MAP_CONTIGS_BWA.out.qc_filecheck, MAP_CONTIGS_BWA.out.assembly)

        // Collect QC File Checks
        ch_qc_filechecks = ch_qc_filechecks
                                .mix(ASSEMBLE_CONTIGS_SKESA.out.qc_filecheck)
                                .mix(MAP_CONTIGS_BWA.out.qc_filecheck)
    } else {
        // Defaulting to SPAdes assembler
        // PROCESS: Run SPAdes to assemble contigs with cleaned paired reads and cleaned singletons
        ASSEMBLE_CONTIGS_SPADES (
            ch_cleaned_reads
        )
        ch_versions = ch_versions.mix(ASSEMBLE_CONTIGS_SPADES.out.versions)

        // ch_contigs = ASSEMBLE_CONTIGS_SPADES.out.contigs.map{ meta, file -> [ meta, [file] ] }
        ch_contigs = qcfilecheck(ASSEMBLE_CONTIGS_SPADES.out.qc_filecheck, ASSEMBLE_CONTIGS_SPADES.out.contigs)

        // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
        FILTER_CONTIGS_BIOPYTHON (
            ch_contigs
        )
        ch_versions = ch_versions.mix(FILTER_CONTIGS_BIOPYTHON.out.versions)

        // PROCESS: Use BWA/Samtools/Pilon to correct contigs with cleaned PE reads
        POLISH_ASSEMBLY_BWA_PILON (
            ch_cleaned_reads.join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
        )
        ch_versions = ch_versions.mix(POLISH_ASSEMBLY_BWA_PILON.out.versions)

        // Collect output files
        ch_bam_files = qcfilecheck(POLISH_ASSEMBLY_BWA_PILON.out.qc_filecheck, POLISH_ASSEMBLY_BWA_PILON.out.bam)
        ch_assembly_file  = qcfilecheck(POLISH_ASSEMBLY_BWA_PILON.out.qc_filecheck, POLISH_ASSEMBLY_BWA_PILON.out.assembly)

        // ch_bam_files = POLISH_ASSEMBLY_BWA_PILON.out.bam
        // ch_assembly_file = POLISH_ASSEMBLY_BWA_PILON.out.assembly
        // Collect QC File Checks
        ch_qc_filechecks = ch_qc_filechecks
                                .mix(ASSEMBLE_CONTIGS_SPADES.out.qc_filecheck)
                                .mix(POLISH_ASSEMBLY_BWA_PILON.out.qc_filecheck)
    }

    emit:
    bam_files     = ch_bam_files            // channel: [ val(meta), [{paired,single}.bam] ]
    assembly_file = ch_assembly_file        // channel: [ val(meta), [assembly.fna] ]
    qc_filecheck  = ch_qc_filechecks
    versions      = ch_versions
}
