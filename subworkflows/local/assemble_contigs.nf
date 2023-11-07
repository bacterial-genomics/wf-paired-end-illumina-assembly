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

// Convert params.assembler to lowercase
def toLower(it) {
    it.toString().toLowerCase()
}

params.assembler = toLower(params.assembler)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ASSEMBLE_CONTIGS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLE_CONTIGS {

    take:
    ch_cleaned_reads

    main:
    ch_versions = Channel.empty()
    ch_qc_filechecks = Channel.empty()

    if ( params.assembler == "skesa" ) {
        // PROCESS: Run SKESA to assemble contigs with cleaned paired reads and cleaned singletons
        ASSEMBLE_CONTIGS_SKESA (
            ch_cleaned_reads
        )

        // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
        FILTER_CONTIGS_BIOPYTHON (
            ch_cleaned_reads
                .join(ASSEMBLE_CONTIGS_SKESA.out.contigs)
        )

        // PROCESS: Use BWA/Samtools/Pilon to correct contigs with cleaned PE reads
        MAP_CONTIGS_BWA (
            ch_cleaned_reads
                .join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
        )

        // Collect output files
        bam      = MAP_CONTIGS_BWA.out.bam
        assembly = MAP_CONTIGS_BWA.out.assembly

        // Collect version info
        ch_versions = ch_versions
                        .mix(ASSEMBLE_CONTIGS_SKESA.out.versions)
                        .mix(FILTER_CONTIGS_BIOPYTHON.out.versions)
                        .mix(MAP_CONTIGS_BWA.out.versions)

        // Collect QC File Checks
        ch_qc_filechecks = ch_qc_filechecks
                            .concat(
                                ASSEMBLE_CONTIGS_SKESA.out.qc_raw_assembly_filecheck,
                                MAP_CONTIGS_BWA.out.qc_filtered_asm_filecheck,
                                MAP_CONTIGS_BWA.out.qc_pe_alignment_filecheck,
                                MAP_CONTIGS_BWA.out.qc_corrected_asm_filecheck,
                                MAP_CONTIGS_BWA.out.qc_se_alignment_filecheck
                            )
    } else {
        // Defaulting to SPAdes assembler
        // PROCESS: Run SPAdes to assemble contigs with cleaned paired reads and cleaned singletons
        ASSEMBLE_CONTIGS_SPADES (
            ch_cleaned_reads
        )

        // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
        FILTER_CONTIGS_BIOPYTHON (
            ch_cleaned_reads
                .join(ASSEMBLE_CONTIGS_SPADES.out.contigs)
        )

        // PROCESS: Use BWA/Samtools/Pilon to correct contigs with cleaned PE reads
        POLISH_ASSEMBLY_BWA_PILON (
            ch_cleaned_reads
                .join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
        )

        // Collect output files
        bam      = POLISH_ASSEMBLY_BWA_PILON.out.bam
        assembly = POLISH_ASSEMBLY_BWA_PILON.out.assembly

        // Collect version info
        ch_versions = ch_versions
                        .mix(ASSEMBLE_CONTIGS_SPADES.out.versions)
                        .mix(FILTER_CONTIGS_BIOPYTHON.out.versions)
                        .mix(POLISH_ASSEMBLY_BWA_PILON.out.versions)

        // Collect QC File Checks
        ch_qc_filechecks = ch_qc_filechecks
                                .concat(
                                    ASSEMBLE_CONTIGS_SPADES.out.qc_raw_assembly_filecheck,
                                    POLISH_ASSEMBLY_BWA_PILON.out.qc_filtered_asm_filecheck,
                                    POLISH_ASSEMBLY_BWA_PILON.out.qc_pe_alignment_filecheck,
                                    POLISH_ASSEMBLY_BWA_PILON.out.qc_polished_asm_filecheck,
                                    POLISH_ASSEMBLY_BWA_PILON.out.qc_corrected_asm_filecheck,
                                    POLISH_ASSEMBLY_BWA_PILON.out.qc_se_alignment_filecheck
                                )
    }

    emit:
    bam           = bam
    assembly      = assembly
    versions      = ch_versions
    qc_filechecks = ch_qc_filechecks
}
