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

assembler = toLower(params.assembler)

// Check QC filechecks for a failure
def checkQCFilechecks(it) {
    it.flatten().map{
        file ->
            // Obtain file contents
            getData = file.getText()

            // Add header if needed
            if ( !getData.split('\n').first().contains('Sample name') ) {
                file.write("Sample name\tQC step\tOutcome (Pass/Fail)\n")
                file.append(getData)
            }

            // Move to QC log directory
            file.copyTo(params.qc_filecheck_log_dir)

            // Check file contents for failure
            if ( getData.contains('FAIL') ) {
                error("${file.getBaseName().split('\\.').last().replace('_', ' ')} check failed!")
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
    ch_cleaned_reads    // channel: [ val(meta), [read1], [read2], [single] ]

    main:
    ch_versions = Channel.empty()
    ch_qc_filechecks = Channel.empty()

    if ( assembler == "skesa" ) {
        // PROCESS: Run SKESA to assemble contigs with cleaned paired reads and cleaned singletons
        ASSEMBLE_CONTIGS_SKESA (
            ch_cleaned_reads
        )
        ch_versions = ch_versions.mix(ASSEMBLE_CONTIGS_SKESA.out.versions)
        checkQCFilechecks(ASSEMBLE_CONTIGS_SKESA.out.qc_filecheck)

        // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
        FILTER_CONTIGS_BIOPYTHON (
            ch_cleaned_reads.join(ASSEMBLE_CONTIGS_SKESA.out.contigs)
        )
        ch_versions = ch_versions.mix(FILTER_CONTIGS_BIOPYTHON.out.versions)

        // PROCESS: Use BWA/Samtools/Pilon to correct contigs with cleaned PE reads
        MAP_CONTIGS_BWA (
            ch_cleaned_reads.join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
        )
        ch_versions = ch_versions.mix(MAP_CONTIGS_BWA.out.versions)
        checkQCFilechecks(MAP_CONTIGS_BWA.out.qc_filecheck)

        // Collect output files
        bam      = MAP_CONTIGS_BWA.out.bam
        assembly = MAP_CONTIGS_BWA.out.assembly

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
        checkQCFilechecks(ASSEMBLE_CONTIGS_SPADES.out.qc_filecheck)

        // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
        FILTER_CONTIGS_BIOPYTHON (
            ch_cleaned_reads.join(ASSEMBLE_CONTIGS_SPADES.out.contigs)
        )
        ch_versions = ch_versions.mix(FILTER_CONTIGS_BIOPYTHON.out.versions)

        // PROCESS: Use BWA/Samtools/Pilon to correct contigs with cleaned PE reads
        POLISH_ASSEMBLY_BWA_PILON (
            ch_cleaned_reads.join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
        )
        ch_versions = ch_versions.mix(POLISH_ASSEMBLY_BWA_PILON.out.versions)
        checkQCFilechecks(POLISH_ASSEMBLY_BWA_PILON.out.qc_filecheck)

        // Collect output files
        bam      = POLISH_ASSEMBLY_BWA_PILON.out.bam
        assembly = POLISH_ASSEMBLY_BWA_PILON.out.assembly

        // Collect QC File Checks
        ch_qc_filechecks = ch_qc_filechecks
                                .mix(ASSEMBLE_CONTIGS_SPADES.out.qc_filecheck)
                                .mix(POLISH_ASSEMBLY_BWA_PILON.out.qc_filecheck)
    }

    emit:
    bam           = bam             // channel: [ val(meta), [paired.bam], [single.bam] ]
    assembly      = assembly        // channel: [ val(meta), [assembly.fna] ]
    versions      = ch_versions
    qc_filecheck  = ch_qc_filechecks
}
