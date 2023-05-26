/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSkesa.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.kraken1_db, params.kraken2_db, params.blast_db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet or directory not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// CONFIGS: Import configs for this workflow
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { INFILE_HANDLING_UNIX } from "../modules/local/infile_handling_unix/main"
include { REMOVE_PHIX_BBDUK } from "../modules/local/remove_phix_bbduk/main"
include { TRIM_READS_TRIMMOMATIC } from "../modules/local/trim_reads_trimmomatic/main"
//include { TRIM_READS_FASTP } from "../modules/local/trim_reads_fastp/main"
include { OVERLAP_PAIRED_READS_FLASH } from "../modules/local/overlap_paired_reads_flash/main"
//include { OVERLAP_PAIRED_READS_PEAR } from "../modules/local/overlap_paired_reads_pear/main"
include { READ_CLASSIFY_KRAKEN_ONE; READ_CLASSIFY_KRAKEN_TWO; } from "../modules/local/read_classify_kraken/main"
//include { READ_CLASSIFY_CENTRIFUGE } from "../modules/local/read_classify_centrifuge/main"
include { ASSEMBLE_SKESA } from "../modules/local/assemble_skesa/main"
include { FILTER_CONTIGS_BIOPYTHON } from "../modules/local/filter_contigs_biopython/main"
include { MAP_CONTIGS_BWA } from "../modules/local/map_contigs_bwa/main"
//include { POLISH_ASSEMBLY_UNICYCLER } from "../modules/local/polish_assembly_unicycler/main"
include { EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS } from "../modules/local/extract_read_alignment_depths_bedtools/main"
include { MLST_MLST } from "../modules/local/mlst_mlst/main"
//include { MLST_SRST2 } from "../modules/local/mlst_srst2/main"
include { ANNOTATE_PROKKA } from "../modules/local/annotate_prokka/main"
//include { ANNOTATE_BACTA } from "../modules/local/annotate_bacta/main"
include { EXTRACT_16S_BIOPYTHON } from "../modules/local/extract_16S_biopython/main"
include { EXTRACT_16S_BARRNAP } from "../modules/local/extract_16S_barrnap/main"
//include { 16S_EXTRACT_RNAMMER } from "../modules/local/16S_extract_rnammer/main"
include { ALIGN_16S_BLAST } from "../modules/local/align_16S_blast/main"
include { BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON } from "../modules/local/best_16S_blastn_bitscore_taxon_python/main"
include { QA_ASSEMBLY_QUAST } from "../modules/local/qa_assembly_quast/main"
//include { QA_ASSEMBLY_BUSCO } from "../modules/local/qa_assembly_busco/main"
//include { QA_ASSEMBLY_CAT } from "../modules/local/qa_assembly_cat/main"
//include { QA_ASSEMBLY_CHECKM2 } from "../modules/local/qa_assembly_checkm2/main"
include { CALCULATE_COVERAGE_UNIX } from "../modules/local/calculate_coverage_unix/main"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SKESA {

    // SETUP: Define empty channels to concatenate certain outputs
    ch_versions = Channel.empty()
    ch_alnstats_summary = Channel.empty()
    ch_blast_summary = Channel.empty()
    ch_ssu_species = Channel.empty()
    ch_genome_cov_summary = Channel.empty()
    ch_mlst_summary = Channel.empty()
    ch_assembly_summary = Channel.empty()
    ch_cleaned_summary = Channel.empty()
    ch_qc_filecheck = Channel.empty()

    // Check input for samplesheet or pull inputs from directory
    INPUT_CHECK (
        ch_input
    )

    // Check input files meet size criteria
    INFILE_HANDLING_UNIX (
        INPUT_CHECK.out.raw_reads
    )

    // Collect version info
    ch_versions = ch_versions.mix(INFILE_HANDLING_UNIX.out.versions)

    // PROCESS: Run bbduk to remove PhiX reads
    REMOVE_PHIX_BBDUK (
        INFILE_HANDLING_UNIX.out.input
    )

    // Collect version info
    ch_versions = ch_versions.mix(REMOVE_PHIX_BBDUK.out.versions)

    // PROCESS: Run trimmomatic to clip adapters and do quality trimming
    TRIM_READS_TRIMMOMATIC (
        REMOVE_PHIX_BBDUK.out.phix_removed
    )

    // Collect version info
    ch_versions = ch_versions.mix(TRIM_READS_TRIMMOMATIC.out.versions)

    // PROCESS: Run flash to merge overlapping sister reads into singleton reads
    OVERLAP_PAIRED_READS_FLASH (
        INFILE_HANDLING_UNIX.out.input.join(TRIM_READS_TRIMMOMATIC.out.trimmo)
    )

    // Collect version info
    ch_versions = ch_versions.mix(OVERLAP_PAIRED_READS_FLASH.out.versions)

    // PROCESS: Run kraken1 on paired cleaned reads
    READ_CLASSIFY_KRAKEN_ONE (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads
    )

    // Collect version info
    ch_versions = ch_versions.mix(READ_CLASSIFY_KRAKEN_ONE.out.versions)

    // PROCESS: Run kraken2 on paired cleaned reads
    READ_CLASSIFY_KRAKEN_TWO (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads
    )

    // Collect version info
    ch_versions = ch_versions.mix(READ_CLASSIFY_KRAKEN_TWO.out.versions)

    // PROCESS: Run SKESA to assemble contigs with cleaned paired reads and cleaned singletons
    ASSEMBLE_SKESA (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads
    )

    // Collect version info
    ch_versions = ch_versions.mix(ASSEMBLE_SKESA.out.versions)

    // PROCESS: Filter contigs based on length, coverage, GC skew, and compositional complexity
    FILTER_CONTIGS_BIOPYTHON (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads.join(ASSEMBLE_SKESA.out.contigs)
    )

    // Collect version info
    ch_versions = ch_versions.mix(FILTER_CONTIGS_BIOPYTHON.out.versions)

    // PROCESS: Use BWA/Samtools/Pilon to correct contigs with cleaned PE reads
    MAP_CONTIGS_BWA (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads.join(FILTER_CONTIGS_BIOPYTHON.out.uncorrected_contigs)
    )

    // Collect version info
    ch_versions = ch_versions.mix(MAP_CONTIGS_BWA.out.versions)

    // PROCESS: Run Bedtools to extract coverage from the pre-computed BAM alignment file
    EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS (
        MAP_CONTIGS_BWA.out.bam
    )

    // Collect all Summary Stats and concatenate into one file
    ch_alnstats_summary = ch_alnstats_summary.mix(EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.summary_alnstats)
    ch_alnstats_summary.collectFile(name: 'Summary.Illumina.CleanedReads-AlnStats.tab', storeDir: "${params.outdir}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.versions)

    // PROCESS: Run MLST to find MLST for each polished assembly
    MLST_MLST (
        MAP_CONTIGS_BWA.out.bam.join(MAP_CONTIGS_BWA.out.assembly)
    )

    // Collect all MLST Summaries and concatenate into one file
    ch_mlst_summary = ch_mlst_summary.mix(MLST_MLST.out.summary_mlst)
    ch_mlst_summary.collectFile(name: 'Summary.MLST.tab', storeDir: "${params.outdir}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(MLST_MLST.out.versions)

    // PROCESS: Annotate the polished assembly using Prokka
    ANNOTATE_PROKKA (
        MAP_CONTIGS_BWA.out.bam.join(MAP_CONTIGS_BWA.out.assembly)
    )

    // Collect version info
    ch_versions = ch_versions.mix(ANNOTATE_PROKKA.out.versions)

    // PROCESS: Attempt to extract 16S rRNA gene records from annotation file
    EXTRACT_16S_BIOPYTHON (
        ANNOTATE_PROKKA.out.annotation.join(MAP_CONTIGS_BWA.out.assembly)
    )

    // Collect version info
    ch_versions = ch_versions.mix(EXTRACT_16S_BIOPYTHON.out.versions)

    // PROCESS: Extract 16S rRNA gene sequences with Barrnap if missing from 16S_EXTRACT_BIOPYTHON
    EXTRACT_16S_BARRNAP (
        ANNOTATE_PROKKA.out.annotation.join(MAP_CONTIGS_BWA.out.assembly).join(EXTRACT_16S_BIOPYTHON.out.extracted_rna)
    )

    // Collect version info
    ch_versions = ch_versions.mix(EXTRACT_16S_BARRNAP.out.versions)

    // PROCESS: Run Blast on predicted 16S ribosomal RNA genes
    ALIGN_16S_BLAST (
        EXTRACT_16S_BARRNAP.out.extracted_base.join(MAP_CONTIGS_BWA.out.assembly)
    )

    // Collect version info
    ch_versions = ch_versions.mix(ALIGN_16S_BLAST.out.versions)

    // PROCESS: Filter Blast output for best alignment, based on bitscore
    BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON (
        ALIGN_16S_BLAST.out.blast_tsv.join(MAP_CONTIGS_BWA.out.assembly)
    )

    // Collect all BLAST Summaries and concatenate into one file
    ch_blast_summary = ch_blast_summary.mix(BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.blast_summary)
    ch_blast_summary.collectFile(name: 'Summary.16S.tab', storeDir: "${params.outdir}/qa")

    // Collect all BLAST Top Species Summaries and concatenate into one file
    ch_ssu_species = ch_ssu_species.mix(BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.ssu_species)
    ch_ssu_species.collectFile(name: '16S-top-species.tsv', storeDir: "${params.outdir}/ssu")

    // Collect version info
    ch_versions = ch_versions.mix(BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.versions)

    // PROCESS: Run QUAST on the polished assembly for quality assessment and
    //  report the number of cleaned basepairs used to form the assembly
    QA_ASSEMBLY_QUAST (
        OVERLAP_PAIRED_READS_FLASH.out.gzip_reads.join(MAP_CONTIGS_BWA.out.assembly)
    )

    // Collect all Assembly Summaries and concatenate into one file
    ch_assembly_summary = ch_assembly_summary.mix(QA_ASSEMBLY_QUAST.out.summary_assemblies)
    ch_assembly_summary.collectFile(name: 'Summary.Assemblies.tab', keepHeader: true, storeDir: "${params.outdir}/qa")

    // Collect all Cleaned Read/Base Summaries and concatenate into one file
    ch_cleaned_summary = ch_cleaned_summary.mix(QA_ASSEMBLY_QUAST.out.summary_reads)
    ch_cleaned_summary.collectFile(name: 'Summary.Illumina.CleanedReads-Bases.tab', storeDir: "${params.outdir}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(QA_ASSEMBLY_QUAST.out.versions)

    // PROCESS: Calculate genome assembly depth of coverage
    CALCULATE_COVERAGE_UNIX (
        QA_ASSEMBLY_QUAST.out.qa_summaries.join(EXTRACT_READ_ALIGNMENT_DEPTHS_BEDTOOLS.out.summary_stats)
    )

    // Collect all Genome Coverage Summaries and concatenate into one file
    ch_genome_cov_summary = ch_genome_cov_summary.mix(CALCULATE_COVERAGE_UNIX.out.genome_coverage)
    ch_genome_cov_summary.collectFile(name: 'Summary.Illumina.GenomeCoverage.tab', storeDir: "${params.outdir}/qa")

    // Collect version info
    ch_versions = ch_versions.mix(CALCULATE_COVERAGE_UNIX.out.versions)

    // PATTERN: Collate method for version information
    ch_versions.unique().collectFile(name: 'software_versions.yml', storeDir: params.logpath)

    // Collect all QC File Checks and concatenate into one file
    ch_qc_filecheck = ch_qc_filecheck.concat(
        INFILE_HANDLING_UNIX.out.qc_input_filecheck,
        REMOVE_PHIX_BBDUK.out.qc_phix_genome_filecheck,
        REMOVE_PHIX_BBDUK.out.qc_phix_removed_filecheck,
        TRIM_READS_TRIMMOMATIC.out.qc_adapters_filecheck,
        TRIM_READS_TRIMMOMATIC.out.qc_removed_adapters_filecheck,
        OVERLAP_PAIRED_READS_FLASH.out.qc_nonoverlap_filecheck,
        ASSEMBLE_SKESA.out.qc_raw_assembly_filecheck,
        MAP_CONTIGS_BWA.out.qc_filtered_asm_filecheck,
        MAP_CONTIGS_BWA.out.qc_pe_alignment_filecheck,
        MAP_CONTIGS_BWA.out.qc_corrected_asm_filecheck,
        MAP_CONTIGS_BWA.out.qc_se_alignment_filecheck,
        ANNOTATE_PROKKA.out.qc_annotated_filecheck,
        EXTRACT_16S_BARRNAP.out.qc_ssu_extracted_filecheck,
        EXTRACT_16S_BARRNAP.out.qc_ssu_renamed_filecheck,
        ALIGN_16S_BLAST.out.qc_blastn_filecheck,
        BEST_16S_BLASTN_BITSCORE_TAXON_PYTHON.out.qc_filtered_blastn_filecheck
    )
    ch_qc_filecheck.collectFile(name: 'Summary.QC_File_Checks.tab', storeDir: "${params.outdir}/qa", sort: {it.getSimpleName()})

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
