/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQ_READPREPROCESSING_FASTQ } from '../../subworkflows/nf-core/fastq_readpreprocessing_fastq'
include { getGenomeAttribute            } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow READPREPROCESSING {

    take:
    input_fastqs    // channel: fastq files read in from --input (samplesheet) and/or --fastqs
    parameters      // pipeline parameters

    main:

    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()

    // Prepare host fasta and bowtie2 index for host removal
    if (parameters.run_bowtie2_host_removal) {
        // Load igenomes fasta and bowtie2 index
        igenomes_fasta   = getGenomeAttribute('fasta')
        igenomes_index   = getGenomeAttribute('bowtie2')

        // create bowtie2 fasta and index channels based on input parameters
        if (parameters.genome) {
            ch_bowtie2_fasta    = file(igenomes_fasta, checkIfExists: true)
            ch_bowtie2_index    = file(igenomes_index, checkIfExists: true)
        } else if (parameters.host_bowtie2_index) {
            if (!parameters.host_fasta) {
                error "[nf-microbe/toolkit]: --run_bowtie2_host_removal = true but no iGenomes (--genome) or custom fasta (--host_fasta) provided"
            }
            ch_bowtie2_fasta    = file(parameters.host_fasta, checkIfExists: true)
            ch_bowtie2_index    = file(parameters.host_bowtie2_index, checkIfExists: true)
        } else {
            if (!parameters.host_fasta) {
                error "[nf-microbe/toolkit]: --run_bowtie2_host_removal = true but no iGenomes (--genome) or custom fasta (--host_fasta) provided"
            }
            ch_bowtie2_fasta    = file(parameters.host_fasta, checkIfExists: true)
            ch_bowtie2_index    = null
        }
    } else {
        ch_bowtie2_fasta    = null
        ch_bowtie2_index    = null
    }

    //
    // SUBWORKFLOW: Read preprocessing
    //
    FASTQ_READPREPROCESSING_FASTQ(
        input_fastqs,                           // channel: [ [ meta.id, meta.single_end, meta.run ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
        parameters.run_fastp,                   // boolean: run fastp
        parameters.perform_run_merging,         // boolean: merge runs within samples
        parameters.run_bowtie2_host_removal,    // boolean: run host removal
        ch_bowtie2_fasta,                       // string: /path/to/fasta.fasta.gz (OPTIONAL)
        ch_bowtie2_index                        // string: /path/to/host_index.bt2 (OPTIONAL)
    )

    emit:
    preprocessed_fastq_gz   = FASTQ_READPREPROCESSING_FASTQ.out.preprocessed_fastq_gz
    multiqc_files           = ch_multiqc_files.mix(FASTQ_READPREPROCESSING_FASTQ.out.multiqc_files)
    versions                = FASTQ_READPREPROCESSING_FASTQ.out.versions
    workdirs_to_clean       = FASTQ_READPREPROCESSING_FASTQ.out.workdirs_to_clean
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
