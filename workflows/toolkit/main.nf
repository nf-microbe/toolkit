/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'

include { CUSTOM_CLEANWORKDIRS      } from '../../modules/nf-core/custom/cleanworkdirs/main'
include { MULTIQC                   } from '../../modules/nf-core/multiqc/main'

include { FASTQ_READPREPROCESSING_FASTQ } from '../../subworkflows/nf-core/fastq_readpreprocessing_fastq'
include { getGenomeAttribute            } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'
include { methodsDescriptionText        } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'
include { paramsSummaryMultiqc          } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TOOLKIT {

    take:
    input_fastqs    // channel: fastq files read in from --input (samplesheet) and --fastqs
    input_fastas    // channel: fasta files read in from --input (samplesheet) and --fastas
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
    ch_preprocessed_fastq_gz    = FASTQ_READPREPROCESSING_FASTQ.out.preprocessed_fastq_gz
    ch_multiqc_files            = ch_multiqc_files.mix(FASTQ_READPREPROCESSING_FASTQ.out.multiqc_files)
    ch_versions                 = FASTQ_READPREPROCESSING_FASTQ.out.versions
    ch_workdirs_to_clean        = FASTQ_READPREPROCESSING_FASTQ.out.workdirs_to_clean

    //
    // MODULE: Clean intermediate files
    //
    if (parameters.remove_intermediate_files) {
        //
        // MODULE: Clean up intermediate working directories
        //
        ch_workdirs_to_clean_unique = ch_workdirs_to_clean.unique()
        CUSTOM_CLEANWORKDIRS(ch_workdirs_to_clean_unique)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${parameters.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = parameters.multiqc_config ?
        Channel.fromPath(parameters.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = parameters.multiqc_logo ?
        Channel.fromPath(parameters.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = parameters.multiqc_methods_description ?
        file(parameters.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
