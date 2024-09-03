/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap              } from 'plugin/nf-validation'

include { CUSTOM_CLEANWORKDIRS          } from '../../modules/nf-core/custom/cleanworkdirs/main'
include { MULTIQC                       } from '../../modules/nf-core/multiqc/main'

include { getGenomeAttribute            } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'
include { methodsDescriptionText        } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'

include { FASTQ_READASSEMBLY_FASTA      } from '../../subworkflows/nf-core/fastq_readassembly_fasta'
include { paramsSummaryMultiqc          } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

include { READPREPROCESSING             } from '../readpreprocessing'


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

    //
    // WORKFLOW: Read preprocessing
    //
    READPREPROCESSING(
        input_fastqs,   // channel: [ [ meta.id, meta.single_end, meta.run ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
        parameters      // pipeline parameters
    )
    ch_preprocessed_fastq_gz    = READPREPROCESSING.out.preprocessed_fastq_gz
    ch_multiqc_files            = ch_multiqc_files.mix(READPREPROCESSING.out.multiqc_files)
    ch_versions                 = READPREPROCESSING.out.versions
    ch_workdirs_to_clean        = READPREPROCESSING.out.workdirs_to_clean

    //
    // WORKFLOW: Read assembly
    //
    FASTQ_READASSEMBLY_FASTA(
        ch_preprocessed_fastq_gz,           // channel: [ [ meta.id, meta.single_end, meta.run ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
        parameters.run_megahit_single,      // boolean: false
        parameters.run_megahit_coassembly,  // boolean: false
        parameters.run_spades_single,       // boolean: false
        parameters.run_spades_coassembly,   // boolean: false
        parameters.use_spades_scaffolds,    // boolean: false
        parameters.run_penguin_single,      // boolean: false
        parameters.run_penguin_coassembly,  // boolean: false
    )
    ch_assemblies_fasta_gz  = FASTQ_READASSEMBLY_FASTA.out.assemblies_fasta_gz
    ch_assembly_graph_gz    = FASTQ_READASSEMBLY_FASTA.out.assembly_graph_gz
    ch_spades_logs          = FASTQ_READASSEMBLY_FASTA.out.spades_logs
    ch_multiqc_files        = ch_multiqc_files.mix(FASTQ_READASSEMBLY_FASTA.out.multiqc_files)
    ch_versions             = FASTQ_READASSEMBLY_FASTA.out.versions

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
