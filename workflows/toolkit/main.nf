/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Import plugins
include { paramsSummaryMap              } from 'plugin/nf-validation'

// Import modules
include { CAT_FASTQ as RUNMERGE_FASTQ   } from '../../modules/nf-core/cat/fastq'
include { CUSTOM_CLEANWORKDIRS          } from '../../modules/nf-core/custom/cleanworkdirs'
include { FASTQC as FASTQC_RAW          } from '../../modules/nf-core/fastqc'
include { FASTQC as FASTQC_PREPROCESSED } from '../../modules/nf-core/fastqc'
include { FASTP                         } from '../../modules/nf-core/fastp'
include { LOGAN_CONTIGAWSCLI            } from '../../modules/nf-core/logan/contigawscli'
include { SEQKIT_SEQ as COBRA_SEQ       } from '../../modules/nf-core/seqkit/seq'
include { SEQKIT_FX2TAB as COBRA_FX2TAB } from '../../modules/nf-core/seqkit/fx2tab'
include { SRA_SRATOOLS                  } from '../../modules/nf-core/sra/sratools'
include { MULTIQC                       } from '../../modules/nf-core/multiqc'

// Import functions from subworkflows
include { getGenomeAttribute            } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'
include { getWorkDirs; rmEmptyFastQs    } from '../../subworkflows/nf-core/utils_nfmicrobe_functions'
include { methodsDescriptionText        } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'

// Import subworkflows
include { ACCESSION_LOGAN_FASTA                 } from '../../subworkflows/nf-core/accession_logan_fasta'
include { FASTA_ASSEMBLYQC_FASTA                } from '../../subworkflows/nf-core/fasta_assemblyqc_fasta'
include { FASTA_CHECKV_FASTATSV                 } from '../../subworkflows/nf-core/fasta_checkv_fastatsv'
include { FASTA_MGECLASSIFICATION_FASTATSV      } from '../../subworkflows/nf-core/fasta_mgeclassification_fastatsv'
include { FASTQ_BOWTIE2_FASTQ                   } from '../../subworkflows/nf-core/fastq_bowtie2_fastq'
include { FASTQ_READASSEMBLY_FASTA              } from '../../subworkflows/nf-core/fastq_readassembly_fasta'
include { FASTQFASTAGFA_ASSEMBLYEXTENSION_FASTA } from '../../subworkflows/nf-core/fastqfastagfa_assemblyextension_fasta'
include { paramsSummaryMultiqc                  } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

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

    if (parameters.sra_accessions) {
        // Load SRA accessions one-by-one from file
        ch_sra_accessions = Channel.fromPath(parameters.sra_accessions)
            .splitCsv(header:false)
            .map { row -> [ [ id: row[0], group: row[0]], row[0] ] }
    }

    /*
    -------------------------------------------------
        READ DOWNLOAD
    -------------------------------------------------
    */
    if (ch_sra_accessions && parameters.download_sra_fastqs) {
        //
        // MODULE: Download SRA FastQ files
        //
        SRA_SRATOOLS(
            ch_sra_accessions.map { meta, acc -> [ [ id: "sra_" + meta.id, group: "sra_" + meta.id ], acc ] }
        )
        ch_versions     = ch_versions.mix(SRA_SRATOOLS.out.versions)
        input_fastqs    = input_fastqs
            .mix(
                SRA_SRATOOLS.out.fastq.map { meta, fastq ->
                    meta.single_end = fastq.size() == 2 ? false : true
                    [ meta, fastq ]
                }
            )
    }

    /*
    -------------------------------------------------
        RAW READ QC
    -------------------------------------------------
    */
    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW (
        input_fastqs
    )
    ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{ zip -> zip[1] })
    ch_versions         = ch_versions.mix(FASTQC_RAW.out.versions.first())

    /*
    -------------------------------------------------
        READ PREPROCESSING
    -------------------------------------------------
    */
    if (parameters.run_fastp) {
        //
        // MODULE: Run fastp on raw reads
        //
        FASTP(
            input_fastqs,
            [],
            false,
            false,
            false
        )
        ch_preprocessed_prefilt_fastq_gz    = FASTP.out.reads
        ch_multiqc_files                    = ch_multiqc_files.mix(FASTP.out.json.collect{ json -> json[1] })
        ch_versions                         = ch_versions.mix(FASTP.out.versions)

        // REMOVE EMPTY FASTQ FILES
        ch_preprocessed_fastq_gz = rmEmptyFastQs(ch_preprocessed_prefilt_fastq_gz, false)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_preprocessing_workdirs = getWorkDirs(
            input_fastqs,
            ch_preprocessed_fastq_gz,
            false
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
            ch_pre_preprocessing_workdirs.map { meta, dir -> [ meta, dir, 'PREPROCESSING' ] }
        )
    } else {
        ch_preprocessed_fastq_gz = input_fastqs
    }

    /*
    -------------------------------------------------
        RUN MERGING
    -------------------------------------------------
    */
    if (parameters.perform_run_merging) {
        // prepare reads for concatenating within runs
        ch_reads_forcat = ch_preprocessed_fastq_gz
            .map { meta, reads -> [ meta - meta.subMap('run'), reads ] }
            .groupTuple(sort: 'deep')
            .branch {
                meta, reads ->
                    cat:      reads.size() >= 2 // SE: [ [ meta ], [ S1_R1, S2_R1 ] ]; PE: [ [ meta ], [ [ S1_R1, S1_R2 ], [ S2_R1, S2_R2 ] ] ]
                    skip_cat: true              // Can skip merging if only single lanes
            }

        //
        // MODULE: Concatenate reads across runs, within samples
        //
        RUNMERGE_FASTQ(
            ch_reads_forcat.cat.map { meta, reads -> [ meta, reads.flatten() ] }
        )
        ch_runmerged_fastq_gz   = RUNMERGE_FASTQ.out.reads
        ch_versions             = ch_versions.mix(RUNMERGE_FASTQ.out.versions)

        // Ensure we don't have nests of nests so that structure is in form expected for assembly
        ch_reads_forcat_skipped = ch_reads_forcat.skip_cat
            .map { meta, reads ->
                def new_reads = meta.single_end ? reads[0] : reads.flatten()
                [ meta, new_reads ]
            }

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_runmerge_workdirs = getWorkDirs(
            ch_preprocessed_fastq_gz.map { meta, reads ->
                if (meta.single_end) {
                    [ meta - meta.subMap('run'), reads ]
                } else {
                    [ meta - meta.subMap('run'), reads[0] ]
                }
            },
            ch_runmerged_fastq_gz,
            false
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
            ch_pre_runmerge_workdirs.map { meta, dir -> [ meta, dir, 'RUNMERGING' ] }
        )

        // Combine single run and multi-run-merged data
        ch_runmerged_fastq_gz = ch_runmerged_fastq_gz.mix(ch_reads_forcat_skipped)
    } else {
        ch_runmerged_fastq_gz = ch_preprocessed_fastq_gz
    }

    /*
    -------------------------------------------------
        HOST READ REMOVAL
    -------------------------------------------------
    */
    // Load igenomes fasta and bowtie2 index
    igenomes_fasta   = getGenomeAttribute('fasta')
    igenomes_index   = getGenomeAttribute('bowtie2')

    if (parameters.run_bowtie2_host_removal) {
        // create host fasta and index channels based on input parameters
        if (parameters.genome) {
            ch_host_fasta_gz    = file(igenomes_fasta, checkIfExists: true)
            ch_bowtie2_index    = file(igenomes_index, checkIfExists: true)
        } else if (parameters.host_bowtie2_index) {
            if (!parameters.host_fasta) {
                error "[nf-microbe/toolkit]: --run_bowtie2_host_removal = true but no iGenomes (--genome) or custom fasta (--host_fasta) provided"
            }
            ch_host_fasta_gz    = file(parameters.host_fasta, checkIfExists: true)
            ch_bowtie2_index    = file(parameters.host_bowtie2_index, checkIfExists: true)
        } else {
            if (!parameters.host_fasta) {
                error "[nf-microbe/toolkit]: --run_bowtie2_host_removal = true but no iGenomes (--genome) or custom fasta (--host_fasta) provided"
            }
            ch_host_fasta_gz    = file(parameters.host_fasta, checkIfExists: true)
            ch_bowtie2_index    = null
        }
    }

    if (parameters.run_bowtie2_host_removal) {
        //
        // SUBWORKFLOW: Remove host reads using Bowtie2
        //
        FASTQ_BOWTIE2_FASTQ(
            ch_runmerged_fastq_gz,
            ch_host_fasta_gz,
            ch_bowtie2_index
        )
        ch_hostremoved_prefilt_fastq_gz = FASTQ_BOWTIE2_FASTQ.out.fastq_gz
        ch_multiqc_files                = ch_multiqc_files.mix(FASTQ_BOWTIE2_FASTQ.out.bt2_log.collect{ bt2_log -> bt2_log[1] })
        ch_versions                     = ch_versions.mix(FASTQ_BOWTIE2_FASTQ.out.versions)

        // REMOVE EMPTY FASTQ FILES FROM CHANNEL
        ch_hostremoved_fastq_gz = rmEmptyFastQs(ch_hostremoved_prefilt_fastq_gz, false)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_hostremove_workdirs = getWorkDirs(
            ch_runmerged_fastq_gz.map { meta, reads -> [ meta, reads[0] ] },
            ch_hostremoved_prefilt_fastq_gz,
            false
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
            ch_pre_hostremove_workdirs.map { meta, dir -> [ meta, dir, 'HOSTREMOVAL' ] }
        )
    } else {
        ch_hostremoved_fastq_gz = ch_runmerged_fastq_gz
    }

    if (parameters.run_fastp ||
        parameters.perform_run_merging ||
        parameters.run_bowtie2_host_removal) {
        //
        // MODULE: Run FastQC on preprocessed reads
        //
        FASTQC_PREPROCESSED (
            ch_hostremoved_fastq_gz
        )
        ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_PREPROCESSED.out.zip.collect{ zip -> zip[1]})
        ch_versions         = ch_versions.mix(FASTQC_PREPROCESSED.out.versions.first())
    }

    /*
    -------------------------------------------------
        LOGAN MULTIPLIER DOWNLOAD
    -------------------------------------------------
    */
    if (ch_sra_accessions && (parameters.download_logan_mult_contigs || parameters.download_logan_mult_unitigs) ) {
        //
        // SUBWORKFLOW: Download Logan multiplied assemblies
        //
        ACCESSION_LOGAN_FASTA(
            ch_sra_accessions,
            parameters.download_logan_mult_contigs,
            parameters.download_logan_mult_unitigs
        )
        ch_versions             = ch_versions.mix(ACCESSION_LOGAN_FASTA.out.versions)
        ch_hostremoved_fastq_gz = ch_hostremoved_fastq_gz
            .mix(
                ACCESSION_LOGAN_FASTA.out.logan_mult_fasta_gz.map { meta, fasta ->
                    [ meta + [ single_end: true ], fasta ]
                }
            )
    }

    /*
    -------------------------------------------------
        ASSEMBLY DOWNLOAD
    -------------------------------------------------
    */
    if (ch_sra_accessions && parameters.download_logan_contigs) {
        //
        // MODULE: Download Logan contigs
        //
        LOGAN_CONTIGAWSCLI(
            ch_sra_accessions.map { meta, acc ->
                [ [ id: "logan_assemblies_" + meta.id, group: "logan_assemblies_" + meta.id, assembler: 'minia', assembly_method: 'single' ], acc ] }
        )
        ch_versions     = ch_versions.mix(LOGAN_CONTIGAWSCLI.out.versions)
        input_fastas    = input_fastas.mix(LOGAN_CONTIGAWSCLI.out.fasta)
    }

    /*
    -------------------------------------------------
        READ ASSEMBLY
    -------------------------------------------------
    */
    if (parameters.run_megahit_single ||
        parameters.run_megahit_coassembly ||
        parameters.run_spades_single ||
        parameters.run_spades_coassembly ||
        parameters.run_penguin_single ||
        parameters.run_penguin_coassembly) {
            //
            // SUBWORKFLOW: Read assembly
            //
            FASTQ_READASSEMBLY_FASTA(
                ch_hostremoved_fastq_gz,
                parameters.run_megahit_single,
                parameters.run_megahit_coassembly,
                parameters.run_spades_single,
                parameters.run_spades_coassembly,
                parameters.use_spades_scaffolds,
                parameters.run_penguin_single,
                parameters.run_penguin_coassembly,
            )
            input_fastas            = input_fastas.mix(FASTQ_READASSEMBLY_FASTA.out.assemblies_fasta_gz)
            ch_assembly_graph_gz    = FASTQ_READASSEMBLY_FASTA.out.assembly_graph_gz
            ch_assembly_logs        = FASTQ_READASSEMBLY_FASTA.out.assembly_logs
            ch_hostremoved_fastq_gz = ch_hostremoved_fastq_gz.mix(FASTQ_READASSEMBLY_FASTA.out.coassembly_fastq_gz)
            ch_multiqc_files        = ch_multiqc_files.mix(FASTQ_READASSEMBLY_FASTA.out.multiqc_files)
            ch_versions             = FASTQ_READASSEMBLY_FASTA.out.versions
        } else {
            ch_assembly_graph_gz = []
            ch_assembly_logs     = []
        }

    /*
    -------------------------------------------------
        ASSEMBLY EXTENSION
    -------------------------------------------------
    */
    ch_extension_fastqs = ch_hostremoved_fastq_gz
        .filter { meta, reads -> !meta.single_end }

    ch_extension_fastas = input_fastas
        .filter { meta, reads -> ( meta.assembler == 'megahit' || meta.assembler == 'spades' ) && !meta.single_end}

    if (parameters.run_cobra) {
        //
        // MODULE: Identify sequences longer than 10kb
        //
        COBRA_SEQ(
            ch_extension_fastas
        )
        ch_versions = ch_versions.mix(COBRA_SEQ.out.versions)

        //
        // MODULES: Convert fasta to tabular format
        //
        COBRA_FX2TAB(
            COBRA_SEQ.out.fastx
        )
        ch_versions         = ch_versions.mix(COBRA_FX2TAB.out.versions)
        ch_query_contigs    = COBRA_FX2TAB.out.text
    } else {
        ch_query_contigs    = []
    }

    //
    // SUBWORKFLOW: Extend assemblies with read overlaps and graph resolution
    //
    FASTQFASTAGFA_ASSEMBLYEXTENSION_FASTA(
        ch_extension_fastqs,
        parameters.run_cobra,
        ch_extension_fastas,
        ch_assembly_logs,
        ch_query_contigs,
        parameters.run_phables,
        ch_assembly_graph_gz,
        file("${projectDir}/assets/configs/phables/phables_config.yml", checkIfExists:true),
        parameters.phables_db
    )
    ch_versions     = ch_versions.mix(FASTQFASTAGFA_ASSEMBLYEXTENSION_FASTA.out.versions)
    input_fastas    = input_fastas.map { meta, fasta -> [ meta + [ extension: 'no_extension' ], fasta ] }
    input_fastas    = input_fastas.mix(FASTQFASTAGFA_ASSEMBLYEXTENSION_FASTA.out.fasta_gz)

    /*
    -------------------------------------------------
        ASSEMBLY QC
    -------------------------------------------------
    */
    //
    // SUBWORKFLOW: Assess the quality of assemblies
    //
    FASTA_ASSEMBLYQC_FASTA(
        input_fastas,
        parameters.assembly_min_len,
        parameters.run_seqkit_stats,
        parameters.run_trfinder,
        parameters.use_trfinder_fasta,
        parameters.run_tantan,
        parameters.run_sequence_stats
    )
    ch_versions                     = ch_versions.mix(FASTA_ASSEMBLYQC_FASTA.out.versions)
    ch_filtered_assembly_fasta_gz   = FASTA_ASSEMBLYQC_FASTA.out.filtered_fasta_gz
    ch_filtered_assembly_faa_gz     = FASTA_ASSEMBLYQC_FASTA.out.pyrodigalgv_faa_gz
    ch_workdirs_to_clean            = ch_workdirs_to_clean.mix(FASTA_ASSEMBLYQC_FASTA.out.workdirs_to_clean)


    /*
    -------------------------------------------------
        MGE IDENTIFICATION
    -------------------------------------------------
    */
    //
    // SUBWORKFLOW: identify MGEs in assemblies
    //
    FASTA_MGECLASSIFICATION_FASTATSV(
        ch_filtered_assembly_fasta_gz,
        ch_filtered_assembly_faa_gz,
        parameters.run_genomad,
        parameters.genomad_db,
        parameters.run_pyhmmer_virus,
        parameters.run_pyhmmer_plasmid,
        parameters.run_pyhmmer_busco
    )
    ch_versions                 = ch_versions.mix(FASTA_MGECLASSIFICATION_FASTATSV.out.versions)

    /*
    -------------------------------------------------
        VIRUS COMPLETENESS
    -------------------------------------------------
    */
    if (parameters.run_checkv) {
        FASTA_CHECKV_FASTATSV(
            ch_filtered_assembly_fasta_gz,
            parameters.checkv_db,
            "${projectDir}/assets/db/ncbi_info.tsv"
        )
        ch_versions                 = ch_versions.mix(FASTA_CHECKV_FASTATSV.out.versions)
        ch_checkv_contamination_tsv = FASTA_CHECKV_FASTATSV.out.contamination_tsv
        ch_checkv_completeness_tsv  = FASTA_CHECKV_FASTATSV.out.completeness_tsv
        ch_checkv_genbank_hits_tsv  = FASTA_CHECKV_FASTATSV.out.genbank_hits_tsv
    } else {
        ch_checkv_contamination_tsv = []
        ch_checkv_completeness_tsv  = []
        ch_checkv_genbank_hits_tsv  = []
    }


    /*
    -------------------------------------------------
        VIRUS ANNOTATION
    -------------------------------------------------
    */

    /*
    -------------------------------------------------
        INTERMEDIATE FILE CLEANUP
    -------------------------------------------------
    */
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

    /*
    -------------------------------------------------
        REPORT GENERATION
    -------------------------------------------------
    */
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
    multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions        = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
