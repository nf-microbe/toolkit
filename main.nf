#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT PLUGINS / FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// PLUGINS
include { samplesheetToList                     } from 'plugin/nf-schema'

// FUNCTIONS
include { validateInputSamplesheet              } from './subworkflows/local/utils_nfmicrobe_toolkit_pipeline'
include { rmEmptyFastAs; rmEmptyFastQs          } from './subworkflows/nf-core/utils_nfmicrobe_functions'

// MODULES
include { CAT_FASTQ as CAT_FASTQ_RUNMERGE       } from './modules/nf-core/cat/fastq'
include { CAT_FASTQ as CAT_FASTQ_COASSEMBLY     } from './modules/nf-core/cat/fastq'
// include { CSVTK_CONCAT                          } from './modules/nf-core/csvtk/concat'
include { FASTP                                         } from './modules/nf-core/fastp'
include { LOGAN_CONTIGAWSCLIMULTIPLIER                  } from './modules/nf-core/logan/contigawsclimultiplier'
include { LOGAN_UNITIGAWSCLIMULTIPLIER                  } from './modules/nf-core/logan/unitigawsclimultiplier'
include { MEGAHIT as MEGAHIT_COASSEMBLY                 } from './modules/nf-core/megahit'
include { MEGAHIT as MEGAHIT_SINGLE                     } from './modules/nf-core/megahit'
include { PLASS_PENGUIN as PLASS_PENGUIN_COASSEMBLY     } from './modules/nf-core/plass/penguin'
include { PLASS_PENGUIN as PLASS_PENGUIN_SINGLE         } from './modules/nf-core/plass/penguin'
// include { PYHMMER_BUSCOCLASSIFY                 } from './modules/nf-core/pyhmmer/buscoclassify'
// include { PYHMMER_PLASMIDCLASSIFY               } from './modules/nf-core/pyhmmer/plasmidclassify'
// include { PYHMMER_VIRUSCLASSIFY                 } from './modules/nf-core/pyhmmer/virusclassify'
// include { SEQKIT_REPLACE                        } from './modules/nf-core/seqkit/replace'
// include { SEQKIT_SEQ                            } from './modules/nf-core/seqkit/seq'
// include { SEQKIT_SPLIT2                         } from './modules/nf-core/seqkit/split2'
// include { SEQKIT_STATS                          } from './modules/nf-core/seqkit/stats'
// include { SEQUENCESTATS                         } from './modules/nf-core/sequencestats'
include { SPADES as SPADES_COASSEMBLY                   } from './modules/nf-core/spades'
include { SPADES as SPADES_SINGLE                       } from './modules/nf-core/spades'
include { SRA_SRATOOLS                                  } from './modules/nf-core/sra/sratools'
include { SRA_ASPERACLI                                 } from './modules/nf-core/sra/asperacli'
// include { TANTAN                                } from './modules/nf-core/tantan'
// include { TRFINDER                              } from './modules/nf-core/trfinder'
// include { TRTRIMMER                             } from './modules/nf-core/trtrimmer'
// include { VTDB_COMBINEDATA                      } from './modules/local/vtdb/combinedata'
// include { VTDB_CLASSIFICATIONFILTER             } from './modules/local/vtdb/classificationfilter'
// include { VTDB_COMPLETENESSFILTER               } from './modules/local/vtdb/completenessfilter'
// include { VTDB_COMPOSITIONFILTER                } from './modules/local/vtdb/compositionfilter'
// include { VTDB_FILTERSEQUENCES                  } from './modules/local/vtdb/filtersequences'

// SUBWORKFLOWS
// include { FASTA_CHECKV_TSV                      } from './subworkflows/nf-core/fasta_checkv_tsv'
// include { FASTA_GENOMAD_FAATSV                  } from './subworkflows/nf-core/fasta_genomad_faatsv'
// include { FASTA_SEQHASHER_FASTA                 } from './subworkflows/nf-core/fasta_seqhasher_fasta'
// include { FASTA_VCLUST_FASTA                    } from './subworkflows/nf-core/fasta_vclust_fasta'
include { FASTQ_BOWTIE2_FASTQ                   } from './subworkflows/nf-core/fastq_bowtie2_fastq'
include { FASTQ_VIROMEQC_TSV                    } from './subworkflows/nf-core/fastq_viromeqc_tsv'
include { FASTQFASTA_COBRA_FASTA                } from './subworkflows/nf-core/fastqfasta_cobra_fasta'
// include { FASTQGFA_PHABLES_FASTA                } from './subworkflows/nf-core/fastqgfa_phables_fasta'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    main:

    ch_versions             = Channel.empty()
    ch_input_fastqs_prefilt = Channel.empty()
    ch_input_fastas_prefilt = Channel.empty()

    /*
    -------------------------------------------------
        LOAD INPUTS
    -------------------------------------------------
    */
    if (params.input) {
        //
        // Create channel from input file provided through params.input
        //
        ch_samplesheet = Channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2, fasta ->
                    meta.run        = meta.run == null ? 1 : meta.run
                    meta.group      = meta.group == null ? meta.id : meta.group
                    meta.single_end = fastq_2 ? false : true
                    if (meta.single_end) {
                        return [ meta, [ fastq_1 ], fasta ]
                    } else {
                        return [ meta, [ fastq_1, fastq_2 ], fasta ]
                    }
            }
            .multiMap { meta, fastqs, fastas ->
                fastqs: [ meta, fastqs ]
                fastas: [ meta, fastas ]
            }

        // validate samplesheet
        ch_samplesheet.fastqs
            .map { meta, fastq ->
                [ meta.id, meta, fastq ]
            }
            .groupTuple()
            .map { samplesheet -> validateInputSamplesheet(samplesheet) }

        ch_input_fastqs_prefilt = ch_input_fastqs_prefilt.mix(ch_samplesheet.fastqs)
        ch_input_fastas_prefilt = ch_input_fastas_prefilt.mix(ch_samplesheet.fastas)
    }

    if (params.fastqs) {
        // Prepare FastQ channel
        ch_input_fastqs_prefilt = ch_input_fastqs_prefilt.mix(Channel
            .fromFilePairs(params.fastqs, size:-1)
            .ifEmpty { error("Cannot find any reads matching: ${params.fastqs}\nNB: Path needs to be enclosed in quotes!") }
            .map { meta, fastq ->
                def meta_new = [:]
                meta_new.id           = meta
                meta_new.run          = 1
                meta_new.group        = meta
                meta_new.single_end   = fastq.size() == 1 ? true : false
                if ( meta_new.single_end ) {
                    return [ meta_new, [ fastq[0] ] ]
                } else {
                    return [ meta_new, [ fastq[0], fastq[1] ] ]
                }
            }
        )
    }

    if (params.fastas) {
        // Prepare FastA channel
        ch_input_fastas_prefilt = ch_input_fastas_prefilt.mix(Channel
            .fromPath(params.fastas)
            .ifEmpty { error("Cannot find any FastA files matching: ${params.fastas}\nNB: Path needs to be enclosed in quotes!") }
            .map { fasta ->
                def meta_new = [:]
                meta_new.id                 = fasta.getSimpleName()
                meta_new.run                = 1
                meta_new.group              = fasta.getSimpleName()
                return [ meta_new, fasta ]
            }
        )
    }

    // Filter out empty input channels
    ch_input_fastqs = ch_input_fastqs_prefilt.filter { meta, fastqs -> fastqs[0] }
    ch_input_fastas = ch_input_fastas_prefilt.filter { meta, fastas -> fastas }

    ch_sra_accessions       = null

    if (params.sra_accessions) {
        // Load SRA accessions one-by-one from file
        ch_sra_accessions = Channel.fromPath(params.sra_accessions)
            .splitCsv(header:false)
            .map { row -> [ [ id: row[0], run: 1, group: row[0] ], row[0] ] }
    }

    /*
    -------------------------------------------------
        READ DOWNLOAD
    -------------------------------------------------
    */
    if (ch_sra_accessions && params.sra_download_method == 'sratools') {
        //
        // MODULE: Download SRA FastQ files
        //
        SRA_SRATOOLS(
            ch_sra_accessions.map { meta, acc -> [ [ id: "sra_" + meta.id, run: meta.run, group: "sra_" + meta.id ], acc ] }
        )
        ch_versions     = ch_versions.mix(SRA_SRATOOLS.out.versions)
        ch_input_fastqs = ch_input_fastqs
            .mix(
                SRA_SRATOOLS.out.fastq.map { meta, fastq ->
                    meta.single_end = fastq.size() == 2 ? false : true
                    [ meta, fastq ]
                }
            )
    }

    if (ch_sra_accessions && params.sra_download_method == 'asperacli') {
        //
        // MODULE: Download SRA FastQ files
        //
        SRA_ASPERACLI(
            ch_sra_accessions.map { meta, acc -> [ [ id: "sra_" + meta.id, run: meta.run, group: "sra_" + meta.id ], acc ] }
        )
        ch_versions     = ch_versions.mix(SRA_ASPERACLI.out.versions)
        ch_input_fastqs = ch_input_fastqs
            .mix(
                SRA_ASPERACLI.out.fastq.map { meta, fastq ->
                    meta.single_end = fastq.size() == 2 ? false : true
                    [ meta, fastq ]
                }
            )
    }

    /*
    -------------------------------------------------
        READ PREPROCESSING
    -------------------------------------------------
    */
    if (params.run_fastp) {
        //
        // MODULE: Run fastp on raw reads
        //
        FASTP(
            ch_input_fastqs,
            [],
            false,
            false,
            false
        )
        ch_fastp_prefilt_fastq_gz   = FASTP.out.reads
        ch_versions                 = ch_versions.mix(FASTP.out.versions)

        // REMOVE EMPTY FASTQ FILES
        ch_fastp_fastq_gz   = rmEmptyFastQs(ch_fastp_prefilt_fastq_gz, false)
    } else {
        ch_fastp_fastq_gz   = ch_input_fastqs
    }

    /*
    -------------------------------------------------
        HOST READ REMOVAL
    -------------------------------------------------
    */
    if (params.run_bowtie2_host_removal) {
        // Load igenomes fasta and bowtie2 index
        igenomes_fasta   = params.genomes[params.igenomes_host_key]['fasta']
        igenomes_index   = params.genomes[params.igenomes_host_key]['bowtie2']

        // create host fasta and index channels based on input parameters
        if (params.igenomes_host_key) {
            ch_host_fasta_gz    = Channel.value([ [ id:'bowtie2_fasta' ], file(igenomes_fasta, checkIfExists: true) ])
            ch_bowtie2_index    = Channel.value([ [ id:'bowtie2_fasta' ], file(igenomes_index, checkIfExists: true) ])
        } else if (params.host_bowtie2_index) {
            if (!params.host_fasta) {
                error "[nf-microbe/toolkit]: --run_bowtie2_host_removal = true but no iGenomes (--igenomes_host_key) or custom fasta (--host_fasta) provided"
            }
            ch_host_fasta_gz    = Channel.value([ [ id:'bowtie2_fasta' ], file(params.host_fasta, checkIfExists: true) ])
            ch_bowtie2_index    = Channel.value([ [ id:'bowtie2_fasta' ], file(params.host_bowtie2_index, checkIfExists: true) ])
        } else {
            if (!params.host_fasta) {
                error "[nf-microbe/toolkit]: --run_bowtie2_host_removal = true but no iGenomes (--igenomes_host_key) or custom fasta (--host_fasta) provided"
            }
            ch_host_fasta_gz    = Channel.value([ [ id:'bowtie2_fasta' ], file(params.host_fasta, checkIfExists: true) ])
            ch_bowtie2_index    = null
        }

        //
        // SUBWORKFLOW: Remove host reads using Bowtie2
        //
        FASTQ_BOWTIE2_FASTQ(
            ch_fastp_fastq_gz,
            ch_host_fasta_gz,
            ch_bowtie2_index
        )
        ch_hostremoved_prefilt_fastq_gz = FASTQ_BOWTIE2_FASTQ.out.fastq_gz
        ch_versions                     = ch_versions.mix(FASTQ_BOWTIE2_FASTQ.out.versions)

        // REMOVE EMPTY FASTQ FILES FROM CHANNEL
        ch_hostremoved_fastq_gz = rmEmptyFastQs(ch_hostremoved_prefilt_fastq_gz, false)
    } else {
        ch_hostremoved_fastq_gz = ch_fastp_fastq_gz
    }

    //
    // TODO: MODULE: Count reads after host removal
    //

    /*
    -------------------------------------------------
        VIRUS ENRICHMENT ESTIMATION
    -------------------------------------------------
    */
    if (params.run_viromeqc) {
        //
        // SUBWORKFLOW: Estimate viral enrichment
        //
        FASTQ_VIROMEQC_TSV(
            ch_hostremoved_fastq_gz,
            params.viromeqc_db
        )
        ch_versions = ch_versions.mix(FASTQ_VIROMEQC_TSV.out.versions)
    }

    /*
    -------------------------------------------------
        RUN MERGING
    -------------------------------------------------
    */
    if (params.perform_run_merging) {
        // prepare reads for concatenating within runs
        ch_reads_forcat = ch_hostremoved_fastq_gz
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
        CAT_FASTQ_RUNMERGE(
            ch_reads_forcat.cat.map { meta, reads -> [ meta + [ run: 'merged' ], reads.flatten() ] }
        )
        ch_runmerged_fastq_gz   = CAT_FASTQ_RUNMERGE.out.reads
        ch_versions             = ch_versions.mix(CAT_FASTQ_RUNMERGE.out.versions)

        // Ensure we don't have nests of nests so that structure is in form expected for assembly
        ch_reads_forcat_skipped = ch_reads_forcat.skip_cat
            .map { meta, reads ->
                def new_reads = meta.single_end ? reads[0] : reads.flatten()
                [ meta, new_reads ]
            }

        // Combine single run and multi-run-merged data
        ch_runmerged_fastq_gz = ch_runmerged_fastq_gz.mix(ch_reads_forcat_skipped.map { meta, reads -> [ meta + [ run: 1 ], reads ] })
    } else {
        ch_runmerged_fastq_gz = ch_hostremoved_fastq_gz
    }

    //
    // TODO: MODULE: Count reads after run merging
    //

    /*
    -------------------------------------------------
        LOGAN CONTIG/UNITIG DOWNLOAD
    -------------------------------------------------
    */
    if (ch_sra_accessions && (params.download_logan_contigs)) {
        //
        // MODULE: Download logan contigs
        //
        LOGAN_CONTIGAWSCLIMULTIPLIER(
            ch_sra_accessions.map { meta, accession -> [ [ id: "logan_contig_" + meta.id, run: 1, group: "logan_contig_" + meta.id, single_end: true ], accession] }
        )
        ch_versions                     = ch_versions.mix(LOGAN_CONTIGAWSCLIMULTIPLIER.out.versions.first())
        ch_logan_contig_fa_gz           = LOGAN_CONTIGAWSCLIMULTIPLIER.out.raw_fasta.map { meta, fasta -> [ meta + [ assembler: 'minia', assembly_method: 'single' ], fasta ] }
        ch_logan_contig_filt_fa_gz      = LOGAN_CONTIGAWSCLIMULTIPLIER.out.filtered_fasta
        ch_logan_contig_mult_fa_gz      = LOGAN_CONTIGAWSCLIMULTIPLIER.out.multiplied_fasta
        ch_logan_contig_filt_mult_fa_gz = LOGAN_CONTIGAWSCLIMULTIPLIER.out.filt_mult_fasta

        ch_input_fastas = ch_input_fastas.mix(ch_logan_contig_fa_gz)

        ch_runmerged_fastq_gz   = ch_runmerged_fastq_gz
            .mix(ch_logan_contig_filt_fa_gz)
            .mix(ch_logan_contig_mult_fa_gz)
            .mix(ch_logan_contig_filt_mult_fa_gz)

        if (params.assemble_logan_contigs) {
            ch_runmerged_fastq_gz  = ch_runmerged_fastq_gz
                .mix(ch_logan_contig_fa_gz)
        }
    }

    if (ch_sra_accessions && (params.download_logan_unitigs)) {
        //
        // MODULE: Download logan contigs
        //
        LOGAN_UNITIGAWSCLIMULTIPLIER(
            ch_sra_accessions.map { meta, accession -> [ [ id: "logan_unitig_" + meta.id, run: 1, group: "logan_unitig_" + meta.id, single_end: true ], accession] }
        )
        ch_versions                     = ch_versions.mix(LOGAN_UNITIGAWSCLIMULTIPLIER.out.versions.first())
        ch_logan_unitig_fa_gz           = LOGAN_UNITIGAWSCLIMULTIPLIER.out.raw_fasta
        ch_logan_unitig_filt_fa_gz      = LOGAN_UNITIGAWSCLIMULTIPLIER.out.filtered_fasta
        ch_logan_unitig_mult_fa_gz      = LOGAN_UNITIGAWSCLIMULTIPLIER.out.multiplied_fasta
        ch_logan_unitig_filt_mult_fa_gz = LOGAN_UNITIGAWSCLIMULTIPLIER.out.filt_mult_fasta

        ch_runmerged_fastq_gz   = ch_runmerged_fastq_gz
            .mix(ch_logan_unitig_fa_gz)
            .mix(ch_logan_unitig_filt_fa_gz)
            .mix(ch_logan_unitig_mult_fa_gz)
            .mix(ch_logan_unitig_filt_mult_fa_gz)
    }

    /*
    -------------------------------------------------
        ASSEMBLY DOWNLOAD
    -------------------------------------------------
    */


    /*
    -------------------------------------------------
        READ ASSEMBLY
    -------------------------------------------------
    */
    ch_coassembly_fastq_gz          = Channel.empty()
    ch_assemblies_prefilt_fasta_gz  = Channel.empty()
    ch_assembly_graph_gfa_gz        = Channel.empty()
    ch_assembly_logs                = Channel.empty()

    if (params.run_megahit_coassembly || params.run_spades_coassembly || params.run_penguin_coassembly) {
        // group and set group as new id
        ch_cat_coassembly_fastq_gz = ch_runmerged_fastq_gz
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple( by: 0, sort:'deep' )
            .map { group, meta, reads ->
                def meta_new                = [:]
                meta_new.id                 = "group-$group"
                meta_new.run                = 1
                meta_new.group              = group
                meta_new.single_end         = meta.single_end[0]
                if ( meta_new.single_end ) {
                    return [ meta_new, reads.collect { fastqs -> fastqs } ]
                } else {
                    return [ meta_new, reads.flatten() ]
                }
            }
            .branch {
                meta, reads ->
                    coassembly: meta.single_end && reads.size() >= 2 || !meta.single_end && reads.size() >= 4
                    skip_coassembly: true   // Can skip coassembly if there is not multiple samples
            }

        //
        // MODULE: Combine reads within groups for coassembly
        //
        CAT_FASTQ_COASSEMBLY(
            ch_cat_coassembly_fastq_gz.coassembly
        )
        ch_coassembly_fastq_gz  = CAT_FASTQ_COASSEMBLY.out.reads
        ch_versions             = ch_versions.mix(CAT_FASTQ_COASSEMBLY.out.versions)
    }

    if (params.run_megahit_single) {
        //
        // MODULE: Assemble reads individually with MEGAHIT
        //
        MEGAHIT_SINGLE(
            ch_runmerged_fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'megahit', assembly_method: 'single'], fasta ] },
            "${projectDir}/bin/fastg2gfa"
        ).contigs
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(MEGAHIT_SINGLE.out.contigs)
        ch_assembly_graph_gfa_gz        = ch_assembly_graph_gfa_gz.mix(MEGAHIT_SINGLE.out.gfa)
        ch_assembly_logs                = ch_assembly_logs.mix(MEGAHIT_SINGLE.out.log)
        ch_versions                     = ch_versions.mix(MEGAHIT_SINGLE.out.versions)
    }

    if (params.run_megahit_coassembly) {
        //
        // MODULE: Co-assemble reads with MEGAHIT
        //
        MEGAHIT_COASSEMBLY(
            ch_coassembly_fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'megahit', assembly_method: 'coassembly' ], fasta ] },
            "${projectDir}/bin/fastg2gfa"
        ).contigs
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(MEGAHIT_COASSEMBLY.out.contigs)
        ch_assembly_graph_gfa_gz        = ch_assembly_graph_gfa_gz.mix(MEGAHIT_COASSEMBLY.out.gfa)
        ch_assembly_logs                = ch_assembly_logs.mix(MEGAHIT_COASSEMBLY.out.log)
        ch_versions                     = ch_versions.mix(MEGAHIT_COASSEMBLY.out.versions)
    }

    if (params.run_spades_single) {
        // prepare reads for metaspades input
        ch_spades_single_input = ch_runmerged_fastq_gz
            .map { meta, fastq ->
                [ meta + [ assembler: 'spades', assembly_method: 'single' ], fastq, [], [] ]
            }
            .branch { meta, fastq, extra1, extra2 ->
                single_end: meta.single_end
                paired_end: true
            }

        //
        // MODULE: Assemble reads individually with metaSPAdes
        //
        SPADES_SINGLE(
            ch_spades_single_input.paired_end,
            [],
            []
        )
        if (params.use_spades_scaffolds) {
            ch_spades_single_fasta_gz   = SPADES_SINGLE.out.scaffolds
        } else {
            ch_spades_single_fasta_gz   = SPADES_SINGLE.out.contigs
        }
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_spades_single_fasta_gz)
        ch_assembly_graph_gfa_gz        = ch_assembly_graph_gfa_gz.mix(SPADES_SINGLE.out.gfa)
        ch_assembly_logs                = ch_assembly_logs.mix(SPADES_SINGLE.out.log)
        ch_versions                     = ch_versions.mix(SPADES_SINGLE.out.versions)
    }

    if (params.run_spades_coassembly) {
        // prepare reads for metaspades input
        ch_metaspades_coassembly_input = ch_coassembly_fastq_gz
            .map { meta, fastq ->
                [ meta + [ assembler: 'spades', assembly_method: 'coassembly' ], fastq, [], [] ]
            }
            .branch { meta, fastq, extra1, extra2 ->
                single_end: meta.single_end
                paired_end: true
            }

        //
        // MODULE: Co-assemble reads with metaSPAdes
        //
        SPADES_COASSEMBLY(
            ch_metaspades_coassembly_input.paired_end,
            [],
            []
        )
        if (params.use_spades_scaffolds) {
            ch_spades_co_fasta_gz   = SPADES_COASSEMBLY.out.scaffolds
        } else {
            ch_spades_co_fasta_gz   = SPADES_COASSEMBLY.out.contigs
        }
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_spades_co_fasta_gz)
        ch_assembly_graph_gfa_gz        = ch_assembly_graph_gfa_gz.mix(SPADES_COASSEMBLY.out.gfa)
        ch_assembly_logs                = ch_assembly_logs.mix(SPADES_COASSEMBLY.out.log)
        ch_versions                     = ch_versions.mix(SPADES_COASSEMBLY.out.versions)
    }

    if (params.run_penguin_single) {
        //
        // MODULE: Assemble reads individually with PenguiN
        //
        PLASS_PENGUIN_SINGLE(
            ch_runmerged_fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'penguin', assembly_method: 'single' ], fasta ] }
        )
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(PLASS_PENGUIN_SINGLE.out.contigs)
        ch_versions                     = ch_versions.mix(PLASS_PENGUIN_SINGLE.out.versions)
    }

    if (params.run_penguin_coassembly) {
        //
        // MODULE: Co-assemble reads with PenguiN
        //
        PLASS_PENGUIN_COASSEMBLY(
            ch_cat_coassembly_fastq_gz.coassembly.map { meta, fastq -> [ meta + [ assembler: 'penguin', assembly_method: 'coassembly' ], fastq ] }
        )
        ch_versions                     = ch_versions.mix(PLASS_PENGUIN_COASSEMBLY.out.versions)
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(PLASS_PENGUIN_COASSEMBLY.out.contigs)
    }

    // REMOVE EMPTY FASTA FILES FROM CHANNEL
    ch_assemblies_fasta_gz  = rmEmptyFastAs(ch_assemblies_prefilt_fasta_gz.mix(ch_input_fastas), false).map {
        meta, fasta -> [ meta, fasta ]
    }

    /*
    -------------------------------------------------
        ASSEMBLY EXTENSION
    -------------------------------------------------
    */
    // Combine single run and multi-run-merged data
    ch_runmerged_fastq_gz = ch_runmerged_fastq_gz.mix(ch_coassembly_fastq_gz)

    // Identify FastQ and FastA files eligible for extension
    ch_extension_fastqs = ch_runmerged_fastq_gz.filter { meta, reads -> !meta.single_end }
    ch_extension_fastas = ch_assemblies_fasta_gz.filter { meta, reads -> ( meta.assembler == 'megahit' || meta.assembler == 'spades' ) && !meta.single_end}

    if (params.run_cobra) {
        //
        // SUBWORKFLOW: Extend assemblies with cobra
        //
        FASTQFASTA_COBRA_FASTA(
            ch_extension_fastqs,
            ch_extension_fastas.map{ meta, fasta -> [ meta + [ extension: 'cobra' ], fasta ] },
            ch_assembly_logs.map{ meta, fasta -> [ meta + [ extension: 'cobra' ], fasta ] },
        )
        ch_assemblies_fasta_gz  = ch_assemblies_fasta_gz.mix(rmEmptyFastAs(FASTQFASTA_COBRA_FASTA.out.fasta_gz, false))
        ch_versions             = ch_versions.mix(FASTQFASTA_COBRA_FASTA.out.versions)
    }

    // if (params.run_phables) {
    //     //
    //     // SUBWORKFLOW: Resolve genome graphs with phables
    //     //
    //     FASTQGFA_PHABLES_FASTA(
    //         ch_extension_fastqs,
    //         rmEmptyFastAs(ch_assembly_graph_gfa_gz, false).map{ meta, fastq -> [ meta + [ extension: 'phables' ], fastq ] },
    //         "${projectDir}/assets/configs/phables/phables_config.yml",
    //         params.phables_db
    //     )
    //     ch_assemblies_fasta_gz  = ch_assemblies_fasta_gz.mix(rmEmptyFastAs(FASTQGFA_PHABLES_FASTA.out.fasta_gz, false))
    //     ch_versions             = ch_versions.mix(FASTQGFA_PHABLES_FASTA.out.versions)
    // }

    // /*
    // -------------------------------------------------
    //     ASSEMBLY QC/FILTERING
    // -------------------------------------------------
    // */
    // if (params.assembly_min_len > 0) {
    //     //
    //     // MODULE: Filter assemblies by length
    //     //
    //     SEQKIT_SEQ(
    //         ch_assemblies_fasta_gz
    //     )
    //     ch_seqkit_seq_prefilt_fasta_gz  = SEQKIT_SEQ.out.fastx
    //     ch_versions                     = ch_versions.mix(SEQKIT_SEQ.out.versions)

    //     // REMOVE EMPTY FASTA FILES FROM CHANNEL
    //     ch_seqkit_seq_fasta_gz = rmEmptyFastAs(ch_seqkit_seq_prefilt_fasta_gz, false)
    // } else {
    //     ch_seqkit_seq_fasta_gz = ch_assemblies_fasta_gz
    // }

    // if (params.run_seqkit_stats) {
    //     //
    //     // MODULE: Calculate assembly statistics
    //     //
    //     SEQKIT_STATS(
    //         ch_seqkit_seq_fasta_gz
    //     )
    //     ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions)
    // }

    // if (params.prepend_sample_id) {
    //     //
    //     // MODULE: Prepend assembly metadata to assembly sequences
    //     //
    //     SEQKIT_REPLACE(
    //         ch_seqkit_seq_fasta_gz
    //     )
    //     ch_versions                     = ch_versions.mix(SEQKIT_REPLACE.out.versions)
    //     ch_renamed_assembly_fasta_gz    = SEQKIT_REPLACE.out.fastx
    // } else {
    //     ch_renamed_assembly_fasta_gz    = ch_seqkit_seq_fasta_gz
    // }

    // /*
    // -------------------------------------------------
    //     ASSEMBLY SPLIT
    // -------------------------------------------------
    // */
    // if (params.assembly_split_size > 0) {
    //     //
    //     // MODULE: Split assemblies into chunks for downstream processing
    //     //
    //     SEQKIT_SPLIT2(
    //         ch_renamed_assembly_fasta_gz
    //     )
    //     ch_versions         = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
    //     ch_split_fasta_gz   = SEQKIT_SPLIT2.out.reads
    //         .map{ meta, fasta ->
    //             fasta.collect { file ->
    //                 return [ meta + [ chunk: file.getBaseName().tokenize('.')[1]], file ]
    //             }
    //         }
    //         .flatMap { chunks -> chunks }
    // } else {
    //     ch_split_fasta_gz   = ch_renamed_assembly_fasta_gz.map { meta, fasta -> [ meta + [ chunk: '1' ], fasta ] }
    // }

    // if (params.run_trfinder) {
    //     //
    //     // MODULE: Identify sequences with terminal repeats
    //     //
    //     TRFINDER(
    //         ch_split_fasta_gz
    //     )
    //     ch_trfinder_prefilt_fasta_gz    = TRFINDER.out.fasta
    //     ch_trfinder_tsv                 = TRFINDER.out.stats
    //     ch_versions                     = ch_versions.mix(TRFINDER.out.versions.first())

    //     if (params.use_trfinder_fasta) {
    //         // REMOVE EMPTY FASTA FILES FROM CHANNEL
    //         ch_trfinder_fasta_gz = rmEmptyFastAs(ch_trfinder_prefilt_fasta_gz, false)
    //     } else {
    //         ch_trfinder_fasta_gz = ch_seqkit_seq_fasta_gz
    //     }
    // } else {
    //     ch_trfinder_fasta_gz    = ch_split_fasta_gz
    //     ch_trfinder_tsv         = []
    // }

    // /*
    // -------------------------------------------------
    //     MGE IDENTIFICATION
    // -------------------------------------------------
    // */
    // if (params.run_genomad || params.run_pyhmmer_virus || params.run_pyhmmer_plasmid || params.run_pyhmmer_busco) {
    //     //
    //     // SUBWORKFLOW: Classify MGEs with geNomad
    //     //
    //     FASTA_GENOMAD_FAATSV(
    //         ch_split_fasta_gz,
    //         params.genomad_db
    //     )
    //     ch_versions                 = ch_versions.mix(FASTA_GENOMAD_FAATSV.out.versions.first())
    //     ch_genomad_genes_tsv        = FASTA_GENOMAD_FAATSV.out.genes_tsv
    //     ch_genomad_features_tsv     = FASTA_GENOMAD_FAATSV.out.features_tsv
    //     ch_genomad_scores_tsv       = FASTA_GENOMAD_FAATSV.out.scores_tsv
    //     ch_genomad_taxonomy_tsv     = FASTA_GENOMAD_FAATSV.out.taxonomy_tsv
    //     ch_genomad_proteins_faa_gz  = FASTA_GENOMAD_FAATSV.out.proteins_faa_gz
    // } else {
    //     ch_genomad_genes_tsv    = []
    //     ch_genomad_features_tsv = []
    //     ch_genomad_scores_tsv   = []
    //     ch_genomad_taxonomy_tsv = []
    // }

    // if (params.run_pyhmmer_busco) {
    //     //
    //     // MODULE: Classify MGEs with pyHMMER
    //     //
    //     PYHMMER_BUSCOCLASSIFY(
    //         ch_genomad_proteins_faa_gz,
    //         "${projectDir}/assets/hmms/busco_hmms/archaea_buscos.hmm",
    //         "${projectDir}/assets/hmms/busco_hmms/score_cutoffs/archaea_odb10.cutoffs",
    //         "${projectDir}/assets/hmms/busco_hmms/bacteria_buscos.hmm",
    //         "${projectDir}/assets/hmms/busco_hmms/score_cutoffs/bacteria_odb10.cutoffs"
    //     )
    //     ch_versions         = ch_versions.mix(PYHMMER_BUSCOCLASSIFY.out.versions)
    //     ch_busco_hmms_tsv   = PYHMMER_BUSCOCLASSIFY.out.markers
    // } else {
    //     ch_busco_hmms_tsv   = []
    // }

    // if (params.run_pyhmmer_plasmid) {
    //     //
    //     // MODULE: Classify MGEs with pyHMMER
    //     //
    //     PYHMMER_PLASMIDCLASSIFY(
    //         ch_genomad_proteins_faa_gz,
    //         "${projectDir}/assets/hmms/plasmid_hmms/plasmid_hallmarks.hmm"
    //     )
    //     ch_versions         = ch_versions.mix(PYHMMER_PLASMIDCLASSIFY.out.versions)
    //     ch_plasmid_hmms_tsv = PYHMMER_PLASMIDCLASSIFY.out.markers
    // } else {
    //     ch_plasmid_hmms_tsv = []
    // }

    // if (params.run_pyhmmer_virus) {
    //     //
    //     // MODULE: Classify MGEs with pyHMMER
    //     //
    //     PYHMMER_VIRUSCLASSIFY(
    //         ch_genomad_proteins_faa_gz,
    //         "${projectDir}/assets/hmms/virus_hmms/DJR_MCP_virus_hallmarks.hmm",
    //         "${projectDir}/assets/hmms/virus_hmms/inovirus_MCP_virus_hallmarks.hmm",
    //         "${projectDir}/assets/hmms/virus_hmms/pleolipoviridae_virus_hallmarks.hmm"
    //     )
    //     ch_versions         = ch_versions.mix(PYHMMER_VIRUSCLASSIFY.out.versions)
    //     ch_virus_hmms_tsv   = PYHMMER_VIRUSCLASSIFY.out.markers
    // } else {
    //     ch_virus_hmms_tsv   = []
    // }

    // /*
    // -------------------------------------------------
    //     VIRUS COMPLETENESS
    // -------------------------------------------------
    // */
    // if (params.run_checkv) {
    //     FASTA_CHECKV_TSV(
    //         ch_split_fasta_gz,
    //         params.checkv_db,
    //         "${projectDir}/assets/db/ncbi_info.tsv"
    //     )
    //     ch_versions                 = ch_versions.mix(FASTA_CHECKV_TSV.out.versions)
    //     ch_checkv_contamination_tsv = FASTA_CHECKV_TSV.out.contamination_tsv
    //     ch_checkv_completeness_tsv  = FASTA_CHECKV_TSV.out.completeness_tsv
    //     ch_checkv_genbank_hits_tsv  = FASTA_CHECKV_TSV.out.genbank_hits_tsv
    // } else {
    //     ch_checkv_contamination_tsv = []
    //     ch_checkv_completeness_tsv  = []
    //     ch_checkv_genbank_hits_tsv  = []
    // }

    // /*
    // -------------------------------------------------
    //     MGE QC
    // -------------------------------------------------
    // */

    // if (params.run_tantan) {
    //     //
    //     // MODULE: identify low-complexity regions
    //     //
    //     TANTAN(
    //         ch_split_fasta_gz
    //     )
    //     ch_tantan_bed   = TANTAN.out.bed
    //     ch_versions     = ch_versions.mix(TANTAN.out.versions.first())
    // } else {
    //     ch_tantan_bed = []
    // }

    // if (params.run_sequence_stats) {
    //     // join inputs by meta.id for sequence stats
    //     ch_seq_stats_inputs = ch_split_fasta_gz
    //         .join(ch_genomad_proteins_faa_gz)
    //         .multiMap { meta, fasta, faa ->
    //             fasta:  [ meta, fasta ]
    //             faa:    [ meta, faa ]
    //         }

    //     //
    //     // MODULE: Calculate sequence statistics
    //     //
    //     SEQUENCESTATS(
    //         ch_seq_stats_inputs.fasta,
    //         ch_seq_stats_inputs.faa
    //     )
    //     ch_versions             = ch_versions.mix(SEQUENCESTATS.out.versions)
    //     ch_sequence_stats_tsv   = SEQUENCESTATS.out.stats
    // } else {
    //     ch_sequence_stats_tsv   = []
    // }

    // // /*
    // // -------------------------------------------------
    // //     MGE FILTERING
    // // -------------------------------------------------
    // // */

    // if (!params.combined_data &&
    //     (
    //         params.run_trfinder ||
    //         params.run_genomad ||
    //         params.run_pyhmmer_buscoclassify ||
    //         params.run_pyhmmer_plasmidclassify ||
    //         params.run_pyhmmer_virusclassify ||
    //         params.run_checkv ||
    //         params.run_tantan ||
    //         params.run_sequence_stats
    //     )
    // ) {
    //     // join all sequence data files for input
    //     ch_combinedata_input    = ch_split_fasta_gz
    //         .join(ch_trfinder_tsv)
    //         .join(ch_genomad_scores_tsv)
    //         .join(ch_genomad_genes_tsv)
    //         .join(ch_genomad_taxonomy_tsv)
    //         .join(ch_busco_hmms_tsv)
    //         .join(ch_plasmid_hmms_tsv)
    //         .join(ch_virus_hmms_tsv)
    //         .join(ch_checkv_completeness_tsv)
    //         .join(ch_checkv_contamination_tsv)
    //         .join(ch_tantan_bed)
    //         .join(ch_sequence_stats_tsv)
    //         .multiMap { it ->
    //             fasta:          [ it[0], it[1] ]
    //             trfinder:       [ it[0], it[2] ]
    //             genomad_scores: [ it[0], it[3] ]
    //             genomad_genes:  [ it[0], it[4] ]
    //             genomad_taxa:   [ it[0], it[5] ]
    //             busco_hmms:     [ it[0], it[6] ]
    //             plasmid_hmms:   [ it[0], it[7] ]
    //             virus_hmms:     [ it[0], it[8] ]
    //             completeness:   [ it[0], it[9] ]
    //             contamination:  [ it[0], it[10] ]
    //             tantan:         [ it[0], it[11] ]
    //             nuc_stats:      [ it[0], it[12] ]
    //         }

    //     //
    //     // MODULE: Combine all sequence data into one TSV file
    //     //
    //     VTDB_COMBINEDATA(
    //         ch_combinedata_input.fasta,
    //         ch_combinedata_input.trfinder,
    //         ch_combinedata_input.genomad_scores,
    //         ch_combinedata_input.genomad_genes,
    //         ch_combinedata_input.genomad_taxa,
    //         ch_combinedata_input.busco_hmms,
    //         ch_combinedata_input.plasmid_hmms,
    //         ch_combinedata_input.virus_hmms,
    //         ch_combinedata_input.completeness,
    //         ch_combinedata_input.contamination,
    //         ch_combinedata_input.tantan,
    //         ch_combinedata_input.nuc_stats
    //     )
    //     ch_combined_data_tsv    = VTDB_COMBINEDATA.out.tsv
    //     ch_versions             = ch_versions.mix(VTDB_COMBINEDATA.out.versions)
    // } else if (params.combined_data) {
    //     // load combined data file from params
    //     ch_combined_data_tsv    = [
    //         [ id:'all_samples' ],
    //         file(params.combined_data, checkIfExists:true)
    //     ]
    // } else {
    //     ch_combined_data_tsv    = Channel.empty()
    // }

    // if (params.classification_filters) {
    //     //
    //     // MODULE: Assess whether sequence passes viral classification filters
    //     //
    //     VTDB_CLASSIFICATIONFILTER(
    //         ch_combined_data_tsv,
    //         params.classification_filters
    //     )
    //     ch_class_data_tsv   = VTDB_CLASSIFICATIONFILTER.out.class_data
    //     ch_versions         = ch_versions.mix(VTDB_CLASSIFICATIONFILTER.out.versions)
    // } else {
    //     ch_class_data_tsv   = ch_combined_data_tsv
    // }

    // if (params.composition_filters) {
    //     //
    //     // MODULE: Assess whether sequence passes composition filters
    //     //
    //     VTDB_COMPOSITIONFILTER(
    //         ch_class_data_tsv,
    //         params.composition_filters
    //     )
    //     ch_compos_data_tsv  = VTDB_COMPOSITIONFILTER.out.compos_data
    //     ch_versions         = ch_versions.mix(VTDB_COMPOSITIONFILTER.out.versions)

    // } else {
    //     ch_compos_data_tsv  = ch_class_data_tsv
    // }

    // if (params.completeness_filters) {
    //     //
    //     // MODULE: Assess whether sequence passes completeness filters
    //     //
    //     VTDB_COMPLETENESSFILTER(
    //         ch_compos_data_tsv,
    //         params.completeness_filters
    //     )
    //     ch_compl_data_tsv   = VTDB_COMPLETENESSFILTER.out.compl_data
    //     ch_versions         = ch_versions.mix(VTDB_COMPLETENESSFILTER.out.versions)
    // } else {
    //     ch_compl_data_tsv   = ch_compos_data_tsv
    // }

    // // combine all completeness data files
    // ch_csvtk_concat_input   = ch_compl_data_tsv
    //     .map { meta, tsv -> [ [ id: "all_samples" ], tsv ] }
    //     .groupTuple(sort: 'deep')

    // //
    // // MODULE: Create a combined TSV file of all filtered data
    // //
    // CSVTK_CONCAT(
    //     ch_csvtk_concat_input,
    //     "tsv",
    //     "tsv"
    // )
    // ch_concat_compl_data_tsv    = CSVTK_CONCAT.out.csv
    // ch_versions                 = ch_versions.mix(CSVTK_CONCAT.out.versions)


    // /*----------------------------------------------------------------------------
    //     Filter viral sequences/proteins based on classification, composition, and completeness
    // ------------------------------------------------------------------------------*/
    // if (params.run_sequence_filtering) {
    //     // join completeness data and sequences
    //     ch_filter_input = ch_split_fasta_gz
    //         .join(ch_genomad_proteins_faa_gz)
    //         .join(ch_compl_data_tsv)
    //         .multiMap { meta, fasta, proteins, compl_data ->
    //             fasta:      [ meta, fasta ]
    //             proteins:   [ meta, proteins ]
    //             compl_data: [ meta, compl_data ]
    //         }

    //     if (!params.sequences_to_keep) {
    //         ch_seqs_to_keep = []
    //     } else {
    //         ch_seqs_to_keep = file(params.sequences_to_keep, checkIfExists:true)
    //     }

    //     //
    //     // MODULE: Filter sequences and proteins
    //     //
    //     VTDB_FILTERSEQUENCES(
    //         ch_filter_input.fasta,
    //         ch_filter_input.proteins,
    //         ch_filter_input.compl_data,
    //         ch_seqs_to_keep
    //         )
    //     ch_mge_filt_fasta_gz    = VTDB_FILTERSEQUENCES.out.fasta
    //     ch_mge_filt_faa_gz      = VTDB_FILTERSEQUENCES.out.proteins
    //     ch_versions                     = ch_versions.mix(VTDB_FILTERSEQUENCES.out.versions)
    // } else {
    //     ch_mge_filt_fasta_gz    = ch_split_fasta_gz
    //     ch_mge_filt_faa_gz      = ch_genomad_proteins_faa_gz
    // }

    // /*
    // -------------------------------------------------
    //     ASSEMBLY DEREPLICATION
    // -------------------------------------------------
    // */
    // if (params.run_trtrimmer) {
    //     //
    //     // MODULE: Trim tandem repeats from assemblies
    //     //
    //     TRTRIMMER(
    //         ch_mge_filt_fasta_gz,
    //         "${projectDir}/bin/tr-trimmer"
    //     )
    //     ch_versions             = ch_versions.mix(TRTRIMMER.out.versions)
    //     ch_trtrimmer_fasta_gz   = TRTRIMMER.out.fasta
    // } else {
    //     ch_trtrimmer_fasta_gz   = ch_mge_filt_fasta_gz
    // }

    // if (params.run_seqhasher) {
    //     //
    //     // SUBWORKFLOW: Hash assemblies and dereplicate
    //     //
    //     FASTA_SEQHASHER_FASTA(
    //         ch_trtrimmer_fasta_gz,
    //         "${projectDir}/bin/seq-hasher"
    //     )
    //     ch_versions         = ch_versions.mix(FASTA_SEQHASHER_FASTA.out.versions)
    //     ch_derep_fasta_gz   = FASTA_SEQHASHER_FASTA.out.unique_seqs_fasta_gz
    // } else {
    //     ch_derep_fasta_gz   = ch_trtrimmer_fasta_gz
    // }

    // /*
    // -------------------------------------------------
    //     MGE TAXONOMY
    // -------------------------------------------------
    // */
    // // TODO: Add Taxmyphage


    // /*
    // -------------------------------------------------
    //     PHAGE LIFESTYLE
    // -------------------------------------------------
    // */
    // // TODO: Add BACPHLIP
    // // TODO: Add integrase detection



    // /*
    // -------------------------------------------------
    //     MGE HOST
    // -------------------------------------------------
    // */
    // // TODO: Add iPHoP
    // // TODO: Add phist


    // /*
    // -------------------------------------------------
    //     FUNCTIONAL ANNOTATION
    // -------------------------------------------------
    // */
    // // TODO:

    // /*
    // -------------------------------------------------
    //     VIRUS CLUSTERING
    // -------------------------------------------------
    // */
    // if (params.run_vclust) {
    //     FASTA_VCLUST_FASTA(
    //         ch_derep_fasta_gz,
    //         "${projectDir}/bin/vclust"
    //     )
    //     ch_versions         = ch_versions.mix(FASTA_VCLUST_FASTA.out.versions)
    //     ch_vclust_fasta_gz  = FASTA_VCLUST_FASTA.out.cluster_reps_fasta_gz
    // } else {
    //     ch_vclust_fasta_gz  = ch_derep_fasta_gz
    // }

    // /*
    // -------------------------------------------------
    //     READ ALIGNMENT
    // -------------------------------------------------
    // */
    // // TODO: Add Coverm/contig

    // /*
    // -------------------------------------------------
    //     MICRODIVERSITY
    // -------------------------------------------------
    // */
    // // TODO: Add instrain/profile and instrain/compare


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
