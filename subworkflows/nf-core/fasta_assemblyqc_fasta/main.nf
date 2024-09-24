// Import modules
include { SEQKIT_SEQ            } from '../../../modules/nf-core/seqkit/seq'
include { SEQKIT_STATS          } from '../../../modules/nf-core/seqkit/stats'
include { TRFINDER              } from '../../../modules/nf-core/trfinder'
include { TANTAN                } from '../../../modules/nf-core/tantan'
include { PYRODIGALGV           } from '../../../modules/nf-core/pyrodigalgv'
include { SEQUENCESTATS         } from '../../../modules/nf-core/sequencestats'

// Import subworkflows
include { getWorkDirs; rmEmptyFastAs    } from '../utils_nfmicrobe_functions'

workflow FASTA_ASSEMBLYQC_FASTA {

    take:
    fasta_gz                // channel: [ val(meta), [ fasta.gz ] ]
    min_len                 // val: int
    run_seqkit_stats        // boolean: false
    run_trfinder            // boolean: false
    use_trfinder_fasta      // boolean: false
    run_tantan              // boolean: false
    run_sequence_stats      // boolean: false

    main:

    ch_versions             = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()

    if (min_len > 0) {
        //
        // MODULE: Filter assemblies by length
        //
        SEQKIT_SEQ(
            fasta_gz
        )
        ch_seqkit_seq_prefilt_fasta_gz  = SEQKIT_SEQ.out.fastx
        ch_versions                     = ch_versions.mix(SEQKIT_SEQ.out.versions.first())

        // REMOVE EMPTY FASTA FILES FROM CHANNEL
        ch_seqkit_seq_fasta_gz = rmEmptyFastAs(ch_seqkit_seq_prefilt_fasta_gz, false)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_seqkit_workdirs = getWorkDirs(
            fasta_gz,
            ch_seqkit_seq_fasta_gz,
            false
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
            ch_pre_seqkit_workdirs.map { meta, dir -> [ meta, dir, 'SEQKIT_SEQ' ] }
        )
    } else {
        ch_seqkit_seq_fasta_gz = fasta_gz
    }

    if (run_seqkit_stats) {
        //
        // MODULE: Calculate assembly statistics
        //
        SEQKIT_STATS(
            ch_seqkit_seq_fasta_gz
        )
        ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())
    }

    if (run_trfinder) {
        //
        // MODULE: Identify sequences with terminal repeats
        //
        TRFINDER(
            ch_seqkit_seq_fasta_gz
        )
        ch_trfinder_prefilt_fasta_gz    = TRFINDER.out.fasta
        ch_trfinder_tsv                 = TRFINDER.out.stats
        ch_versions                     = ch_versions.mix(TRFINDER.out.versions.first())

        if (use_trfinder_fasta) {
            // REMOVE EMPTY FASTA FILES FROM CHANNEL
            ch_trfinder_fasta_gz = rmEmptyFastAs(ch_trfinder_prefilt_fasta_gz, false)

            // IDENTIFY WORKDIRS TO CLEAN
            ch_pre_trfinder_workdirs = getWorkDirs(
                ch_seqkit_seq_fasta_gz,
                ch_trfinder_fasta_gz,
                false
            )
            ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
                ch_pre_trfinder_workdirs.map { meta, dir -> [ meta, dir, 'TRFINDER' ] }
            )
        } else {
            ch_trfinder_fasta_gz = ch_seqkit_seq_fasta_gz
        }
    } else {
        ch_trfinder_fasta_gz    = ch_seqkit_seq_fasta_gz
        ch_trfinder_tsv         = []
    }

    if (run_tantan) {
        //
        // MODULE: identify low-complexity regions
        //
        TANTAN(
            ch_trfinder_fasta_gz
        )
        ch_tantan_bed   = TANTAN.out.bed
        ch_versions     = ch_versions.mix(TANTAN.out.versions.first())
    } else {
        ch_tantan_bed = []
    }

    //
    // MODULE: Predict open reading frames
    //
    PYRODIGALGV(
        ch_trfinder_fasta_gz
    )
    ch_pyrodigalgv_faa_gz   = PYRODIGALGV.out.faa
    ch_versions             = ch_versions.mix(PYRODIGALGV.out.versions.first())

    if (run_sequence_stats) {
        // join inputs by meta.id for nucleotide stats
        ch_seq_stats_inputs = ch_trfinder_fasta_gz
            .join(ch_pyrodigalgv_faa_gz)
            .multiMap { meta, fasta, faa ->
                fasta:  [ meta, fasta ]
                faa:    [ meta, faa ]
            }
        //
        // MODULE: Calculate nucleotide statistics
        //
        SEQUENCESTATS(
            ch_seq_stats_inputs.fasta,
            ch_seq_stats_inputs.faa
        )
        ch_seq_stats_tsv    = SEQUENCESTATS.out.stats
        ch_versions         = ch_versions.mix(SEQUENCESTATS.out.versions.first())
    } else {
        ch_seq_stats_tsv = []
    }

    emit:
    filtered_fasta_gz   = ch_trfinder_fasta_gz  // channel: [ [ meta.id ], fasta.gz ]
    trfinder_tsv        = ch_trfinder_tsv       // channel: [ [ meta.id ], trfinder.tsv ]
    tantan_bed          = ch_tantan_bed         // channel: [ [ meta.id ], tantan.bed ]
    pyrodigalgv_faa_gz  = ch_pyrodigalgv_faa_gz // channel: [ [ meta.id ], faa.gz ]
    seq_stats_tsv       = ch_seq_stats_tsv      // channel: [ [ meta.id ], stats.tsv ]
    workdirs_to_clean   = ch_workdirs_to_clean  // channel: [ [ meta.id, workdir, module ] ]
    versions            = ch_versions           // path: [ versions.yml ]
}

