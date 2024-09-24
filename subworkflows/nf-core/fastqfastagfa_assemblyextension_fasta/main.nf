// Import subworkflows
include { FASTQFASTA_COBRA_FASTA        } from '../fastqfasta_cobra_fasta'
include { FASTQGFA_PHABLES_FASTA        } from '../fastqgfa_phables_fasta'
include { getWorkDirs; rmEmptyFastAs    } from '../utils_nfmicrobe_functions/main'

workflow FASTQFASTAGFA_ASSEMBLYEXTENSION_FASTA {

    take:
    fastq_gz            // channel: [ [ meta.id, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
    run_cobra           // boolean: false
    fasta_gz            // channel: [ [ meta.id, meta.single_end, meta.assembler ], path(fasta)]
    assembly_logs       // channel: [ [ meta.id, meta.single_end, meta.assembler ], path(logs) ]
    query_contigs_tsv   // channel: [ [ meta.id, meta.single_end, meta.assembler ], path(queries) ]
    run_phables         // boolean: false
    gfa_gz              // channel: [ [ meta.id ], gfa.gz ]
    phables_config      // channel: [ phables.config ]
    phables_db          // channel: [ phables_db ]

    main:

    ch_versions                     = Channel.empty()
    ch_extended_prefilt_fasta_gz    = Channel.empty()

    if (run_cobra) {
        //
        // SUBWORKFLOW: Extend assemblies with cobra
        //
        FASTQFASTA_COBRA_FASTA(
            fastq_gz,
            fasta_gz,
            query_contigs_tsv,
            assembly_logs
        )
        ch_extended_prefilt_fasta_gz    = ch_extended_prefilt_fasta_gz.mix(FASTQFASTA_COBRA_FASTA.out.fasta_gz)
        ch_versions                     = ch_versions.mix(FASTQFASTA_COBRA_FASTA.out.versions.first())
    }

    if (run_phables) {
        //
        // SUBWORKFLOW: Extend assemblies with phables
        //
        FASTQGFA_PHABLES_FASTA(
            fastq_gz,
            gfa_gz,
            phables_config,
            phables_db
        )
        ch_extended_prefilt_fasta_gz    = ch_extended_prefilt_fasta_gz.mix(FASTQGFA_PHABLES_FASTA.out.fasta_gz)
        ch_versions                     = ch_versions.mix(FASTQGFA_PHABLES_FASTA.out.versions.first())
    }

    // REMOVE EMPTY FASTA FILES FROM CHANNEL
    ch_extended_fasta_gz = rmEmptyFastAs(ch_extended_prefilt_fasta_gz, false)

    emit:
    fasta_gz    = ch_extended_fasta_gz  // channel: [ val(meta), [ fasta.gz ] ]
    versions    = ch_versions           // channel: [ versions.yml ]
}

