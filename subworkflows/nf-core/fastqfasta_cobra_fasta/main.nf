// Import modules
include { COVERM_CONTIG as COBRA_CONTIG } from '../../../modules/nf-core/coverm/contig/main'
include { COBRAMETA                     } from '../../../modules/nf-core/cobrameta/main'

workflow FASTQFASTA_COBRA_FASTA {

    take:
    fastq_gz            // channel: [ [ meta.id, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
    fasta_gz            // channel: [ [ meta.id, meta.single_end, meta.assembler ], path(fasta)]
    query_contigs_tsv   // channel: [ [ meta.id, meta.single_end, meta.assembler ], path(queries) ]
    assembly_logs       // channel: [ [ meta.id, meta.single_end, meta.assembler ], path(logs) ]

    main:

    ch_versions = Channel.empty()

    // join fastq and fasta by meta.id
    ch_coverm_input = fastq_gz.map { meta, fastq -> [ meta.id, fastq ] }
        .combine(fasta_gz.map { meta, fasta -> [ meta.id, meta, fasta ] }, by:0)
        .multiMap { meta_id, fastq, meta_fasta, fasta ->
            fastq: [ meta_fasta, fastq ]
            fasta: [ meta_fasta, fasta ]
        }

    //
    // MODULE: Align reads to their corresponding assembly
    //
    COBRA_CONTIG(
        ch_coverm_input.fastq,
        ch_coverm_input.fasta
    )
    ch_versions = ch_versions.mix(COBRA_CONTIG.out.versions.first())

    // prepare input for cobra
    ch_cobra_input = fasta_gz
        .combine(query_contigs_tsv, by:0)
        .combine(COBRA_CONTIG.out.tsv, by:0)
        .combine(COBRA_CONTIG.out.bam, by:0)
        .combine(assembly_logs, by:0)
        .multiMap { meta, fasta, query, coverage, bam, log ->
            fasta:      [ meta, fasta ]
            coverage:   [ meta, coverage ]
            query:      [ meta, query ]
            bam:        [ meta, bam ]
            log:        [ meta, log ]
        }

    //
    // MODULE: Extend query contigs
    //
    COBRAMETA(
        ch_cobra_input.fasta.map{ meta, fasta -> [ meta + [ extension: 'cobra' ], fasta ] },
        ch_cobra_input.coverage,
        ch_cobra_input.query,
        ch_cobra_input.bam,
        ch_cobra_input.log,
        ch_cobra_input.fasta.map { meta, fasta -> meta.assembler },
    )

    emit:
    fasta_gz    = COBRAMETA.out.fasta   // channel: [ [ meta.id, meta.assembler, meta.extension ], [ fasta.gz ] ]
    versions    = ch_versions           // channel: [ versions.yml ]
}
