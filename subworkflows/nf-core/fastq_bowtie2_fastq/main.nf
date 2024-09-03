// Import modules
include { BOWTIE2_BUILD } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN } from '../../../modules/nf-core/bowtie2/align/main'

// Run subworkflow
workflow FASTQ_BOWTIE2_FASTQ {

    take:
    fastq_gz        // channel: [ [ meta ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
    fasta_gz        // val: /path/to/fasta.fasta.gz (OPTIONAL)
    bt2_index       // val: /path/to/host_index.bt2 (OPTIONAL)

    main:
    ch_versions = Channel.empty()

    // prepare fasta and bowtie2 index
    if (fasta_gz && bt2_index) {
        ch_bowtie2_fasta    = Channel.value(
            [ [ id:'bowtie2_fasta' ], file(fasta_gz, checkIfExists: true) ]
        )
        ch_bowtie2_index    = Channel.value(
            [ [ id:'bowtie2_index' ], file(bt2_index, checkIfExists: true) ]
        )
    } else {
        ch_bowtie2_fasta    = Channel.value(
            [ [ id:'bowtie2_fasta' ], file(fasta_gz, checkIfExists: true) ]
        )
        ch_bowtie2_index   = null
    }

    //
    // MODULE: Create bowtie2 database from host FASTA
    //
    if (!ch_bowtie2_index) {
        BOWTIE2_BUILD(
            ch_bowtie2_fasta
        )
        ch_bowtie2_index    = BOWTIE2_BUILD.out.index
        ch_versions         = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    //
    // MODULE: Map, generate BAM with all reads and unmapped reads in FASTQ for downstream
    //
    BOWTIE2_ALIGN(
        fastq_gz,
        ch_bowtie2_index,
        ch_bowtie2_fasta,
        true,
        true
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    emit:
    fastq_gz    = BOWTIE2_ALIGN.out.fastq   // channel: [ val(meta), [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    bt2_log     = BOWTIE2_ALIGN.out.log     // channel: [ val(meta), bt2.log ]
    versions    = ch_versions               // channel: [ software_versions.yml ]
}

