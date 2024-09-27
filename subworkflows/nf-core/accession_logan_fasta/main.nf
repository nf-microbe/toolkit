// Import modules
include { LOGAN_CONTIGAWSCLIMULTIPLIER  } from '../../../modules/nf-core/logan/contigawsclimultiplier/main'
include { LOGAN_UNITIGAWSCLIMULTIPLIER  } from '../../../modules/nf-core/logan/unitigawsclimultiplier/main'

workflow ACCESSION_LOGAN_FASTA {

    take:
    accessions              // channel: [ val(meta), [ acc ] ]
    download_contigs        // boolean
    download_unitigs        // boolean


    main:

    ch_versions = Channel.empty()

    ch_logan_mult_fasta_gz  = Channel.empty()

    if (download_contigs) {
        LOGAN_CONTIGAWSCLIMULTIPLIER(
            accessions.map { meta, accession -> [ [ id: "logan_contig_" + meta.id, group: "logan_contig_" + meta.id ], accession] }
        )
        ch_logan_mult_fasta_gz  = ch_logan_mult_fasta_gz.mix(LOGAN_CONTIGAWSCLIMULTIPLIER.out.fasta)
        ch_versions             = ch_versions.mix(LOGAN_CONTIGAWSCLIMULTIPLIER.out.versions.first())
    }

    if (download_unitigs) {
        LOGAN_UNITIGAWSCLIMULTIPLIER(
            accessions.map { meta, accession -> [ [ id: "logan_unitig_" + meta.id, group: "logan_unitig_" + meta.id ], accession] }
        )
        ch_logan_mult_fasta_gz  = ch_logan_mult_fasta_gz.mix(LOGAN_UNITIGAWSCLIMULTIPLIER.out.fasta)
        ch_versions             = ch_versions.mix(LOGAN_UNITIGAWSCLIMULTIPLIER.out.versions.first())
    }

    emit:
    logan_mult_fasta_gz = ch_logan_mult_fasta_gz    // channel: [ [ meta.id, meta.group ], fasta_gz ]
    versions            = ch_versions               // channel: [ versions.yml ]
}

