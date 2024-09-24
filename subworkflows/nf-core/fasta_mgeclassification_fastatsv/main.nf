// Import modules
include { PYHMMER_VIRUSCLASSIFY     } from '../../../modules/nf-core/pyhmmer/virusclassify'
include { PYHMMER_PLASMIDCLASSIFY   } from '../../../modules/nf-core/pyhmmer/plasmidclassify'
include { PYHMMER_BUSCOCLASSIFY     } from '../../../modules/nf-core/pyhmmer/buscoclassify'

// Import subworkflows
include { FASTA_GENOMAD_FASTATSV        } from '../fasta_genomad_fastatsv'
include { getWorkDirs; rmEmptyFastAs    } from '../utils_nfmicrobe_functions'

workflow FASTA_MGECLASSIFICATION_FASTATSV {

    take:
    fasta_gz                    // [ [ meta ], fasta.gz ]
    faa_gz                      // [ [ meta ], faa.gz ]
    run_genomad                 // boolean: false
    genomad_db                  // path: path/to/genomad_db
    run_pyhmmer_virus           // boolean: false
    run_pyhmmer_plasmid         // boolean: false
    run_pyhmmer_busco           // boolean: false

    main:

    ch_versions             = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()

    if (run_genomad) {
        //
        // SUBWORKFLOW: Classify MGEs with geNomad
        //
        FASTA_GENOMAD_FASTATSV(
            fasta_gz,
            genomad_db
        )
        ch_versions = ch_versions.mix(FASTA_GENOMAD_FASTATSV.out.versions.first())
        ch_genomad_genes_tsv    = FASTA_GENOMAD_FASTATSV.out.genes_tsv
        ch_genomad_features_tsv = FASTA_GENOMAD_FASTATSV.out.features_tsv
        ch_genomad_scores_tsv   = FASTA_GENOMAD_FASTATSV.out.scores_tsv
        ch_genomad_taxonomy_tsv = FASTA_GENOMAD_FASTATSV.out.taxonomy_tsv
    } else {
        ch_genomad_genes_tsv    = []
        ch_genomad_features_tsv = []
        ch_genomad_scores_tsv   = []
        ch_genomad_taxonomy_tsv = []
    }

    if (run_pyhmmer_virus) {
        //
        // MODULE: Classify MGEs with pyHMMER
        //
        PYHMMER_VIRUSCLASSIFY(
            faa_gz,
            "${projectDir}/assets/hmms/virus_hmms/DJR_MCP_virus_hallmarks.hmm",
            "${projectDir}/assets/hmms/virus_hmms/inovirus_MCP_virus_hallmarks.hmm",
            "${projectDir}/assets/hmms/virus_hmms/pleolipoviridae_virus_hallmarks.hmm"
        )
        ch_versions                 = ch_versions.mix(PYHMMER_VIRUSCLASSIFY.out.versions.first())
        ch_pyhmmer_virus_hmms_tsv   = PYHMMER_VIRUSCLASSIFY.out.markers
    } else {
        ch_pyhmmer_virus_hmms_tsv = []
    }

    if (run_pyhmmer_plasmid) {
        //
        // MODULE: Classify MGEs with pyHMMER
        //
        PYHMMER_PLASMIDCLASSIFY(
            faa_gz,
            "${projectDir}/assets/hmms/plasmid_hmms/CONJscan_plasmid_hallmarks.hmm",
            "${projectDir}/assets/hmms/plasmid_hmms/Pfam_NCBIfam_plasmid_hallmarks.hmm",
        )
        ch_versions                 = ch_versions.mix(PYHMMER_PLASMIDCLASSIFY.out.versions.first())
        ch_pyhmmer_plasmid_hmms_tsv = PYHMMER_PLASMIDCLASSIFY.out.markers
    } else {
        ch_pyhmmer_plasmid_hmms_tsv = []
    }

    if (run_pyhmmer_busco) {
        //
        // MODULE: Classify MGEs with pyHMMER
        //
        PYHMMER_BUSCOCLASSIFY(
            faa_gz,
            "${projectDir}/assets/hmms/busco_hmms/archaea_buscos.hmm",
            "${projectDir}/assets/hmms/busco_hmms/score_cutoffs/archaea_odb10.cutoffs",
            "${projectDir}/assets/hmms/busco_hmms/bacteria_buscos.hmm",
            "${projectDir}/assets/hmms/busco_hmms/score_cutoffs/bacteria_odb10.cutoffs"
        )
        ch_versions                 = ch_versions.mix(PYHMMER_BUSCOCLASSIFY.out.versions.first())
        ch_pyhmmer_busco_hmms_tsv   = PYHMMER_BUSCOCLASSIFY.out.markers
    } else {
        ch_pyhmmer_busco_hmms_tsv = []
    }

    emit:
    genomad_genes_tsv           = ch_genomad_genes_tsv
    genomad_features_tsv        = ch_genomad_features_tsv
    genomad_scores_tsv          = ch_genomad_scores_tsv
    genomad_taxonomy_tsv        = ch_genomad_taxonomy_tsv
    pyhmmer_virus_hmms_tsv      = ch_pyhmmer_virus_hmms_tsv
    pyhmmer_plasmid_hmms_tsv    = ch_pyhmmer_plasmid_hmms_tsv
    pyhmmer_busco_hmms_tsv      = ch_pyhmmer_busco_hmms_tsv
    versions                    = ch_versions
}

