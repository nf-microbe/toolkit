// Import modules
include { CAT_FASTQ                     } from '../../../modules/nf-core/cat/fastq/main'
include { FASTP                         } from '../../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_RAW          } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_PREPROCESSED } from '../../../modules/nf-core/fastqc/main'

// Import subworkflows
include { getWorkDirs; rmEmptyFastQs    } from '../utils_nfmicrobe_functions/main'
include { FASTQ_BOWTIE2_FASTQ           } from '../fastq_bowtie2_fastq/main'

// Run workflow
workflow FASTQ_READPREPROCESSING_FASTQ {

    take:
    raw_fastq_gz            // channel: [ [ meta.id, meta.single_end, meta.run ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
    run_fastp               // val: false
    perform_run_merging     // val: false
    run_bowtie2_host_removal// val: false
    host_fasta_gz           // val: /path/to/fasta.fasta.gz (OPTIONAL)
    host_bt2_index          // val: /path/to/host_index.bt2 (OPTIONAL)

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW (
        raw_fastq_gz
    )
    ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{ zip -> zip[1] })
    ch_versions         = ch_versions.mix(FASTQC_RAW.out.versions.first())

    if (run_fastp) {
        //
        // MODULE: Run fastp on raw reads
        //
        FASTP(
            raw_fastq_gz,
            [],
            false,
            false,
            false
        )
        ch_fastp_prefilt_fastq_gz   = FASTP.out.reads
        ch_multiqc_files            = ch_multiqc_files.mix(FASTP.out.json.collect{ json -> json[1] })
        ch_versions                 = ch_versions.mix(FASTP.out.versions)

        // REMOVE EMPTY FASTQ FILES
        ch_fastp_fastq_gz = rmEmptyFastQs(ch_fastp_prefilt_fastq_gz, false)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_fastp_workdirs = getWorkDirs(
            raw_fastq_gz,
            ch_fastp_fastq_gz,
            false
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
            ch_pre_fastp_workdirs.map { meta, dir -> [ meta, dir, 'RAW' ] }
        )
    } else {
        ch_fastp_fastq_gz = raw_fastq_gz
    }

    if (perform_run_merging) {
        // prepare reads for concatenating within runs
        ch_reads_forcat = ch_fastp_fastq_gz
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
        CAT_FASTQ(
            ch_reads_forcat.cat.map { meta, reads -> [ meta, reads.flatten() ] }
        )
        ch_cat_reads_fastq_gz   = CAT_FASTQ.out.reads
        ch_versions             = ch_versions.mix(CAT_FASTQ.out.versions)

        // Ensure we don't have nests of nests so that structure is in form expected for assembly
        ch_reads_forcat_skipped = ch_reads_forcat.skip_cat
            .map { meta, reads ->
                def new_reads = meta.single_end ? reads[0] : reads.flatten()
                [ meta, new_reads ]
            }

        // Combine single run and multi-run-merged data
        ch_runmerged_fastq_gz = ch_cat_reads_fastq_gz.mix(ch_reads_forcat_skipped)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_merge_workdirs = getWorkDirs(
            ch_fastp_fastq_gz.map { meta, reads -> [ meta - meta.subMap('run'), reads[0] ] },
            ch_cat_reads_fastq_gz,
            false
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
            ch_pre_merge_workdirs.map { meta, dir -> [ meta, dir, 'FASTP' ] }
        )
    } else {
        ch_runmerged_fastq_gz = ch_fastp_fastq_gz
    }

    if (run_bowtie2_host_removal) {
        //
        // SUBWORKFLOW: Remove host reads using Bowtie2
        //
        FASTQ_BOWTIE2_FASTQ(
            ch_runmerged_fastq_gz,
            host_fasta_gz,
            host_bt2_index
        )
        ch_bt2_prefilt_fastq_gz = FASTQ_BOWTIE2_FASTQ.out.fastq_gz
        ch_multiqc_files        = ch_multiqc_files.mix(FASTQ_BOWTIE2_FASTQ.out.bt2_log.collect{ bt2_log -> bt2_log[1] })
        ch_versions             = ch_versions.mix(FASTQ_BOWTIE2_FASTQ.out.versions)

        // REMOVE EMPTY FASTQ FILES FROM CHANNEL
        ch_bt2_fastq_gz = rmEmptyFastQs(ch_bt2_prefilt_fastq_gz, false)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_bt2_workdirs = getWorkDirs(
            ch_runmerged_fastq_gz.map { meta, reads -> [ meta, reads[0] ] },
            ch_bt2_fastq_gz,
            false
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
            ch_pre_bt2_workdirs.map { meta, dir -> [ meta, dir, 'CAT_RUNMERGE' ] }
        )
    } else {
        ch_bt2_fastq_gz  = ch_runmerged_fastq_gz
    }

    if (perform_run_merging || run_fastp || run_bowtie2_host_removal) {
        //
        // MODULE: Run FastQC on preprocessed reads
        //
        FASTQC_PREPROCESSED (
            ch_bt2_fastq_gz
        )
        ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_PREPROCESSED.out.zip.collect{ zip -> zip[1]})
        ch_versions         = ch_versions.mix(FASTQC_PREPROCESSED.out.versions.first())
    }

    emit:
    preprocessed_fastq_gz   = ch_bt2_fastq_gz       // channel: [ [ meta.id, meta.single_end, meta.run ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    multiqc_files           = ch_multiqc_files      // channel: /path/to/multiqc_files
    versions                = ch_versions           // channel: [ 'software_versions.yml' ]
    workdirs_to_clean       = ch_workdirs_to_clean  // channel: [ [ meta.id ], /path/to/workdir, 'RAW' | 'FASTP' | 'CAT_RUNMERGE' ]
}
