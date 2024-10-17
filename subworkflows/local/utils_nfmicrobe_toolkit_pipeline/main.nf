//
// Subworkflow with functionality specific to the nf-microbe/toolkit pipeline
//


//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def metas = input[1]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }
    // Check that multiple runs of the same sample are placed in the same group
    def grouping_ok = metas.collect{ meta -> meta.group }.unique().size == 1
    if (!grouping_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be placed into the same group: ${metas[0].id}")
    }
    // Check that multiple runs of the same sample are given different run ids
    def runs_ok   = metas.collect{ meta -> meta.run }.unique().size == metas.collect{ meta -> meta.run }.size
    if (!runs_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be given a different run id: ${metas[0].id}")
    }
}
