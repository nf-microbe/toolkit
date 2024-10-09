process CHECKV_ENDTOEND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0':
        'biocontainers/checkv:1.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path ("${prefix}_quality_summary.tsv") , emit: quality_summary
    tuple val(meta), path ("${prefix}_completeness.tsv")    , emit: completeness
    tuple val(meta), path ("${prefix}_contamination.tsv")   , emit: contamination
    tuple val(meta), path ("${prefix}_complete_genomes.tsv"), emit: complete_genomes
    tuple val(meta), path ("${prefix}_aai.tsv")             , emit: aai
    // tuple val(meta), path ("${prefix}/proviruses.fna")      , emit: proviruses
    // tuple val(meta), path ("${prefix}/viruses.fna")         , emit: viruses
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    checkv \\
        end_to_end \\
        ${args} \\
        -t ${task.cpus} \\
        -d ${db} \\
        ${fasta} \\
        ${prefix}

    # rename files by adding prefix
    mv ${prefix}/quality_summary.tsv ${prefix}_quality_summary.tsv
    mv ${prefix}/completeness.tsv ${prefix}_completeness.tsv
    mv ${prefix}/contamination.tsv ${prefix}_contamination.tsv
    mv ${prefix}/complete_genomes.tsv ${prefix}_complete_genomes.tsv
    mv ${prefix}/tmp/aai.tsv ${prefix}_aai.tsv

    rm ${prefix}/proviruses.fna
    rm ${prefix}/viruses.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}_quality_summary.tsv
    touch ${prefix}_completeness.tsv
    touch ${prefix}_contamination.tsv
    touch ${prefix}_complete_genomes.tsv
    touch ${prefix}_aai.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """
}
