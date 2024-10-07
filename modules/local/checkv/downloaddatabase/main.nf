process CHECKV_DOWNLOADDATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0':
        'biocontainers/checkv:1.0.1--pyhdfd78af_0' }"

    output:
    path "checkv_db/*"  , emit: checkv_db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    checkv download_database \\
        ${args} \\
        ./checkv_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    mkdir checkv_db
    touch checkv_db/README.txt
    mkdir checkv_db/genome_db
    touch checkv_db/genome_db/changelog.tsv
    touch checkv_db/genome_db/checkv_error.tsv
    touch checkv_db/genome_db/checkv_info.tsv
    touch checkv_db/genome_db/checkv_reps.faa
    touch checkv_db/genome_db/checkv_reps.fna
    touch checkv_db/genome_db/checkv_reps.tsv
    mkdir checkv_db/hmm_db
    touch checkv_db/hmm_db/checkv_hmms.tsv
    touch checkv_db/hmm_db/genome_lengths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """

}
