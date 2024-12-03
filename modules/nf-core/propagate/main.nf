process PROPAGATE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1' :
        'biocontainers/vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta) , path(fastq)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(genomad)

    output:
    tuple val(meta), path("${prefix}.propagate.tsv")    , emit: results
    tuple val(meta), path("${prefix}.propagate.log")    , emit: log     , optional: true
    path  "versions.yml"                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    # convert genomad to coords
    echo -e "scaffold\tfragment\tstart\tstop" > ${prefix}.coords.tsv
    cat ${genomad} | \\
    grep "|provirus" | \\
    awk -F'\t' '{ print \$1 "\t" \$1 "\t" \$4 }' | sed 's/|provirus[^\t]*//1; s/-/\t/' \\
    >> ${prefix}.coords.tsv

    # remove description (after space) from sequence headers
    sed '/^>/ s/ .*//' ${fasta_name} > ${prefix}.renamed.fasta

    Propagate \\
        -f ${prefix}.renamed.fasta \\
        -v ${prefix}.coords.tsv \\
        -r ${fastq} \\
        -o ${prefix}.propagate \\
        -t ${task.cpus} \\
        ${args}

    mv ${prefix}.propagate/${prefix}.propagate.tsv ./${prefix}.propagate.tsv
    mv ${prefix}.propagate/${prefix}.propagate.log ./${prefix}.propagate.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        propagate: \$(Propagate --version 2>&1 | sed 's/^.*PropagAtE v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.propagate.tsv
    touch ${prefix}.propagate.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        propagate: \$(Propagate --version 2>&1 | sed 's/^.*PropagAtE v//; s/ .*\$//')
    END_VERSIONS
    """
}
