process GENOMAD_ENDTOEND {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.8.0--pyhdfd78af_1':
        'biocontainers/genomad:1.8.0--pyhdfd78af_1' }"

    input:
    tuple val(meta) , path(fasta)
    path  genomad_db

    output:
    tuple val(meta), path("${prefix}_plasmid_genes.tsv")                , emit: plasmid_genes
    tuple val(meta), path("${prefix}_plasmid_summary.tsv")              , emit: plasmid_summary
    tuple val(meta), path("${prefix}_virus_genes.tsv")                  , emit: virus_genes
    tuple val(meta), path("${prefix}_virus_summary.tsv")                , emit: virus_summary
    tuple val(meta), path("${prefix}_virus.fna.gz")                     , emit: virus_fasta
    tuple val(meta), path("${prefix}_provirus.tsv")                     , emit: provirus
    tuple val(meta), path("${prefix}_genes.tsv")                        , emit: genes
    tuple val(meta), path("${prefix}_features.tsv")                     , emit: features
    tuple val(meta), path("${prefix}_aggregated_classification.tsv")    , emit: scores
    tuple val(meta), path("${prefix}_taxonomy.tsv")                     , emit: taxonomy
    tuple val(meta), path("${prefix}_proteins.faa.gz")                  , emit: faa
    path "versions.yml"                                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def filename = "${fasta}"[0..<"${fasta}".lastIndexOf('.')]
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    genomad \\
        end-to-end \\
        ${fasta} \\
        ./ \\
        ${genomad_db} \\
        --threads ${task.cpus} \\
        ${args}

    # remove provirus outputs to avoid wildcard issues
    if [ -f *_aggregated_classification/*_provirus_aggregated_classification.tsv ]; then
        rm *_aggregated_classification/*_provirus_aggregated_classification.tsv
        rm *_marker_classification/*_provirus_features.tsv
    fi

    # save virus outputs
    gzip -c *_summary/*_virus.fna > ${prefix}_virus.fna.gz
    mv *_summary/*_virus_summary.tsv ${prefix}_virus_summary.tsv
    mv *_summary/*_virus_genes.tsv ${prefix}_virus_genes.tsv
    mv *_find_proviruses/*_provirus.tsv ${prefix}_provirus.tsv

    # save plasmid outputs
    mv *_summary/*_plasmid_summary.tsv ${prefix}_plasmid_summary.tsv
    mv *_summary/*_plasmid_genes.tsv ${prefix}_plasmid_genes.tsv

    # save other important outputs
    mv *_aggregated_classification/*_aggregated_classification.tsv ${prefix}_aggregated_classification.tsv
    mv *_annotate/*_taxonomy.tsv ${prefix}_taxonomy.tsv
    mv *_annotate/*_genes.tsv ${prefix}_genes.tsv
    mv *_marker_classification/*_features.tsv ${prefix}_features.tsv
    gzip -c *_annotate/*_proteins.faa > ${prefix}_proteins.faa.gz

    # clean output directories
    rm -rf ./${filename}_*/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def filename = "${fasta}"[0..<"${fasta}".lastIndexOf('.')]
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_virus.fna.gz
    touch ${prefix}_virus_summary.tsv
    touch ${prefix}_virus_genes.tsv
    touch ${prefix}_provirus.tsv

    touch ${prefix}_plasmid_summary.tsv
    touch ${prefix}_plasmid_genes.tsv

    echo "" | gzip > ${prefix}_proteins.faa.gz
    touch ${prefix}_aggregated_classification.tsv
    touch ${prefix}_taxonomy.tsv
    touch ${prefix}_genes.tsv
    touch ${prefix}_features.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}
