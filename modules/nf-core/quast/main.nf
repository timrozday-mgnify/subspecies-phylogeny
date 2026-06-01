process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321heaaa4ec_4' :
        'quay.io/biocontainers/quast:5.2.0--py39pl5321heaaa4ec_4' }"

    input:
    tuple val(meta), path(consensus)

    output:
    tuple val(meta), path("${prefix}")       , emit: results
    tuple val(meta), path("${prefix}.tsv.gz"), emit: tsv
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    metaquast.py \\
        --output-dir ${prefix} \\
        --threads ${task.cpus} \\
        ${args} \\
        ${consensus.join(' ')}

    gzip -c ${prefix}/report.tsv > ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//;s/ .*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/report.tsv ${prefix}/report.html ${prefix}/quast.log
    touch ${prefix}/transposed_report.txt ${prefix}/icarus.html
    mkdir -p ${prefix}/basic_stats
    touch ${prefix}/basic_stats/cumulative_plot.pdf
    touch ${prefix}/basic_stats/Nx_plot.pdf
    touch ${prefix}/basic_stats/GC_content_plot.pdf

    gzip -c ${prefix}/report.tsv > ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: 5.2.0
    END_VERSIONS
    """
}
