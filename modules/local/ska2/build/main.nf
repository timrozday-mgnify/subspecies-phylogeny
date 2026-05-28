process SKA2_BUILD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska2:0.5.1--h4349ce8_0' :
        'quay.io/biocontainers/ska2:0.5.1--h4349ce8_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.skf"),  emit: skf
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ska build \\
        $args \\
        -o ${prefix} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.skf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """
}
