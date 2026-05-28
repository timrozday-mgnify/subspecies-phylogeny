process SKA2_LO {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska2:0.5.1--h4349ce8_0' :
        'quay.io/biocontainers/ska2:0.5.1--h4349ce8_0' }"

    input:
    path(merged_skf)
    path(reference)

    output:
    path("${prefix}_snps.fas"),    emit: snps
    path("${prefix}_indels.vcf"),  emit: indels
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix     = task.ext.prefix ?: 'lo_output'
    def ref_arg = reference ? "-r ${reference}" : ''
    """
    ska lo \\
        $args \\
        $ref_arg \\
        --threads $task.cpus \\
        ${merged_skf} \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'lo_output'
    """
    touch ${prefix}_snps.fas
    touch ${prefix}_indels.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """
}
