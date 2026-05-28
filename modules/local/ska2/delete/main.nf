process SKA2_DELETE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska2:0.5.1--h4349ce8_0' :
        'quay.io/biocontainers/ska2:0.5.1--h4349ce8_0' }"

    input:
    path(merged_skf)
    path(samples_file)

    output:
    path("*.skf"),        emit: skf
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: 'merged_deleted'
    // ska delete -f expects tab-separated name<TAB>path lines; convert from newline-delimited names
    """
    awk '{print \$0"\\tdummy"}' ${samples_file} > _delete_list.tsv

    ska delete \\
        $args \\
        --skf-file ${merged_skf} \\
        -f _delete_list.tsv \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: 'merged_deleted'
    """
    touch ${prefix}.skf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """
}
