process SKA2_WEED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska2:0.5.1--h4349ce8_0' :
        'quay.io/biocontainers/ska2:0.5.1--h4349ce8_0' }"

    input:
    tuple val(meta), path(merged_skf)

    output:
    tuple val(meta), path("${prefix}.skf"), path("${prefix}.nkmers.txt"), emit: skf
    path "versions.yml",                                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args overrides the dynamic default; otherwise fall back to the per-task
    // min_freq carried on the meta map (e.g. the computed batch/branch threshold).
    def args = task.ext.args ?: (meta.min_freq != null ? "--min-freq ${meta.min_freq}" : '')
    prefix   = task.ext.prefix ?: "weeded_${meta.id}"
    """
    ska weed \\
        $args \\
        -o ${prefix}.skf \\
        ${merged_skf}

    # A min-freq close to 1 can weed out every k-mer, leaving an unusable SKF;
    # surface the count so the workflow can filter the branch before ska map.
    ska nk ${prefix}.skf | sed -n 's/^k-mers=//p' > ${prefix}.nkmers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "weeded_${meta.id}"
    """
    touch ${prefix}.skf
    echo 100 > ${prefix}.nkmers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """
}
