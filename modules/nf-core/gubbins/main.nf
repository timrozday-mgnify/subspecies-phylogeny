process GUBBINS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gubbins:3.4.3--py310h5140242_0' :
        'quay.io/biocontainers/gubbins:3.4.3--py310h5140242_0' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.filtered_polymorphic_sites.fasta"), emit: fasta
    tuple val(meta), path("*.recombination_predictions.gff"),    emit: gff
    tuple val(meta), path("*.summary_of_snp_distribution.vcf"),  emit: vcf
    tuple val(meta), path("*.filtered_polymorphic_sites.phylip"), emit: phylip
    tuple val(meta), path("*.final_tree.tre"),                   emit: tree, optional: true
    tuple val(meta), path("*.per_branch_statistics.csv"),        emit: stats, optional: true
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR="\$PWD/.numba_cache"

    # Cap threads at the number of cores visible to the container (e.g. Docker
    # Desktop may expose fewer cores than Nextflow's task.cpus allocation).
    threads=\$(( ${task.cpus} < \$(nproc) ? ${task.cpus} : \$(nproc) ))

    run_gubbins.py \\
        --threads \$threads \\
        --prefix ${prefix} \\
        $args \\
        ${alignment}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1 | sed 's/^run_gubbins.py //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered_polymorphic_sites.fasta
    touch ${prefix}.recombination_predictions.gff
    touch ${prefix}.summary_of_snp_distribution.vcf
    touch ${prefix}.filtered_polymorphic_sites.phylip
    touch ${prefix}.final_tree.tre
    touch ${prefix}.per_branch_statistics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1 | sed 's/^run_gubbins.py //')
    END_VERSIONS
    """
}
