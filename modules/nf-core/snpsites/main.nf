process SNPSITES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snp-sites:2.5.1--hed695b0_0' :
        'quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.fas"),        emit: fasta
    tuple val(meta), path("*.sites.txt"),  emit: constant_sites
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snp-sites \\
        $alignment \\
        $args \\
        > ${prefix}.fas

    snp-sites -C $alignment > ${prefix}.sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsites: \$(snp-sites -V 2>&1 | sed 's/snp-sites //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fas
    echo '0,0,0,0' > ${prefix}.sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsites: \$(snp-sites -V 2>&1 | sed 's/snp-sites //')
    END_VERSIONS
    """
}
