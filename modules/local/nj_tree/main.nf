process NJ_TREE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/r-ape:5.8--48d6804841ebe369' :
        'community.wave.seqera.io/library/r-ape:5.8--48d6804841ebe369' }"

    input:
    tuple val(meta), path(distances)

    output:
    tuple val(meta), path("*.nwk"), emit: tree
    tuple val(meta), path("*.pdf"), emit: plot
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    nj_tree.R \\
        --input ${distances} \\
        --prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | grep "^R version" | sed 's/R version \\([0-9.]*\\).*/\\1/')
        r-ape: \$(Rscript -e "cat(as.character(packageVersion('ape')), '\\n')" 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nwk ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | grep "^R version" | sed 's/R version \\([0-9.]*\\).*/\\1/')
        r-ape: \$(Rscript -e "cat(as.character(packageVersion('ape')), '\\n')" 2>&1)
    END_VERSIONS
    """
}
