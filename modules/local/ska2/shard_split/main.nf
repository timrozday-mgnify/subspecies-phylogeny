// Split one sample's merged-or-single .skf into n_shards bins by the minimizer
// value of each split k-mer (ska-shard split). Because split-kmer keys are
// canonical, the same k-mer lands in the same bin in every sample, so bins can
// be merged independently and concatenated — see SKA2_SHARDED_MERGE.
//
// Uses the custom ska-minimizer-split image (no bioconda package exists), so the
// `conda` profile is not supported for this module; use docker or singularity.
process SKA2_SHARD_SPLIT {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.1-sif' :
        'ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.1' }"

    input:
    tuple val(meta), path(skf)
    val(n_shards)
    val(minimizer_len)

    output:
    tuple val(meta), path("${prefix}.*.skf"), emit: bins
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    ska-shard split \\
        $args \\
        -n ${n_shards} \\
        -l ${minimizer_len} \\
        -o ${prefix} \\
        ${skf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska-minimizer-split: \$(ska-shard --version | sed 's/ska-shard //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    for i in \$(seq 0 \$(( ${n_shards} - 1 ))); do
        touch ${prefix}.\$i.skf
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska-minimizer-split: \$(ska-shard --version | sed 's/ska-shard //')
    END_VERSIONS
    """
}
