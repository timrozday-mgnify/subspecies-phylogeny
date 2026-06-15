// Concatenate per-bin merged .skf files (ska-shard concat) into one merged.skf.
// The bins partition the split-kmer space, so concatenating their rows is
// equivalent to a single full merge. Inputs are all named merged.skf (one per
// bin), so they are staged into numbered subdirs to avoid name collisions —
// mirrors SKA2_MERGE.
//
// Uses the custom ska-minimizer-split image (no bioconda package exists), so the
// `conda` profile is not supported for this module; use docker or singularity.
process SKA2_SHARD_CONCAT {
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.0-sif' :
        'ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.0' }"

    input:
    path(shard_skfs, stageAs: '?/*')

    output:
    path "merged.skf",   emit: skf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ska-shard concat \\
        $args \\
        -o merged.skf \\
        ${shard_skfs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska-minimizer-split: \$(ska-shard --version | sed 's/ska-shard //')
    END_VERSIONS
    """

    stub:
    """
    touch merged.skf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska-minimizer-split: \$(ska-shard --version | sed 's/ska-shard //')
    END_VERSIONS
    """
}
