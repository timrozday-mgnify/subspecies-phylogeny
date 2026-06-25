// Subset a merged/weeded .skf by retaining roughly one split-kmer every N rows.
// Uses the custom ska-minimizer-split image (no bioconda package exists), so the
// `conda` profile is not supported for this module; use docker or singularity.
process SKA2_SUBSET {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.3-sif' :
        'ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.3' }"

    input:
    tuple val(meta), path(skf), path(nkmers_file)
    val(target_snps)

    output:
    tuple val(meta), path("*.subset.skf"), path("*.subset.nkmers.txt"), path("*.subset.sparsity.txt"), emit: skf
    path "versions.yml",                                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "subset_${meta.id}"
    """
    current_snps=\$(cat ${nkmers_file})
    if ! printf '%s\\n' "\$current_snps" | grep -Eq '^[0-9]+\$'; then
        echo "Could not read numeric k-mer count from ${nkmers_file}" >&2
        exit 1
    fi

    target_snps=${target_snps}
    if [ "\$target_snps" -le 0 ]; then
        echo "target_snps must be a positive integer; got \$target_snps" >&2
        exit 1
    fi

    sparsity=\$(awk -v current="\$current_snps" -v target="\$target_snps" 'BEGIN {
        sparsity = target / current
        if (sparsity > 1) {
            sparsity = 1
        }
        printf "%.12g\\n", sparsity
    }')

    ska-shard subset \\
        $args \\
        --sparsity "\$sparsity" \\
        -o ${prefix}.subset.skf \\
        ${skf}

    printf '%s\\n' "\$current_snps" > ${prefix}.subset.nkmers.txt
    printf '%s\\n' "\$sparsity" > ${prefix}.subset.sparsity.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska-minimizer-split: \$(ska-shard --version | sed 's/ska-shard //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "subset_${meta.id}"
    """
    touch ${prefix}.subset.skf
    echo 100 > ${prefix}.subset.nkmers.txt
    echo 1 > ${prefix}.subset.sparsity.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska-minimizer-split: \$(ska-shard --version | sed 's/ska-shard //')
    END_VERSIONS
    """
}
