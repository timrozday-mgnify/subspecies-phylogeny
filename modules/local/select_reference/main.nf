// Select the most central (medoid) genome from a FastANI all-vs-all matrix.
// The medoid is the genome with the highest mean ANI to all other samples;
// it minimises gaps when used as the ska map reference and reduces misalignment
// errors compared with a more distant reference.
process SELECT_REFERENCE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gubbins:3.4.3--py310h5140242_0' :
        'quay.io/biocontainers/gubbins:3.4.3--py310h5140242_0' }"

    input:
    path(ani_tsv)
    path(fastas, stageAs: 'fastas/*')

    output:
    path("reference.fa"), emit: reference
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 - <<'PYEOF'
    import os, shutil
    from collections import defaultdict

    ani_sums   = defaultdict(float)
    ani_counts = defaultdict(int)

    with open("${ani_tsv}") as fh:
        for line in fh:
            parts = line.strip().split("\\t")
            if len(parts) < 3:
                continue
            q = os.path.basename(parts[0])
            r = os.path.basename(parts[1])
            if q == r:
                continue
            ani_sums[q]   += float(parts[2])
            ani_counts[q] += 1

    if not ani_sums:
        # Single-sample or empty ANI file: pick the first available FASTA
        fastas = sorted(os.listdir("fastas"))
        shutil.copy(os.path.join("fastas", fastas[0]), "reference.fa")
    else:
        best = max(ani_sums, key=lambda q: ani_sums[q] / ani_counts[q])
        shutil.copy(os.path.join("fastas", best), "reference.fa")
    PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch reference.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
