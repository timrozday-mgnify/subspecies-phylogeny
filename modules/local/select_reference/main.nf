// Select the most central (medoid) genome from a ska2 pairwise SNP distance
// matrix. The medoid is the genome with the lowest mean distance (the
// `Mismatches (proportion)` column, matching the ska2 NJ tree) to all other
// samples; it minimises gaps when used as the ska map reference and reduces
// misalignment errors compared with a more distant reference.
process SELECT_REFERENCE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gubbins:3.4.3--py310h5140242_0' :
        'quay.io/biocontainers/gubbins:3.4.3--py310h5140242_0' }"

    input:
    path(distance_tsv)
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

    # Map ska sample name (FASTA basename sans extension) -> staged FASTA filename.
    # Only staged FASTAs are candidates (i.e. those that passed the trusted filter
    # upstream); medoid means are still taken over every sample in the matrix.
    candidates = {}
    for fn in os.listdir("fastas"):
        candidates[os.path.splitext(fn)[0]] = fn

    dist_sums   = defaultdict(float)
    dist_counts = defaultdict(int)

    with open("${distance_tsv}") as fh:
        header = fh.readline().rstrip("\\n").split("\\t")
        try:
            di = header.index("Mismatches (proportion)")
        except ValueError:
            di = 3  # ska distance default column order
        for line in fh:
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) <= di:
                continue
            a, b = parts[0], parts[1]
            try:
                d = float(parts[di])
            except ValueError:
                continue
            # Accumulate distance from each candidate to every other sample.
            if a in candidates:
                dist_sums[a] += d; dist_counts[a] += 1
            if b in candidates:
                dist_sums[b] += d; dist_counts[b] += 1

    if not dist_sums:
        # Single-sample or empty distance file: pick the first available FASTA.
        first = sorted(candidates.values())[0]
        shutil.copy(os.path.join("fastas", first), "reference.fa")
    else:
        best = min(dist_sums, key=lambda s: dist_sums[s] / dist_counts[s])
        shutil.copy(os.path.join("fastas", candidates[best]), "reference.fa")
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
