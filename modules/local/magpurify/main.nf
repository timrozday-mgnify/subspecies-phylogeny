// MAGpurify — per-contig chimera detection using GC-content and tetranucleotide frequency.
// Runs gc-content and tetra-freq unconditionally (no database required).
// Runs phylo-markers only when a database directory is supplied.
// Outputs are contaminant contig lists for the user to inspect — the genome itself is
// not modified (magpurify clean_bin is intentionally not run).
process MAGPURIFY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magpurify:2.1.2--py_1' :
        'quay.io/biocontainers/magpurify:2.1.2--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)
    path(db)               // MAGpurify database directory for phylo-markers; pass [] to skip

    // MAGpurify writes each module's results into a subdirectory:
    //   ${prefix}_magpurify/gc-content/flagged_contigs
    //   ${prefix}_magpurify/tetra-freq/flagged_contigs
    //   ${prefix}_magpurify/phylo-markers/flagged_contigs (optional)
    output:
    tuple val(meta), path("${prefix}_magpurify/gc-content/flagged_contigs")     , emit: gc_contaminants
    tuple val(meta), path("${prefix}_magpurify/tetra-freq/flagged_contigs")     , emit: tetra_contaminants
    tuple val(meta), path("${prefix}_magpurify/phylo-markers/flagged_contigs")  , optional: true, emit: phylo_contaminants
    tuple val(meta), path("${prefix}_magpurify")                                , emit: results
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    def run_phylo  = (db && db.name != 'input.1') ? "true" : "false"
    """
    # MAGpurify 2.x calls Bio.SeqUtils.GC which was renamed to gc_fraction
    # in BioPython >= 1.80.  Inject a sitecustomize shim so the patch runs
    # automatically at Python startup before any magpurify code is imported.
    mkdir -p /tmp/py_patch
    printf 'try:\\n    import Bio.SeqUtils\\n    if not hasattr(Bio.SeqUtils, \"GC\"):\\n        Bio.SeqUtils.GC = Bio.SeqUtils.gc_fraction\\nexcept ImportError:\\n    pass\\n' \\
        > /tmp/py_patch/sitecustomize.py

    PYTHONPATH="/tmp/py_patch:\${PYTHONPATH:-}" \\
        magpurify gc-content ${fasta} ${prefix}_magpurify ${args}

    PYTHONPATH="/tmp/py_patch:\${PYTHONPATH:-}" \\
        magpurify tetra-freq ${fasta} ${prefix}_magpurify ${args}

    if [ "${run_phylo}" = "true" ]; then
        PYTHONPATH="/tmp/py_patch:\${PYTHONPATH:-}" \\
            magpurify phylo-markers --db ${db} ${fasta} ${prefix}_magpurify ${args}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magpurify: \$(magpurify --version 2>&1 | sed 's/magpurify version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def run_phylo = (db && db.name != 'input.1') ? "true" : "false"
    """
    mkdir -p ${prefix}_magpurify/gc-content
    mkdir -p ${prefix}_magpurify/tetra-freq
    touch ${prefix}_magpurify/gc-content/flagged_contigs
    touch ${prefix}_magpurify/tetra-freq/flagged_contigs
    if [ "${run_phylo}" = "true" ]; then
        mkdir -p ${prefix}_magpurify/phylo-markers
        touch ${prefix}_magpurify/phylo-markers/flagged_contigs
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magpurify: 2.1.2
    END_VERSIONS
    """
}
