process BUSCO_BUSCO {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.4.7--pyhdfd78af_0' :
        'quay.io/biocontainers/busco:5.4.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    // Lineage name (e.g. 'bacteria_odb10') or auto mode (e.g. 'auto_prok').
    // Passed to BUSCO via --lineage_dataset or --auto-lineage-prok.
    val lineage
    // Directory containing pre-downloaded BUSCO lineage datasets.
    // Pass [] to let BUSCO download them at runtime.
    path busco_lineages_path

    output:
    tuple val(meta), path("${prefix}-busco.batch_summary.txt"), emit: batch_summary
    tuple val(meta), path("short_summary.*.txt")              , emit: short_summaries, optional: true
    tuple val(meta), path("${prefix}-busco")                  , emit: busco_dir
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    prefix           = task.ext.prefix ?: "${meta.id}"
    def lineage_cmd  = lineage in ['auto', 'auto_prok', 'auto_euk']
        ? lineage.replaceFirst('auto', '--auto-lineage').replaceAll('_', '-')
        : "--lineage_dataset ${lineage}"
    // busco_lineages_path should be the directory containing the 'lineages/' subdirectory.
    // Use .toRealPath() so the bash script receives an absolute path — required because
    // we run BUSCO in a tmpfs directory (cd /tmp) to work around a macOS+Docker
    // multiprocessing.Pool permission error on bind-mounted VirtioFS directories.
    def lineage_dir  = (busco_lineages_path && busco_lineages_path.name != 'input.1')
        ? "--download_path ${busco_lineages_path.toRealPath()} --offline"
        : ''
    """
    # BUSCO 5.4.7 uses multiprocessing.Pool internally.  On macOS+Docker (Colima /
    # VirtioFS), child processes forked by Pool cannot create files on bind-mounted
    # host directories, producing PermissionError on logs/hmmsearch_out.log.
    # Fix: run BUSCO in a container-local tmpfs directory so Pool's child writes
    # to an in-container filesystem; copy the results back afterwards.
    _work="\${PWD}"
    _fasta="\$(realpath ${fasta})"
    _tmp="\$(mktemp -d)"

    cd "\${_tmp}"

    busco \\
        --cpu ${task.cpus} \\
        --in "\${_fasta}" \\
        --out ${prefix}-busco \\
        --mode genome \\
        ${lineage_cmd} \\
        ${lineage_dir} \\
        ${args}

    # Copy output directory and summary files back to the Nextflow work directory
    cp -r ${prefix}-busco "\${_work}/"
    mv ${prefix}-busco/batch_summary.txt "\${_work}/${prefix}-busco.batch_summary.txt" 2>/dev/null || \
        touch "\${_work}/${prefix}-busco.batch_summary.txt"
    find ${prefix}-busco -name 'short_summary.*.txt' -exec cp {} "\${_work}/" \\; 2>/dev/null || true

    cd "\${_work}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version 2>&1 | sed 's/BUSCO //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}-busco
    touch ${prefix}-busco.batch_summary.txt
    touch "short_summary.specific.${lineage}.${prefix}-busco.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: 5.4.7
    END_VERSIONS
    """
}
