// CheckM v1 — lineage-specific marker-gene assessment of genome completeness and contamination.
// Runs lineage_wf on a single genome and exports a tab-delimited QA report.
process CHECKM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm-genome:1.2.3--pyhdfd78af_1' :
        'quay.io/biocontainers/checkm-genome:1.2.3--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)
    path(db)               // CheckM database root directory

    output:
    tuple val(meta), path("${prefix}_checkm_output"), emit: checkm_output
    tuple val(meta), path("${prefix}_checkm_qa.tsv"), emit: checkm_tsv
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # Point CheckM at the mounted database directory
    export CHECKM_DATA_PATH="${db}"

    # CheckM lineage_wf expects a directory of genome files
    mkdir -p input_bins
    cp ${fasta} input_bins/${prefix}.fa

    checkm lineage_wf \\
        --threads ${task.cpus} \\
        -x fa \\
        ${args} \\
        input_bins/ \\
        ${prefix}_checkm_output

    checkm qa \\
        --out_format 2 \\
        --tab_table \\
        -f ${prefix}_checkm_qa.tsv \\
        --threads ${task.cpus} \\
        ${prefix}_checkm_output/lineage.ms \\
        ${prefix}_checkm_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$(checkm --version 2>&1 | sed 's/.*CheckM //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_checkm_output
    touch ${prefix}_checkm_qa.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: 1.2.3
    END_VERSIONS
    """
}
