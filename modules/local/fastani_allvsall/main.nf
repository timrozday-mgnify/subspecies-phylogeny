// All-vs-all FastANI wrapper: stages all genomes so paths in the list files
// are local (work-dir relative). The nf-core FASTANI module cannot be used
// for this because it creates list files from workflow-level absolute paths,
// which are not accessible inside the container when inputs are remote URLs.
process FASTANI_ALLVSALL {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.34--hb66fcc3_7' :
        'quay.io/biocontainers/fastani:1.34--hb66fcc3_7' }"

    input:
    path(genomes)

    output:
    path("fastani.txt"), emit: ani
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    printf '%s\\n' ${genomes} | sort > genome_list.txt

    fastANI \\
        --ql genome_list.txt \\
        --rl genome_list.txt \\
        --threads ${task.cpus} \\
        -o fastani.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | head -1 | sed 's/version //')
    END_VERSIONS
    """

    stub:
    """
    touch fastani.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | head -1 | sed 's/version //')
    END_VERSIONS
    """
}
