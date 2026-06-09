// Batch-aware ska merge. Splits input SKFs into chunks of batch_size, merges each
// chunk, then progressively accumulates: merge(accumulated, chunk) → accumulated.
// This bounds peak memory to batch_size × per-sample SKF size rather than the full
// cohort, making large datasets tractable without sacrificing correctness
// (ska merge is associative).
process SKA2_BATCHED_MERGE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska2:0.5.1--h4349ce8_0' :
        'quay.io/biocontainers/ska2:0.5.1--h4349ce8_0' }"

    input:
    path(skf_files)
    val(batch_size)

    output:
    path("merged.skf"),  emit: skf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    printf '%s\\n' ${skf_files} > skf_filelist.txt
    total=\$(wc -l < skf_filelist.txt)
    echo "[SKA2_BATCHED_MERGE] total SKFs: \$total  batch_size: ${batch_size}"

    if [ "\$total" -le "${batch_size}" ]; then
        # All files fit in one batch — single merge
        echo "[SKA2_BATCHED_MERGE] batch 1/1 (\$total files) -> merged.skf"
        xargs -x -a skf_filelist.txt ska merge $args -o merged
    else
        mkdir -p chunks
        split -l ${batch_size} skf_filelist.txt chunks/chunk_
        n_batches=\$(ls chunks/chunk_* | wc -l)
        echo "[SKA2_BATCHED_MERGE] \$n_batches batches of up to ${batch_size} files each"

        # Merge the first chunk into the accumulator
        first=\$(ls chunks/chunk_* | head -n 1)
        first_count=\$(wc -l < "\$first")
        echo "[SKA2_BATCHED_MERGE] batch 1/\$n_batches (\$first_count files) -> accumulated.skf"
        xargs -x -a "\$first" ska merge $args -o accumulated
        rm "\$first"

        # Append each remaining chunk into the accumulator
        batch_num=2
        for chunk in \$(ls chunks/chunk_*); do
            chunk_count=\$(wc -l < "\$chunk")
            echo "[SKA2_BATCHED_MERGE] batch \$batch_num/\$n_batches (\$chunk_count files) -> accumulated.skf"
            xargs -x -a "\$chunk" ska merge $args accumulated.skf -o accumulated_new
            mv accumulated_new.skf accumulated.skf
            rm "\$chunk"
            batch_num=\$((batch_num + 1))
        done

        mv accumulated.skf merged.skf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """

    stub:
    """
    touch merged.skf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska2: \$(ska --version 2>&1 | sed 's/ska //')
    END_VERSIONS
    """
}
