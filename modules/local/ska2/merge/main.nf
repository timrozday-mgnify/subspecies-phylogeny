process SKA2_MERGE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska2:0.5.1--h4349ce8_0' :
        'quay.io/biocontainers/ska2:0.5.1--h4349ce8_0' }"

    input:
    path(skf_files, stageAs: '?/*')

    output:
    path("merged.skf"),  emit: skf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Bound every `ska merge` invocation to at most this many files on its command
    // line. `xargs` without an explicit `-n`/`-s` chooses its own batch size and,
    // for file lists too long for one command line, silently runs the command
    // *multiple times* — and since every invocation here writes to the same
    // `-o merged`, each one overwrote the previous one's output, silently
    // dropping every sample not in the last batch (observed in production with a
    // 1879-file cohort: two independent "SKA done" runs in one task, final
    // merged.skf left with only the last batch's ~700 samples). Folding batches
    // into an accumulator instead of relying on xargs's auto-split avoids this
    // regardless of cohort size. Overridable per-test via task.ext.merge_batch_size;
    // 100 matches the already-tuned params.ska_merge_batch_size used elsewhere.
    def batch_size = task.ext.merge_batch_size ?: 100
    """
    printf '%s\\n' ${skf_files} > skf_filelist.txt
    total=\$(wc -l < skf_filelist.txt)

    if [ "\$total" -eq 1 ]; then
        cp "\$(cat skf_filelist.txt)" merged.skf
    elif [ "\$total" -le ${batch_size} ]; then
        # Single invocation is safe: argv is bounded by ${batch_size} files.
        xargs -x -a skf_filelist.txt ska merge $args -o merged
    else
        # Fold the file list in batches so no single `ska merge` invocation ever
        # receives more than ${batch_size} files on its command line.
        mkdir -p chunks
        split -l ${batch_size} skf_filelist.txt chunks/chunk_

        first=\$(ls chunks/chunk_* | head -n 1)
        xargs -x -a "\$first" ska merge $args -o accumulated
        rm "\$first"

        for chunk in \$(ls chunks/chunk_* 2>/dev/null); do
            xargs -x -a "\$chunk" ska merge $args accumulated.skf -o accumulated_new
            mv accumulated_new.skf accumulated.skf
            rm "\$chunk"
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
