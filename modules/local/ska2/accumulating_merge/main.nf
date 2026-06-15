// Accumulating ska merge with progressive weeding. Folds the input batch SKFs one
// at a time into a running accumulator, weeding the accumulator after each merge so
// it never holds positions that provably cannot reach the final min_freq across the
// full cohort. The weed threshold tightens as the accumulator grows toward `total`,
// keeping peak memory bounded throughout — the key advantage over a single final merge.
//
// The per-iteration min-freq is computed dynamically from the running genome count.
// Setting task.ext.args (e.g. '--min-freq 0.8') overrides the dynamic computation
// and applies that fixed weed argument after every merge instead.
process SKA2_ACCUMULATING_MERGE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska2:0.5.1--h4349ce8_0' :
        'quay.io/biocontainers/ska2:0.5.1--h4349ce8_0' }"

    input:
    path(skf_files, stageAs: '?/*')   // staged in numbered subdirs to avoid name collisions
    val(total)        // total genome count across all batches
    val(batch_size)   // nominal batch size, for the n_acc estimate
    val(min_freq)     // target final min-freq (0 = accumulate without weeding)

    output:
    path("merged.skf"),  emit: skf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // When set, ext.args overrides the dynamic per-iteration --min-freq.
    def args = task.ext.args ?: ''
    """
    printf '%s\\n' ${skf_files} > skf_filelist.txt
    n_files=\$(wc -l < skf_filelist.txt)
    echo "[SKA2_ACCUMULATING_MERGE] folding \$n_files batch SKFs; total=${total} batch_size=${batch_size} min_freq=${min_freq}"

    # Read file list into a positional array
    set -- \$(cat skf_filelist.txt)

    # Accumulator starts as the first batch SKF
    cp "\$1" accumulated.skf
    shift
    k=1

    while [ "\$#" -gt 0 ]; do
        next="\$1"; shift
        k=\$((k + 1))

        ska merge accumulated.skf "\$next" -o accumulated_new
        mv accumulated_new.skf accumulated.skf

        # n_acc = min(k * batch_size, total) — provably >= the true running count,
        # so the derived threshold under-weeds rather than over-weeds (safe).
        n_acc=\$(( k * ${batch_size} ))
        if [ "\$n_acc" -gt "${total}" ]; then n_acc=${total}; fi

        # Note: ska weed -o takes a literal filename (does not append .skf), unlike ska merge.
        if [ -n "$args" ]; then
            # User override: apply the fixed weed args every iteration
            echo "[SKA2_ACCUMULATING_MERGE] iter \$k: ska weed $args (override)"
            ska weed $args -o accumulated_weeded.skf accumulated.skf
            mv accumulated_weeded.skf accumulated.skf
        else
            # Dynamic threshold: mf = max(0, ceil(min_freq*total) - (total - n_acc)) / n_acc
            mf=\$(awk -v t=${total} -v na=\$n_acc -v f=${min_freq} 'BEGIN {
                mc = int(f * t); if (f * t > mc) mc++;
                mn = mc - (t - na); if (mn < 0) mn = 0;
                printf "%.6f", mn / na
            }')
            do_weed=\$(awk -v m=\$mf 'BEGIN { print (m > 0) ? 1 : 0 }')
            if [ "\$do_weed" -eq 1 ]; then
                echo "[SKA2_ACCUMULATING_MERGE] iter \$k: n_acc=\$n_acc ska weed --min-freq \$mf"
                ska weed --min-freq \$mf -o accumulated_weeded.skf accumulated.skf
                mv accumulated_weeded.skf accumulated.skf
            else
                echo "[SKA2_ACCUMULATING_MERGE] iter \$k: n_acc=\$n_acc weed skipped (threshold 0)"
            fi
        fi
    done

    mv accumulated.skf merged.skf

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
