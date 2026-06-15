// Two-level batched merge subworkflow.
// Nextflow's .collate() fans out parallel batch merges (SKA2_MERGE_BATCH), then an
// accumulating merge (SKA2_ACCUMULATING_MERGE) folds the batch outputs one at a time
// into a running accumulator, weeding after each step. Peak memory per task is bounded
// by the batch size rather than the full cohort.
//
// When min_freq > 0:
//   - each batch SKF is weeded (SKA2_WEED_BATCH) before the fold, removing positions
//     that provably cannot reach min_freq in the full cohort;
//   - the accumulating merge tightens its weed threshold as the accumulator grows.
// Both weeds are conservative (never remove a position that could meet the final
// threshold) — they are memory pre-filters ahead of the real ska_align_min_freq filter.
//
// When min_freq = 0 the accumulating merge is a plain sequential fold (no weeding):
// same memory profile as a single merge but sequential. Memory savings come from the
// weeding, i.e. the large-cohort min_freq > 0 case.
//
// Batch weed formula: batch_min_freq = max(0, (min_freq × total − (total − batch_size)) / batch_size)
include { SKA2_MERGE              as SKA2_MERGE_BATCH } from '../../../modules/local/ska2/merge/main'
include { SKA2_WEED               as SKA2_WEED_BATCH  } from '../../../modules/local/ska2/weed/main'
include { SKA2_ACCUMULATING_MERGE                     } from '../../../modules/local/ska2/accumulating_merge/main'

workflow SKA2_BATCHED_MERGE {
    take:
    ch_skf      // channel: path(*.skf) — one item per sample (not pre-collected)
    batch_size  // val: int — number of SKFs per batch
    min_freq    // val: double — pre-merge weed threshold (0 = disabled)

    main:
    ch_versions = Channel.empty()

    // Total genome count (value channel — broadcasts, reused below)
    ch_total = ch_skf.count()

    // Group individual SKFs into batches and merge each batch in parallel
    ch_batches = ch_skf.collate(batch_size)

    SKA2_MERGE_BATCH(ch_batches)
    ch_versions = ch_versions.mix(SKA2_MERGE_BATCH.out.versions.first())

    if (min_freq > 0) {
        // Per-batch weed threshold from the total genome count. Uses nominal
        // batch_size — slightly conservative for a smaller last batch, but never
        // discards a position that could meet the final threshold.
        ch_batch_mf = ch_total.map { total ->
            def min_count    = Math.ceil(min_freq * total) as int
            def worst_others = total - batch_size
            def min_needed   = Math.max(0, min_count - worst_others)
            Math.min(1.0, (min_needed as double) / batch_size)
        }

        ch_weed_input = SKA2_MERGE_BATCH.out.skf
            .combine(ch_batch_mf)
            .map { skf, bmf -> [ [id: 'batch', min_freq: bmf], skf ] }

        SKA2_WEED_BATCH(ch_weed_input)
        ch_versions = ch_versions.mix(SKA2_WEED_BATCH.out.versions.first())

        ch_for_final = SKA2_WEED_BATCH.out.skf
            .map { meta, skf -> skf }
            .collect()
    } else {
        ch_for_final = SKA2_MERGE_BATCH.out.skf.collect()
    }

    // Sequentially fold the batch SKFs into one accumulator, weeding after each step.
    SKA2_ACCUMULATING_MERGE(ch_for_final, ch_total, batch_size, min_freq)
    ch_versions = ch_versions.mix(SKA2_ACCUMULATING_MERGE.out.versions)

    emit:
    skf      = SKA2_ACCUMULATING_MERGE.out.skf   // path(merged.skf)
    versions = ch_versions
}
