// Two-level batched merge subworkflow using the existing SKA2_MERGE module.
// Nextflow's .collate() fans out parallel batch merges (SKA2_MERGE_BATCH),
// then a final merge (SKA2_MERGE_FINAL) combines the batch outputs.
// Peak memory per task is bounded by the batch size rather than the full cohort.
include { SKA2_MERGE as SKA2_MERGE_BATCH } from '../../../modules/local/ska2/merge/main'
include { SKA2_MERGE as SKA2_MERGE_FINAL } from '../../../modules/local/ska2/merge/main'

workflow SKA2_BATCHED_MERGE {
    take:
    ch_skf      // channel: path(*.skf) — one item per sample (not pre-collected)
    batch_size  // val: int — number of SKFs per batch

    main:
    ch_versions = Channel.empty()

    // Group individual SKFs into batches and merge each batch in parallel
    ch_batches = ch_skf.collate(batch_size)

    SKA2_MERGE_BATCH(ch_batches)
    ch_versions = ch_versions.mix(SKA2_MERGE_BATCH.out.versions.first())

    // Collect all batch-level SKFs and run a single final merge
    SKA2_MERGE_FINAL(SKA2_MERGE_BATCH.out.skf.collect())
    ch_versions = ch_versions.mix(SKA2_MERGE_FINAL.out.versions)

    emit:
    skf      = SKA2_MERGE_FINAL.out.skf   // path(merged.skf)
    versions = ch_versions
}
