include { SKA2_BUILD         } from '../../../modules/local/ska2/build/main'
include { SKA2_BATCHED_MERGE } from '../ska2_batched_merge/main'

workflow SKA2_PHYLOGENY {
    take:
    ch_input    // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    SKA2_BUILD(ch_input)
    ch_versions = ch_versions.mix(SKA2_BUILD.out.versions.first())

    // Pass individual SKF items (not collected) so the batched-merge subworkflow
    // can group them via .collate(batch_size) before merging.
    ch_skf = SKA2_BUILD.out.skf
        .map { meta, skf -> skf }

    SKA2_BATCHED_MERGE(ch_skf, params.ska_merge_batch_size as Integer)
    ch_versions = ch_versions.mix(SKA2_BATCHED_MERGE.out.versions)

    emit:
    skf        = SKA2_BUILD.out.skf            // channel: [ val(meta), path(*.skf) ] — one per sample
    merged_skf = SKA2_BATCHED_MERGE.out.skf    // channel: path(merged.skf)
    versions   = ch_versions
}
