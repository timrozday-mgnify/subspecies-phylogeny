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

    // Pre-merge weed threshold: defaults to the loosest --min-freq used downstream
    // by ska align / ska weed (params.ska_align_min_freq), so the accumulating
    // merge drops positions that no downstream branch can ever need. The minimum
    // across all requested align branches keeps this conservative. Set
    // params.ska_merge_min_freq > 0 to override with a fixed value, or rely on
    // params.skip_alignment to leave merged.skf unweeded.
    def merge_min_freq = params.ska_merge_min_freq as Double
    if (merge_min_freq <= 0 && !params.skip_alignment) {
        merge_min_freq = params.ska_align_min_freq
            .tokenize(',')
            .collect { it.trim() as Double }
            .min()
    }

    SKA2_BATCHED_MERGE(ch_skf, params.ska_merge_batch_size as Integer, merge_min_freq)
    ch_versions = ch_versions.mix(SKA2_BATCHED_MERGE.out.versions)

    emit:
    skf        = SKA2_BUILD.out.skf            // channel: [ val(meta), path(*.skf) ] — one per sample
    merged_skf = SKA2_BATCHED_MERGE.out.skf    // channel: path(merged.skf)
    versions   = ch_versions
}
