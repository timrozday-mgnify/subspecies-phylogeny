include { SKA2_BUILD         } from '../../../modules/local/ska2/build/main'
include { SKA2_SHARDED_MERGE } from '../ska2_sharded_merge/main'

workflow SKA2_PHYLOGENY {
    take:
    ch_input    // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    SKA2_BUILD(ch_input)
    ch_versions = ch_versions.mix(SKA2_BUILD.out.versions.first())

    // Pass individual SKF items (not collected) so the sharded-merge subworkflow
    // can split each one by minimizer before merging per bin.
    ch_skf = SKA2_BUILD.out.skf
        .map { meta, skf -> skf }

    // Per-bin weed threshold: defaults to the loosest --min-freq used downstream
    // by ska align / ska weed (params.ska_align_min_freq), so each merged shard
    // drops positions that no downstream branch can ever need. The minimum across
    // all requested align branches keeps this safe. Set params.ska_merge_min_freq
    // > 0 to override with a fixed value, or rely on params.skip_alignment to
    // leave merged.skf unweeded.
    def merge_min_freq = params.ska_merge_min_freq as Double
    if (merge_min_freq <= 0 && !params.skip_alignment) {
        merge_min_freq = params.ska_align_min_freq
            .tokenize(',')
            .collect { it.trim() as Double }
            .min()
    }

    // Sharded merge: split each sample by minimizer into params.ska_shard_count
    // bins, merge each bin across samples (~1/n of the key space per task), weed
    // each merged bin at merge_min_freq, then concatenate the bins into merged.skf.
    SKA2_SHARDED_MERGE(
        ch_skf,
        params.ska_shard_count as Integer,
        params.ska_shard_minimizer_len as Integer,
        merge_min_freq
    )
    ch_versions = ch_versions.mix(SKA2_SHARDED_MERGE.out.versions)

    emit:
    skf        = SKA2_BUILD.out.skf            // channel: [ val(meta), path(*.skf) ] — one per sample
    merged_skf = SKA2_SHARDED_MERGE.out.skf    // channel: path(merged.skf)
    versions   = ch_versions
}
