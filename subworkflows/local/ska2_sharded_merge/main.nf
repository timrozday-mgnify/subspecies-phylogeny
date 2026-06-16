// Sharded merge subworkflow — a memory-bounded alternative to SKA2_BATCHED_MERGE.
//
// Each per-sample SKF is split into n_shards bins by hashing the full split-kmer
// flank (SKA2_SHARD_SPLIT). Because split-kmer keys are canonical, the same k-mer
// goes to the same bin in every sample, so each bin can be merged across samples
// independently (SKA2_MERGE_SHARD) — every per-bin merge touches only ~1/n_shards
// of the key space, bounding peak memory and running fully in parallel. The merged
// bins are then concatenated back into one merged.skf (SKA2_SHARD_CONCAT).
//
// When min_freq > 0 each merged bin is weeded (SKA2_WEED_SHARD) before concat.
// A merged bin already holds every sample's column for its k-mers, so weeding it
// at the final min_freq is exact (not merely conservative) and matches the
// downstream ska align / ska weed filter — consistent with the batched-merge
// published output.
include { SKA2_SHARD_SPLIT                      } from '../../../modules/local/ska2/shard_split/main'
include { SKA2_MERGE        as SKA2_MERGE_SHARD } from '../../../modules/local/ska2/merge/main'
include { SKA2_WEED         as SKA2_WEED_SHARD  } from '../../../modules/local/ska2/weed/main'
include { SKA2_SHARD_CONCAT                     } from '../../../modules/local/ska2/shard_concat/main'

workflow SKA2_SHARDED_MERGE {
    take:
    ch_skf    // channel: path(*.skf) — one per sample (bare path, like SKA2_BATCHED_MERGE)
    n_shards  // val: int — number of bins
    min_freq  // val: double — final weed threshold (0 = no weed)

    main:
    ch_versions = Channel.empty()

    // Synthesise a per-sample id from the filename so shard outputs are distinct.
    ch_split_in = ch_skf.map { skf -> [ [id: skf.baseName], skf ] }
    SKA2_SHARD_SPLIT(ch_split_in, n_shards)
    ch_versions = ch_versions.mix(SKA2_SHARD_SPLIT.out.versions.first())

    // Regroup shards by bin index across all samples (parse trailing `.<i>.skf`),
    // carrying the sample id alongside each shard file.
    ch_bins = SKA2_SHARD_SPLIT.out.bins
        .flatMap { meta, files ->
            files.collect { f ->
                def idx = (f.name =~ /\.(\d+)\.skf$/)[0][1] as Integer
                tuple(idx, tuple(meta.id, f))
            }
        }
        .groupTuple()
        // Sort each bin's files by sample id so EVERY bin is merged in the same
        // sample order. groupTuple collects in nondeterministic arrival order, and
        // `ska merge` writes its `names` columns in input order — without this, bins
        // would get inconsistent sample orderings and SKA2_SHARD_CONCAT (which stacks
        // their variant matrices column-for-column) would reject them.
        .map { idx, pairs -> pairs.sort { it[0] }.collect { it[1] } }

    // Merge each bin across samples (~1/n_shards of the key space per task).
    SKA2_MERGE_SHARD( ch_bins )
    ch_versions = ch_versions.mix(SKA2_MERGE_SHARD.out.versions.first())

    if (min_freq > 0) {
        ch_weed_in = SKA2_MERGE_SHARD.out.skf
            .map { skf -> [ [id: 'shard', min_freq: min_freq], skf ] }

        SKA2_WEED_SHARD(ch_weed_in)
        ch_versions = ch_versions.mix(SKA2_WEED_SHARD.out.versions.first())

        ch_for_concat = SKA2_WEED_SHARD.out.skf
            .map { meta, skf -> skf }
            .collect()
    } else {
        ch_for_concat = SKA2_MERGE_SHARD.out.skf.collect()
    }

    SKA2_SHARD_CONCAT(ch_for_concat)
    ch_versions = ch_versions.mix(SKA2_SHARD_CONCAT.out.versions)

    emit:
    skf      = SKA2_SHARD_CONCAT.out.skf   // path(merged.skf) — same contract as SKA2_BATCHED_MERGE
    versions = ch_versions
}
