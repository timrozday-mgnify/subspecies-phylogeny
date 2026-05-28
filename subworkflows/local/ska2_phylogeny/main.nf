include { SKA2_BUILD } from '../../../modules/local/ska2/build/main'
include { SKA2_MERGE } from '../../../modules/local/ska2/merge/main'

workflow SKA2_PHYLOGENY {
    take:
    ch_input    // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    SKA2_BUILD(ch_input)
    ch_versions = ch_versions.mix(SKA2_BUILD.out.versions.first())

    ch_skf = SKA2_BUILD.out.skf
        .map { meta, skf -> skf }
        .collect()

    SKA2_MERGE(ch_skf)
    ch_versions = ch_versions.mix(SKA2_MERGE.out.versions)

    emit:
    skf        = SKA2_BUILD.out.skf    // channel: [ val(meta), path(*.skf) ] — one per sample
    merged_skf = SKA2_MERGE.out.skf   // channel: path(*.skf)
    versions   = ch_versions
}
