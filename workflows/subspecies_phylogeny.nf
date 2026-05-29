include { FASTANI_ALLVSALL            } from '../modules/local/fastani_allvsall/main'
include { SKA2_PHYLOGENY              } from '../subworkflows/local/ska2_phylogeny/main'
include { SKA2_ALIGN                  } from '../modules/local/ska2/align/main'
include { SKA2_DELETE                 } from '../modules/local/ska2/delete/main'
include { SKA2_DISTANCE               } from '../modules/local/ska2/distance/main'
include { SKA2_LO                     } from '../modules/local/ska2/lo/main'
include { SKA2_WEED                   } from '../modules/local/ska2/weed/main'
include { SKA2_MAP                    } from '../modules/local/ska2/map/main'
include { SELECT_REFERENCE            } from '../modules/local/select_reference/main'
include { NJ_TREE as NJ_TREE_FASTANI  } from '../modules/local/nj_tree/main'
include { NJ_TREE as NJ_TREE_SKA2     } from '../modules/local/nj_tree/main'
include { SNPSITES                    } from '../modules/nf-core/snpsites/main'
include { GUBBINS                     } from '../modules/nf-core/gubbins/main'
include { IQTREE                      } from '../modules/nf-core/iqtree/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/custom/dumpsoftwareversions/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'

workflow SUBSPECIES_PHYLOGENY {
    take:
    ch_input    // channel: [ val(meta), path(fasta) ] — empty when ska_merged_skf is set

    main:
    ch_versions = Channel.empty()

    // -----------------------------------------------------------------------
    // Upstream: either run the full BUILD → MERGE chain or skip straight to
    // alignment using a pre-computed merged SKF file.
    // -----------------------------------------------------------------------
    ch_nj_fastani  = Channel.empty()
    ch_map_reference = Channel.empty()

    if (params.ska_merged_skf) {
        ch_merged_skf = Channel.fromPath(params.ska_merged_skf, checkIfExists: true)
    } else {
        FASTANI_ALLVSALL(
            ch_input.map { meta, fasta -> fasta }.collect()
        )

        NJ_TREE_FASTANI(
            FASTANI_ALLVSALL.out.ani
                .map { f -> [ [id: 'fastani', format: 'fastani'], f ] }
        )
        ch_versions   = ch_versions.mix(NJ_TREE_FASTANI.out.versions)
        ch_nj_fastani = NJ_TREE_FASTANI.out.tree

        // Select the FastANI medoid as the ska map reference (unless overridden below).
        SELECT_REFERENCE(
            FASTANI_ALLVSALL.out.ani,
            ch_input.map { meta, fasta -> fasta }.collect()
        )
        ch_versions      = ch_versions.mix(SELECT_REFERENCE.out.versions)
        ch_map_reference = SELECT_REFERENCE.out.reference

        SKA2_PHYLOGENY(ch_input)
        ch_versions   = ch_versions.mix(SKA2_PHYLOGENY.out.versions)
        ch_merged_skf = SKA2_PHYLOGENY.out.merged_skf
    }

    // User-supplied reference takes priority over the auto-selected medoid.
    // In --ska_merged_skf mode without --ska_map_reference the channel stays
    // empty and the Gubbins track is skipped with a warning.
    ch_ska_map_ref = params.ska_map_reference
        ? Channel.fromPath(params.ska_map_reference, checkIfExists: true)
        : ch_map_reference

    // -----------------------------------------------------------------------
    // Optional SKA2_DELETE: remove specified samples from the merged SKF.
    // Applies to both pipeline-produced and user-supplied merged SKFs so
    // outliers can be excluded before alignment without rebuilding.
    // -----------------------------------------------------------------------
    if (params.ska_delete_samples) {
        ch_delete_file = Channel.fromPath(params.ska_delete_samples, checkIfExists: true)
        SKA2_DELETE(ch_merged_skf, ch_delete_file)
        ch_versions   = ch_versions.mix(SKA2_DELETE.out.versions)
        ch_merged_skf = SKA2_DELETE.out.skf
    }

    // -----------------------------------------------------------------------
    // Optional SKA2_DISTANCE: pairwise SNP distances from merged SKF.
    // -----------------------------------------------------------------------
    ch_distances  = Channel.empty()
    ch_nj_ska2    = Channel.empty()
    if (params.ska_distance) {
        SKA2_DISTANCE(ch_merged_skf)
        ch_versions  = ch_versions.mix(SKA2_DISTANCE.out.versions)
        ch_distances = SKA2_DISTANCE.out.distances

        NJ_TREE_SKA2(
            SKA2_DISTANCE.out.distances
                .map { f -> [ [id: 'ska2_distance', format: 'ska2'], f ] }
        )
        ch_versions = ch_versions.mix(NJ_TREE_SKA2.out.versions)
        ch_nj_ska2  = NJ_TREE_SKA2.out.tree
    }

    // -----------------------------------------------------------------------
    // Optional SKA2_LO: identify SNPs/INDELs left out by the split-kmer graph
    // (proxy for recombination). An optional reference FASTA can be provided
    // via params.ska_lo_reference to anchor coordinates.
    // -----------------------------------------------------------------------
    ch_lo_snps   = Channel.empty()
    ch_lo_indels = Channel.empty()
    if (params.ska_lo) {
        ch_lo_ref = params.ska_lo_reference
            ? Channel.fromPath(params.ska_lo_reference, checkIfExists: true)
            : Channel.value([])
        SKA2_LO(ch_merged_skf, ch_lo_ref)
        ch_versions  = ch_versions.mix(SKA2_LO.out.versions)
        ch_lo_snps   = SKA2_LO.out.snps
        ch_lo_indels = SKA2_LO.out.indels
    }

    // Pre-declare output channels so the emit block is always valid regardless
    // of which optional steps are enabled.
    ch_alignment = Channel.empty()
    ch_snp_sites = Channel.empty()
    ch_gubbins   = Channel.empty()
    ch_phylogeny = Channel.empty()

    if (!params.skip_alignment) {
        // -----------------------------------------------------------------------
        // SKA2_ALIGN: fan out over every requested --min-freq value.
        // params.ska_align_min_freq is a comma-separated string (e.g. "0.9" or
        // "0.5,0.9,1.0"). Each value becomes one independent analysis branch.
        // -----------------------------------------------------------------------
        ch_min_freq = Channel.fromList(
            params.ska_align_min_freq.tokenize(',').collect { it.trim() }
        )

        ch_align_input = ch_merged_skf
            .combine(ch_min_freq)
            .map { skf, mf -> [ [id: mf, min_freq: mf], skf ] }

        SKA2_ALIGN(ch_align_input)
        ch_versions  = ch_versions.mix(SKA2_ALIGN.out.versions.first())
        ch_alignment = SKA2_ALIGN.out.alignment

        // SNPSITES runs on every alignment branch for its published output
        // and to supply constant-sites counts for ascertainment-bias correction.
        SNPSITES(ch_alignment)
        ch_versions  = ch_versions.mix(SNPSITES.out.versions.first())
        ch_snp_sites = SNPSITES.out.fasta

        ch_gubbins = Channel.empty()
        if (!params.skip_gubbins) {
            // Weed each min_freq branch of the merged SKF (same frequency threshold
            // as ska align) then map against the selected reference genome.
            ch_weed_input = ch_merged_skf
                .combine(ch_min_freq)
                .map { skf, mf -> [ [id: mf, min_freq: mf], skf ] }

            SKA2_WEED(ch_weed_input)
            ch_versions = ch_versions.mix(SKA2_WEED.out.versions.first())

            // Pair each weeded SKF with the reference. When ch_ska_map_ref is
            // empty (ska_merged_skf mode without --ska_map_reference), the
            // combine produces no items and SKA2_MAP / GUBBINS simply never run.
            ch_map_input = SKA2_WEED.out.skf
                .combine(ch_ska_map_ref)
                .map { meta, skf, ref -> [ meta, ref, skf ] }

            SKA2_MAP(ch_map_input)
            ch_versions = ch_versions.mix(SKA2_MAP.out.versions.first())

            GUBBINS(SKA2_MAP.out.alignment)
            ch_versions = ch_versions.mix(GUBBINS.out.versions.first())
            ch_gubbins  = GUBBINS.out.fasta
        }

        if (!params.skip_iqtree) {
            // -----------------------------------------------------------------------
            // IQ-TREE: one or two tracks per min_freq combination.
            //   no_gubbins: snp-sites FASTA + ascertainment-bias correction (-fconst)
            //   gubbins:    Gubbins filtered_polymorphic_sites.fasta, no -fconst
            //               (only when skip_gubbins = false)
            // meta.gubbins carries 'no_gubbins' / 'gubbins' into the publishDir closure.
            // meta.constant_sites is only set on the no_gubbins track; its absence
            // makes the ext.args closure emit an empty string for the gubbins track.
            // -----------------------------------------------------------------------
            ch_iqtree_no_gubbins = SNPSITES.out.fasta
                .join(SNPSITES.out.constant_sites)
                .map { meta, fasta, cs_file ->
                    def cs = cs_file.text.trim()
                    [ meta + [gubbins: 'no_gubbins', constant_sites: cs], [fasta], [] ]
                }

            ch_iqtree_input = params.skip_gubbins
                ? ch_iqtree_no_gubbins
                : ch_iqtree_no_gubbins.mix(
                    ch_gubbins.map { meta, fasta -> [ meta + [gubbins: 'gubbins'], [fasta], [] ] }
                  )

            IQTREE(
                ch_iqtree_input,
                [], [], [], [], [], [], [], [], [], [], [], []
            )
            ch_versions  = ch_versions.mix(IQTREE.out.versions.first())
            ch_phylogeny = IQTREE.out.phylogeny
        }
    }

    // -----------------------------------------------------------------------
    // Software versions → MultiQC
    // -----------------------------------------------------------------------
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.collect())

    MULTIQC(
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml
            .collect()
            .map { files -> [ [ id: 'multiqc' ], files, [], [], [], [] ] }
    )

    emit:
    alignment  = ch_alignment   // per min-freq: [ val(meta), path(*.fasta) ]
    snp_sites  = ch_snp_sites   // per min-freq: [ val(meta), path(*.fas) ]
    gubbins    = ch_gubbins     // per min-freq: [ val(meta), path(*.filtered_polymorphic_sites.fasta) ]
    phylogeny  = ch_phylogeny   // per min-freq × gubbins track: [ val(meta), path(*.treefile) ]
    distances  = ch_distances   // path(distances.tsv)              — only when params.ska_distance
    lo_snps    = ch_lo_snps     // path(*_snps.fas)                 — only when params.ska_lo
    lo_indels  = ch_lo_indels   // path(*_indels.vcf)               — only when params.ska_lo
    nj_fastani = ch_nj_fastani  // [ val(meta), path(fastani.nwk) ] — only when !params.ska_merged_skf
    nj_ska2    = ch_nj_ska2     // [ val(meta), path(ska2_dist.nwk) ] — only when params.ska_distance
    multiqc    = MULTIQC.out.report
    versions   = ch_versions
}
