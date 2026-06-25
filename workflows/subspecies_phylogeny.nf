include { GENOME_QC                   } from '../subworkflows/local/genome_qc/main'
include { FASTANI_ALLVSALL            } from '../modules/local/fastani_allvsall/main'
include { SKA2_PHYLOGENY              } from '../subworkflows/local/ska2_phylogeny/main'
include { SKA2_ALIGN                  } from '../modules/local/ska2/align/main'
include { SKA2_DELETE                 } from '../modules/local/ska2/delete/main'
include { SKA2_DISTANCE               } from '../modules/local/ska2/distance/main'
include { SKA2_LO                     } from '../modules/local/ska2/lo/main'
include { SKA2_WEED                   } from '../modules/local/ska2/weed/main'
include { SKA2_SUBSET                 } from '../modules/local/ska2/subset/main'
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
    // Genome QC — runs before any phylogenetic analysis.
    // Results are informational only; ch_input is not filtered.
    // Individual database-dependent tools (CheckM2, CheckM, GUNC, BUSCO) only
    // run when their respective database path parameter is non-null.
    // Disable with --skip_qc.
    // -----------------------------------------------------------------------
    if (!params.skip_qc) {
        GENOME_QC(ch_input)
        ch_versions = ch_versions.mix(GENOME_QC.out.versions)
    }

    // Pre-declare output channels so the emit block is always valid regardless
    // of which optional steps are enabled.
    ch_alignment  = Channel.empty()
    ch_snp_sites  = Channel.empty()
    ch_gubbins    = Channel.empty()
    ch_phylogeny  = Channel.empty()
    ch_distances  = Channel.empty()
    ch_nj_fastani = Channel.empty()
    ch_nj_ska2    = Channel.empty()
    ch_lo_snps    = Channel.empty()
    ch_lo_indels  = Channel.empty()

    if (!params.skip_phylo) {

    if (params.ska_gubbins_subset) {
        if (params.ska_gubbins_subset_target_snps == null) {
            throw new IllegalArgumentException("--ska_gubbins_subset requires --ska_gubbins_subset_target_snps")
        }
        if ((params.ska_gubbins_subset_target_snps as Integer) <= 0) {
            throw new IllegalArgumentException("--ska_gubbins_subset_target_snps must be a positive integer")
        }
    }

    // -----------------------------------------------------------------------
    // Upstream: either run the full BUILD → MERGE chain or skip straight to
    // alignment using a pre-computed merged SKF file.
    // -----------------------------------------------------------------------
    ch_map_reference = Channel.empty()

    // The Gubbins track needs a reference for ska map. It is resolved in this
    // order: (1) --ska_map_reference override, (2) the ska2 SNP-distance medoid of
    // the input genomes (lowest mean pairwise distance). (2) needs the input
    // genomes themselves to emit a reference FASTA, so in --ska_merged_skf mode it
    // only applies when the samplesheet is still supplied via --input.
    gubbins_track_active = !params.skip_gubbins && !params.skip_alignment

    // -----------------------------------------------------------------------
    // Obtain the merged SKF: reuse a pre-computed one, or run BUILD → MERGE.
    // -----------------------------------------------------------------------
    if (params.ska_merged_skf) {
        ch_merged_skf = Channel.fromPath(params.ska_merged_skf, checkIfExists: true)
    } else {
        SKA2_PHYLOGENY(ch_input)
        ch_versions   = ch_versions.mix(SKA2_PHYLOGENY.out.versions)
        ch_merged_skf = SKA2_PHYLOGENY.out.merged_skf
    }

    // -----------------------------------------------------------------------
    // Optional SKA2_DELETE: remove specified samples from the merged SKF before
    // distances / reference selection, so the reference is chosen from retained
    // samples. Applies to both pipeline-produced and user-supplied merged SKFs.
    // -----------------------------------------------------------------------
    if (params.ska_delete_samples) {
        ch_delete_file = Channel.fromPath(params.ska_delete_samples, checkIfExists: true)
        SKA2_DELETE(ch_merged_skf, ch_delete_file)
        ch_versions   = ch_versions.mix(SKA2_DELETE.out.versions)
        ch_merged_skf = SKA2_DELETE.out.skf
    }

    // -----------------------------------------------------------------------
    // ska2 pairwise SNP distances → ska map reference (medoid) and/or QC NJ tree.
    // Computed when the user asks for them (--ska_distance) or when needed to
    // auto-select the reference.
    // -----------------------------------------------------------------------
    // Candidate genomes for the reference. With --ska_map_ref_trusted_only, only
    // trusted=true genomes are eligible; otherwise all input genomes are.
    def ch_ref_fastas = params.ska_map_ref_trusted_only
        ? ch_input.filter { meta, fasta -> meta.trusted }.map { meta, fasta -> fasta }.collect()
        : ch_input.map { meta, fasta -> fasta }.collect()

    // Auto-select a reference only when the Gubbins track is active, no explicit
    // override was given, and input genomes are available to emit it from.
    ref_auto_needed = gubbins_track_active && !params.ska_map_reference && params.input

    ch_distances = Channel.empty()
    ch_nj_ska2   = Channel.empty()
    if (params.ska_distance || ref_auto_needed) {
        SKA2_DISTANCE(ch_merged_skf)
        ch_versions  = ch_versions.mix(SKA2_DISTANCE.out.versions)
        ch_distances = SKA2_DISTANCE.out.distances

        if (params.ska_distance) {
            NJ_TREE_SKA2(
                SKA2_DISTANCE.out.distances
                    .map { f -> [ [id: 'ska2_distance', format: 'ska2'], f ] }
            )
            ch_versions = ch_versions.mix(NJ_TREE_SKA2.out.versions)
            ch_nj_ska2  = NJ_TREE_SKA2.out.tree
        }

        if (ref_auto_needed) {
            // Pick the ska2 SNP-distance medoid (lowest mean pairwise distance).
            SELECT_REFERENCE(SKA2_DISTANCE.out.distances, ch_ref_fastas)
            ch_versions      = ch_versions.mix(SELECT_REFERENCE.out.versions)
            ch_map_reference = SELECT_REFERENCE.out.reference
        }
    }

    // -----------------------------------------------------------------------
    // Optional FastANI all-vs-all QC (ANI matrix + NJ tree). No longer used for
    // reference selection; enable with --run_fastani. Needs the input genomes.
    // -----------------------------------------------------------------------
    if (params.run_fastani && params.input) {
        FASTANI_ALLVSALL(
            ch_input.map { meta, fasta -> fasta }.collect()
        )
        ch_versions = ch_versions.mix(FASTANI_ALLVSALL.out.versions)

        NJ_TREE_FASTANI(
            FASTANI_ALLVSALL.out.ani
                .map { f -> [ [id: 'fastani', format: 'fastani'], f ] }
        )
        ch_versions   = ch_versions.mix(NJ_TREE_FASTANI.out.versions)
        ch_nj_fastani = NJ_TREE_FASTANI.out.tree
    }

    // User-supplied reference takes priority over the auto-selected medoid.
    // In --ska_merged_skf mode without --ska_map_reference the channel stays
    // empty and the Gubbins track is skipped with a warning.
    ch_ska_map_ref = params.ska_map_reference
        ? Channel.fromPath(params.ska_map_reference, checkIfExists: true)
        : ch_map_reference

    // Warn loudly when the Gubbins track will silently produce nothing because
    // no reference could be resolved: merged-SKF mode, no --ska_map_reference
    // override, and no --input genomes to compute a distance medoid from.
    if (params.ska_merged_skf && !params.ska_map_reference && !params.input && gubbins_track_active) {
        log.warn(
            "SKA2_MAP and GUBBINS will be skipped: --ska_merged_skf is set without " +
            "--ska_map_reference, and no --input samplesheet was provided to auto-select " +
            "a ska2 SNP-distance medoid reference. To run the Gubbins track, set " +
            "--ska_map_reference <ref.fasta>, or pass --input <samplesheet.csv>."
        )
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

            // A high enough --min-freq can weed out every k-mer, leaving an SKF
            // with no SNPs that would otherwise crash ska map / Gubbins downstream.
            ch_weeded_skf = SKA2_WEED.out.skf
                .filter { meta, skf, nkmers_file ->
                    def keep = nkmers_file.text.trim().toInteger() >= 10
                    if (!keep) {
                        log.warn(
                            "Skipping ska map/Gubbins for min_freq=${meta.min_freq}: " +
                            "After ska weed there are too few SNPs to map (<10 remaining)."
                        )
                    }
                    keep
                }

            ch_gubbins_skf = ch_weeded_skf
                .map { meta, skf, nkmers_file -> [ meta, skf ] }
            if (params.ska_gubbins_subset) {
                SKA2_SUBSET(ch_weeded_skf, params.ska_gubbins_subset_target_snps as Integer)
                ch_versions = ch_versions.mix(SKA2_SUBSET.out.versions.first())

                ch_gubbins_skf = SKA2_SUBSET.out.skf
                    .map { meta, skf, nkmers_file, sparsity_file -> [ meta, skf ] }
            }

            // Pair each weeded SKF with the reference. When ch_ska_map_ref is
            // empty (ska_merged_skf mode without --ska_map_reference), the
            // combine produces no items and SKA2_MAP / GUBBINS simply never run.
            ch_map_input = ch_gubbins_skf
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

    } // end !params.skip_phylo

    // -----------------------------------------------------------------------
    // Software versions → MultiQC (always run)
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
    distances  = ch_distances   // path(distances.tsv)              — when ska_distance or reference auto-select
    lo_snps    = ch_lo_snps     // path(*_snps.fas)                 — only when params.ska_lo
    lo_indels  = ch_lo_indels   // path(*_indels.vcf)               — only when params.ska_lo
    nj_fastani = ch_nj_fastani  // [ val(meta), path(fastani.nwk) ] — only when params.run_fastani
    nj_ska2    = ch_nj_ska2     // [ val(meta), path(ska2_dist.nwk) ] — only when params.ska_distance
    multiqc    = MULTIQC.out.report
    versions   = ch_versions
}
