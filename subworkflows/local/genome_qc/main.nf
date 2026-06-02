include { QUAST          } from '../../../modules/nf-core/quast/main'
include { MAGPURIFY      } from '../../../modules/local/magpurify/main'
include { CHECKM2_PREDICT } from '../../../modules/nf-core/checkm2/predict/main'
include { CHECKM         } from '../../../modules/local/checkm/main'
include { GUNC_RUN       } from '../../../modules/nf-core/gunc/run/main'
include { BUSCO_BUSCO    } from '../../../modules/nf-core/busco/busco/main'

workflow GENOME_QC {
    take:
    ch_input    // channel: [ val(meta), path(fasta) ] — may be empty in merged-SKF mode

    main:
    ch_versions = Channel.empty()

    // -----------------------------------------------------------------------
    // QUAST — assembly statistics; no database required, always runs.
    // -----------------------------------------------------------------------
    QUAST(
        ch_input.map { meta, fasta -> [ meta, [ fasta ] ] }
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    // -----------------------------------------------------------------------
    // MAGpurify — per-contig GC-content and tetranucleotide-frequency checks;
    // no database required.  phylo-markers module runs when magpurify_db is set.
    // -----------------------------------------------------------------------
    ch_magpurify_db = params.magpurify_db
        ? Channel.value(file(params.magpurify_db, checkIfExists: true))
        : Channel.value([])

    MAGPURIFY(
        ch_input,
        ch_magpurify_db
    )
    ch_versions = ch_versions.mix(MAGPURIFY.out.versions.first())

    // -----------------------------------------------------------------------
    // CheckM2 — ML-based completeness and contamination; requires database.
    // -----------------------------------------------------------------------
    if (params.checkm2_db) {
        ch_checkm2_db = Channel.value(
            [ [ id: 'checkm2_db' ], file(params.checkm2_db, checkIfExists: true) ]
        )
        CHECKM2_PREDICT(
            ch_input.map { meta, fasta -> [ meta, [ fasta ] ] },
            ch_checkm2_db
        )
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())
    }

    // -----------------------------------------------------------------------
    // CheckM — HMM-based completeness and contamination (legacy); requires database.
    // -----------------------------------------------------------------------
    if (params.checkm_db) {
        ch_checkm_db = Channel.value(
            file(params.checkm_db, checkIfExists: true)
        )
        CHECKM(ch_input, ch_checkm_db)
        ch_versions = ch_versions.mix(CHECKM.out.versions.first())
    }

    // -----------------------------------------------------------------------
    // GUNC — chimeric genome detection using gene-level taxonomy; requires database.
    // -----------------------------------------------------------------------
    if (params.gunc_db) {
        ch_gunc_db = Channel.value(
            file(params.gunc_db, checkIfExists: true)
        )
        GUNC_RUN(
            ch_input.map { meta, fasta -> [ meta, [ fasta ] ] },
            ch_gunc_db
        )
        ch_versions = ch_versions.mix(GUNC_RUN.out.versions.first())
    }

    // -----------------------------------------------------------------------
    // BUSCO — lineage-specific single-copy ortholog completeness.
    // Always runs.  When busco_lineage_db is set, that pre-downloaded dataset
    // is used (fast, reproducible, works offline).  Without it, BUSCO downloads
    // the lineage at runtime using --auto-lineage-prok (~20–60 min/genome,
    // requires internet access inside the container).
    // -----------------------------------------------------------------------
    ch_busco_db = params.busco_lineage_db
        ? Channel.value(file(params.busco_lineage_db, checkIfExists: true))
        : Channel.value([])

    BUSCO_BUSCO(
        ch_input,
        params.busco_lineage ?: 'auto_prok',
        ch_busco_db
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    emit:
    versions = ch_versions
}
