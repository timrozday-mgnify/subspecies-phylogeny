#!/usr/bin/env nextflow

include { SUBSPECIES_PHYLOGENY } from './workflows/subspecies_phylogeny'

workflow {
    main:
    // Parse the samplesheet whenever it is provided. In the normal (build) mode
    // it is required. In --ska_merged_skf mode it is optional: the build/merge
    // steps are skipped, but if the genomes are still supplied they are used to
    // auto-select the FastANI medoid as the ska map reference for the Gubbins
    // track (unless overridden by --ska_map_reference).
    ch_input = params.input
        ? Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def trusted = row.containsKey('trusted') ? (row.trusted?.toLowerCase() in ['true', 'yes', '1']) : false
                def meta    = [ id: row.sample, trusted: trusted ]
                def fasta   = (row.fasta.startsWith('/') || row.fasta =~ /^[a-z]+:\/\//)
                    ? file(row.fasta, checkIfExists: true)
                    : file("${workflow.projectDir}/${row.fasta}", checkIfExists: true)
                [ meta, fasta ]
            }
        : Channel.empty()

    SUBSPECIES_PHYLOGENY(ch_input)
}
