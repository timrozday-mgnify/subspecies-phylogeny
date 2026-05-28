#!/usr/bin/env nextflow

include { SUBSPECIES_PHYLOGENY } from './workflows/subspecies_phylogeny'

workflow {
    main:
    // When resuming from a pre-computed merged SKF, no samplesheet is needed.
    ch_input = params.ska_merged_skf
        ? Channel.empty()
        : Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta  = [ id: row.sample ]
                def fasta = (row.fasta.startsWith('/') || row.fasta =~ /^[a-z]+:\/\//)
                    ? file(row.fasta, checkIfExists: true)
                    : file("${workflow.projectDir}/${row.fasta}", checkIfExists: true)
                [ meta, fasta ]
            }

    SUBSPECIES_PHYLOGENY(ch_input)
}
