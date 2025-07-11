#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import the nf-core RSeQC bamstat module
include { RSEQC_BAMSTAT } from '../../../modules/nf-core/rseqc/bamstat/main'
// include { ISOFOX } from '../../../modules/local/isofox/main'

workflow RSEQC_ANALYSIS {
    /*
     * Input channel: BAM files for QC
     */
    take:
        ch_bam

    main:
        // Run RSeQC bamstat module
        RSEQC_BAMSTAT(ch_bam)

        // Collect version info
        versions_ch = Channel.empty().mix(RSEQC_BAMSTAT.out.versions)

    emit:
        bamstat = RSEQC_BAMSTAT.out.bamstat
        versions = versions_ch
}