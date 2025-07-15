#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { RSEQC_BAMSTAT }         from '../../../modules/nf-core/rseqc/bamstat/main'
include { RSEQC_READDUPLICATION } from '../../../modules/nf-core/rseqc/readduplication/main'

workflow RSEQC_ANALYSIS {
    /*
     * Input channel: [meta, bam, bai]
     */
    take:
        ch_tumor_rna_bam // channel: [meta, bam, bai]

    main:
        // Ensure meta is always a valid map with id
        // ch_bam_for_rseqc = ch_tumor_rna_bam.map { meta, bam, bai ->
        //     def sample_id = meta?.id ?: meta?.sample_id ?: meta?.subject_id ?: meta?.group_id ?: 'unknown'
        //     if (!meta) meta = [id: sample_id]
        //     else meta.id = sample_id
        //     [meta, bam]
        // }

        ch_bam_for_rseqc = ch_tumor_rna_bam.map { meta, bam, bai ->
            def sample_id = meta?.id ?: meta?.sample_id ?: meta?.subject_id ?: meta?.group_id ?: 'unknown'
            def key_val = meta?.key ?: sample_id
            meta.key = key_val
            meta.id = sample_id
            [meta, bam]
        }

        // Run modules
        RSEQC_BAMSTAT(ch_bam_for_rseqc)
        RSEQC_READDUPLICATION(ch_bam_for_rseqc)

        // Collect version info
        ch_versions = Channel.empty()
            .mix(RSEQC_BAMSTAT.out.versions, RSEQC_READDUPLICATION.out.versions)

        // Collect QC outputs for MultiQC
        ch_qc_reports = Channel.empty()
            .mix(
                RSEQC_BAMSTAT.out.bamstat,
                RSEQC_READDUPLICATION.out.seq_xls,
                RSEQC_READDUPLICATION.out.pos_xls,
                RSEQC_READDUPLICATION.out.pdf,
                RSEQC_READDUPLICATION.out.rscript
            )
        // Downstream, always treat ch_qc_reports as tuple (meta, file)

    emit:
        bamstat    = RSEQC_BAMSTAT.out.bamstat
        seq_xls    = RSEQC_READDUPLICATION.out.seq_xls
        pos_xls    = RSEQC_READDUPLICATION.out.pos_xls
        pdf        = RSEQC_READDUPLICATION.out.pdf
        rscript    = RSEQC_READDUPLICATION.out.rscript
        versions   = ch_versions
        qc_reports = ch_qc_reports
}