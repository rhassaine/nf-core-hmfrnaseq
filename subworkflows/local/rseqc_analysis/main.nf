#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Constants
import Utils

// Import modules
include { RSEQC_BAMSTAT }         from '../../../modules/nf-core/rseqc/bamstat/main'
include { RSEQC_READDUPLICATION } from '../../../modules/nf-core/rseqc/readduplication/main'

workflow RSEQC_ANALYSIS {
    take:
    ch_inputs         // channel: [mandatory] [ meta ]
    ch_tumor_rna_bam  // channel: [mandatory] [ meta, bam, bai ]

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_tumor_rna_bam
        .map { meta, tumor_bam, tumor_bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tumor_bai, meta, Constants.INPUT.BAI_RNA_TUMOR),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.RSEQC_DIR)
            runnable: tumor_bam && !has_existing
            skip: true
                return meta
        }

    // Create process input channel for RSeQC modules
    // RSeQC modules expect [meta, bam] not [meta, bam, bai]
    // channel: [ meta_rseqc, tumor_bam ]
    ch_rseqc_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai ->
            def meta_sample = Utils.getTumorRnaSample(meta)
            def meta_rseqc = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: Utils.getTumorRnaSampleName(meta),
                subject_id: meta.subject_id,
                group_id: meta.group_id,
            ]

            return [meta_rseqc, tumor_bam]
        }

    // Run RSeQC modules
    RSEQC_BAMSTAT(ch_rseqc_inputs)
    ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions)

    RSEQC_READDUPLICATION(ch_rseqc_inputs)
    ch_versions = ch_versions.mix(RSEQC_READDUPLICATION.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, rseqc_outputs ]
    ch_bamstat_out = WorkflowOncoanalyser.restoreMeta(RSEQC_BAMSTAT.out.bamstat, ch_inputs)
    ch_readdup_seq_out = WorkflowOncoanalyser.restoreMeta(RSEQC_READDUPLICATION.out.seq_xls, ch_inputs)
    ch_readdup_pos_out = WorkflowOncoanalyser.restoreMeta(RSEQC_READDUPLICATION.out.pos_xls, ch_inputs)
    ch_readdup_pdf_out = WorkflowOncoanalyser.restoreMeta(RSEQC_READDUPLICATION.out.pdf, ch_inputs)
    ch_readdup_r_out = WorkflowOncoanalyser.restoreMeta(RSEQC_READDUPLICATION.out.rscript, ch_inputs)

    // Collect QC outputs for MultiQC
    ch_qc_reports = Channel.empty()
        .mix(
            ch_bamstat_out,
            ch_readdup_seq_out,
            ch_readdup_pos_out,
            ch_readdup_pdf_out,
            ch_readdup_r_out
        )

    // Add skipped samples with empty outputs
    ch_qc_reports_final = ch_qc_reports
        .mix(ch_inputs_sorted.skip.map { meta -> [meta, []] })

    emit:
    qc_reports = ch_qc_reports_final  // channel: [ meta, qc_files ]
    versions   = ch_versions          // channel: [ versions.yml ]
}