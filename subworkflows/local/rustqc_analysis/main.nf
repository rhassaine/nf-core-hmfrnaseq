#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Constants
import Utils
import WorkflowOncoanalyser

include { RUSTQC } from '../../../modules/nf-core/rustqc/main'

workflow RUSTQC_ANALYSIS {
    take:
        ch_inputs         // [meta]
        ch_tumor_rna_bam  // [meta, bam, bai]
        ch_gtf            // [meta2, gtf]

    main:
        // Sort inputs: branch into runnable (has BAM, no existing results) vs skip
        ch_inputs_sorted = ch_tumor_rna_bam
            .map { meta, tumor_bam, tumor_bai ->
                [meta,
                 Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                 Utils.selectCurrentOrExisting(tumor_bai, meta, Constants.INPUT.BAI_RNA_TUMOR)]
            }
            .branch { meta, tumor_bam, tumor_bai ->
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.RSEQC_DIR)
                def has_bam = tumor_bam && !(tumor_bam instanceof List && tumor_bam.isEmpty())
                runnable: has_bam && !has_existing
                skip: true
                    return meta
            }

        // Build RustQC input channel with process-specific meta
        ch_rustqc_inputs = ch_inputs_sorted.runnable
            .map { meta, tumor_bam, tumor_bai ->
                def meta_sample = Utils.getTumorRnaSample(meta)
                def meta_rustqc = [
                    key: meta.group_id,
                    id: "${meta.group_id}_${meta_sample.sample_id}",
                    sample_id: Utils.getTumorRnaSampleName(meta),
                    subject_id: meta.subject_id,
                    group_id: meta.group_id,
                ]
                [meta_rustqc, tumor_bam, tumor_bai]
            }

        // Run RustQC — single-pass QC replacing bamstat, readduplication, and splitbam
        // Note: versions are collected via topics, no manual collection needed
        RUSTQC(ch_rustqc_inputs, ch_gtf)

        // Restore original meta on outputs
        ch_featurecounts = WorkflowOncoanalyser.restoreMeta(RUSTQC.out.featurecounts, ch_inputs)
        ch_rseqc_out     = WorkflowOncoanalyser.restoreMeta(RUSTQC.out.rseqc, ch_inputs)
        ch_samtools_out  = WorkflowOncoanalyser.restoreMeta(RUSTQC.out.samtools, ch_inputs)
        ch_dupradar_out  = WorkflowOncoanalyser.restoreMeta(RUSTQC.out.dupradar, ch_inputs)
        ch_preseq_out    = WorkflowOncoanalyser.restoreMeta(RUSTQC.out.preseq, ch_inputs)
        ch_qualimap_out  = WorkflowOncoanalyser.restoreMeta(RUSTQC.out.qualimap, ch_inputs)

        // Extract biotype_counts.tsv from featurecounts output for rRNA QC gating
        ch_biotype_counts = ch_featurecounts
            .map { meta, files ->
                def biotype_file = files.find { it.name.endsWith('.biotype_counts.tsv') }
                [meta, biotype_file]
            }

        // Collect all QC outputs for MultiQC
        ch_qc_reports = channel.empty()
            .mix(
                ch_featurecounts,
                ch_rseqc_out,
                ch_samtools_out,
                ch_dupradar_out,
                ch_preseq_out,
                ch_qualimap_out,
            )

        // Add skipped samples with empty outputs
        ch_qc_reports_final = ch_qc_reports
            .mix(ch_inputs_sorted.skip.map { meta -> [meta, []] })

    emit:
        qc_reports     = ch_qc_reports_final   // [meta, files] for MultiQC
        biotype_counts = ch_biotype_counts     // [meta, biotype_counts_file] for rRNA QC gate
}
