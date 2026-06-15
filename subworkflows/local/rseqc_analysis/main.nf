#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Constants
import Utils
import WorkflowOncoanalyser

// Import modules
include { RSEQC_BAMSTAT }         from '../../../modules/nf-core/rseqc/bamstat/main'
include { RSEQC_READDUPLICATION } from '../../../modules/nf-core/rseqc/readduplication/main'
include { SAMTOOLS_RRNA_COUNT }   from '../../../modules/local/samtools/rrna_count/main'

workflow RSEQC_ANALYSIS {
    take:
        ch_inputs         // [meta]
        ch_tumor_rna_bam  // [meta, bam, bai]
        ch_bed            // [meta2, bed] (rRNA regions)

    main:
        // Sort inputs
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

        // Common per-sample meta ([meta, bam, bai])
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
                [meta_rseqc, tumor_bam, tumor_bai]
            }

        // Run RSeQC modules
        // Note: versions are collected via topics, no manual collection needed
        RSEQC_BAMSTAT(ch_rseqc_inputs)

        RSEQC_READDUPLICATION(ch_rseqc_inputs)

        // rRNA quantification via samtools interval overlap against the rRNA BED.
        // Same interval method as RSeQC split_bam, but a single samtools pass (no in/ex/junk
        // BAMs written) -> much faster. Emits a MultiQC custom-content YAML + a stats file.
        SAMTOOLS_RRNA_COUNT(ch_rseqc_inputs, ch_bed)

        // Restore meta for outputs
        ch_bamstat_out     = WorkflowOncoanalyser.restoreMeta(RSEQC_BAMSTAT.out.txt, ch_inputs)
        ch_readdup_seq_out = WorkflowOncoanalyser.restoreMeta(RSEQC_READDUPLICATION.out.seq_xls, ch_inputs)
        ch_readdup_pos_out = WorkflowOncoanalyser.restoreMeta(RSEQC_READDUPLICATION.out.pos_xls, ch_inputs)
        ch_readdup_pdf_out = WorkflowOncoanalyser.restoreMeta(RSEQC_READDUPLICATION.out.pdf, ch_inputs)
        ch_readdup_r_out   = WorkflowOncoanalyser.restoreMeta(RSEQC_READDUPLICATION.out.rscript, ch_inputs)

        ch_rrna_stats   = WorkflowOncoanalyser.restoreMeta(SAMTOOLS_RRNA_COUNT.out.stats, ch_inputs)
        ch_rrna_multiqc = WorkflowOncoanalyser.restoreMeta(SAMTOOLS_RRNA_COUNT.out.multiqc, ch_inputs)

        // Collect QC outputs for MultiQC (the rRNA custom-content YAML; not the stats txt)
        ch_qc_reports = channel.empty()
            .mix(
                ch_bamstat_out,
                ch_readdup_seq_out,
                ch_readdup_pos_out,
                ch_readdup_pdf_out,
                ch_readdup_r_out,
                ch_rrna_multiqc
            )

        // Add skipped samples with empty outputs
        ch_qc_reports_final = ch_qc_reports
            .mix(ch_inputs_sorted.skip.map { meta -> [meta, []] })

    emit:
        qc_reports = ch_qc_reports_final
        splitbam_stats = ch_rrna_stats  // [meta, stats_file] (rRNA counts; emit name kept for compatibility)
        // Note: versions are collected via topics, not emitted here
}
