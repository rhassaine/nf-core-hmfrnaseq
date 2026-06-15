#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Constants
import Utils
import WorkflowOncoanalyser

include { SAMTOOLS_RRNA_COUNT } from '../../../modules/local/samtools/rrna_count/main'

// rRNA quantification only. Replaces the standalone RSeQC modules (bam_stat /
// read_duplication / split_bam) — those are redundant with RustQC and the slow part.
// RustQC covers the rest of the post-alignment QC in a single fast pass.
workflow RRNA_COUNT {
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

        // Per-sample meta ([meta, bam, bai]) keyed to the canonical <group_id>_<sample_id>
        ch_rrna_inputs = ch_inputs_sorted.runnable
            .map { meta, tumor_bam, tumor_bai ->
                def meta_sample = Utils.getTumorRnaSample(meta)
                def meta_rrna = [
                    key: meta.group_id,
                    id: "${meta.group_id}_${meta_sample.sample_id}",
                    sample_id: Utils.getTumorRnaSampleName(meta),
                    subject_id: meta.subject_id,
                    group_id: meta.group_id,
                ]
                [meta_rrna, tumor_bam, tumor_bai]
            }

        // rRNA via samtools interval overlap against the rRNA BED (single fast pass)
        SAMTOOLS_RRNA_COUNT(ch_rrna_inputs, ch_bed)

        ch_rrna_stats   = WorkflowOncoanalyser.restoreMeta(SAMTOOLS_RRNA_COUNT.out.stats, ch_inputs)
        ch_rrna_multiqc = WorkflowOncoanalyser.restoreMeta(SAMTOOLS_RRNA_COUNT.out.multiqc, ch_inputs)

        // QC files for MultiQC (the rRNA custom-content YAML; not the stats txt)
        ch_qc_reports_final = channel.empty()
            .mix(ch_rrna_multiqc)
            .mix(ch_inputs_sorted.skip.map { meta -> [meta, []] })

    emit:
        qc_reports = ch_qc_reports_final
        rrna_stats = ch_rrna_stats  // [meta, stats_file]
        // Note: versions are collected via topics, not emitted here
}
