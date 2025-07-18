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
        ch_inputs         // channel: [mandatory] [ meta ]
        ch_tumor_rna_bam  // channel: [meta, bam, bai]

    main:

        // ch_bam_for_rseqc = ch_tumor_rna_bam.map { meta, bam, bai ->
        //     def sample_id = meta?.id ?: meta?.sample_id ?: meta?.subject_id ?: meta?.group_id ?: 'unknown'
        //     def key_val = meta?.key ?: sample_id
        //     meta.key = key_val
        //     meta.id = sample_id
        //     [meta, bam]
        // }

        // ch_bam_for_rseqc = ch_tumor_rna_bam.map { meta, bam, bai ->
        //     def sample_id = meta?.sample_id ?: meta?.id ?: meta?.subject_id ?: meta?.group_id ?: 'unknown'
        //     def group_id = meta?.group_id ?: sample_id
        //     def enriched_meta = [
        //         key: group_id,
        //         id: sample_id,
        //         sample_id: sample_id,
        //         group_id: group_id,
        //         lane: meta?.lane,
        //         library_id: meta?.library_id,
        // // add any other fields you want to propagate
        //     ]
        //     [enriched_meta, bam]
        // }

        // // Run modules
        // RSEQC_BAMSTAT(ch_bam_for_rseqc)
        // RSEQC_READDUPLICATION(ch_bam_for_rseqc)

        // // Collect version info
        // ch_versions = Channel.empty()
        //     .mix(RSEQC_BAMSTAT.out.versions, RSEQC_READDUPLICATION.out.versions)

        // // Collect QC outputs for MultiQC
        // ch_qc_reports = Channel.empty()
        //     .mix(
        //         RSEQC_BAMSTAT.out.bamstat,
        //         RSEQC_READDUPLICATION.out.seq_xls,
        //         RSEQC_READDUPLICATION.out.pos_xls,
        //         RSEQC_READDUPLICATION.out.pdf,
        //         RSEQC_READDUPLICATION.out.rscript
        //     )
        // Downstream, always treat ch_qc_reports as tuple (meta, file)


        // ____________________________________________________________________________________
        /////////// ATTEMPTS AT FIXING METADATA ISSUE 
        // ____________________________________________________________________________________

        // Initialize version channel
        ch_versions = Channel.empty()

        // Select input sources and sort
        // channel: runnable: [ meta, tumor_bam, tumor_bai ]
        // channel: skip:     [ meta ]
        ch_inputs_sorted = ch_tumor_rna_bam
            .map { meta, tumor_bam, tumor_bai ->
                println "[DEBUG] meta: ${meta}"
                def tumor_rna_sample = Utils.getTumorRnaSample(meta)
                println "[DEBUG] tumor_rna_sample: ${tumor_rna_sample}"
                def sample_name = tumor_rna_sample['sample_id']
                println "[DEBUG] sample_name: ${sample_name}"
                return [
                    meta,
                    Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                    Utils.selectCurrentOrExisting(tumor_bai, meta, Constants.INPUT.BAI_RNA_TUMOR),
                ]
            }
            .branch { meta, tumor_bam, tumor_bai ->
                // Add debug here too
                println "[BRANCH DEBUG] meta: ${meta == null ? 'NULL' : meta}"
                println "[BRANCH DEBUG] tumor_bam: ${tumor_bam == null ? 'NULL' : tumor_bam}"
                println "[BRANCH DEBUG] tumor_bai: ${tumor_bai == null ? 'NULL' : tumor_bai}"
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.RSEQ_DIR)
                runnable: tumor_bam && !has_existing
                skip:     true
            }

        // Inspect both branches
        ch_inputs_sorted.runnable.view { println it }           // see runnable tuples
        ch_inputs_sorted.runnable.map { it.size() }.unique().view()
        ch_inputs_sorted.skip.view { println it }               // see skipped metas

        // Create process input channel
        // channel: [ meta_rseqc, tumor_bam, tumor_bai ]
        ch_rseqc_inputs = ch_inputs_sorted.runnable
            .map { meta, tumor_bam, tumor_bai ->

                println "[RSEQ INPUT DEBUG] meta: ${meta}"
                println "[RSEQ INPUT DEBUG] tumor_bam: ${tumor_bam}"
                println "[RSEQ INPUT DEBUG] tumor_bai: ${tumor_bai}"
                def sample_name = Utils.getTumorRnaSampleName(meta)
                println "[RSEQ INPUT DEBUG] sample_name: ${sample_name}"
                def meta_rseqc = [
                    key:       meta.group_id,
                    id:        meta.group_id,
                    sample_id: sample_name,
                ]

                return [meta_rseqc, tumor_bam, tumor_bai]
            }
        
        // Run RSEQC modules
        RSEQC_BAMSTAT        (ch_rseqc_inputs)
        ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions)

        RSEQC_READDUPLICATION(ch_rseqc_inputs)
        ch_versions = ch_versions.mix(RSEQC_READDUPLICATION.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, rseqc_dir ]
        ch_outputs = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(RSEQC_BAMSTAT.out.bamstat, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] }
            )

        // Collect QC outputs for MultiQC
        ch_qc_reports = Channel.empty()
             .mix(
                RSEQC_BAMSTAT.out.bamstat,
                RSEQC_READDUPLICATION.out.seq_xls,
                RSEQC_READDUPLICATION.out.pos_xls,
                RSEQC_READDUPLICATION.out.pdf,
                RSEQC_READDUPLICATION.out.rscript
            )

    emit:
        // bamstat    = RSEQC_BAMSTAT.out.bamstat
        // seq_xls    = RSEQC_READDUPLICATION.out.seq_xls
        // pos_xls    = RSEQC_READDUPLICATION.out.pos_xls
        // pdf        = RSEQC_READDUPLICATION.out.pdf
        // rscript    = RSEQC_READDUPLICATION.out.rscript
        // versions   = ch_versions
        qc_reports = ch_qc_reports

        rseqc_dir  = ch_outputs

        versions   = ch_versions
}