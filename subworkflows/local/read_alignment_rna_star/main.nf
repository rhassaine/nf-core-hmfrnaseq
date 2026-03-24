//
// Align RNA reads with STAR → SAMtools sort → Sambamba merge
// Downstream REDUX_PROCESSING handles fixmate and duplicate marking.
//

import Constants
import Utils
import WorkflowOncoanalyser

include { SAMBAMBA_MERGE  } from '../../../modules/local/sambamba/merge/main'
include { SAMTOOLS_SORT   } from '../../../modules/local/samtools/sort/main'
include { STAR_ALIGN      } from '../../../modules/local/star/align/main'

workflow READ_ALIGNMENT_RNA_STAR {
    take:
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_fastq_inputs     // channel: [mandatory] [ meta_fastq, fwd, rev ]
    genome_star_index   // channel: [mandatory] /path/to/genome_star_index/

    main:
    ch_versions = channel.empty()

    // Branch: samples with existing BAMs skip alignment
    ch_inputs_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_RNA_TUMOR)
            runnable: Utils.hasTumorRnaFastq(meta) && !has_existing
            skip: true
        }

    //
    // MODULE: STAR alignment (per lane)
    //
    ch_star_inputs = ch_fastq_inputs
        .map { meta_fastq, fastq_fwd, fastq_rev ->
            def meta_star = [
                *:meta_fastq,
                read_group: "${meta_fastq.sample_id}.${meta_fastq.library_id}.${meta_fastq.lane}",
            ]
            return [meta_star, fastq_fwd, fastq_rev]
        }

    STAR_ALIGN(
        ch_star_inputs,
        genome_star_index,
    )

    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    //
    // MODULE: SAMtools coordinate sort (all BAMs)
    //
    ch_sort_inputs = STAR_ALIGN.out.bam
        .map { meta_star, bam ->
            def meta_sort = [
                *:meta_star,
                prefix: meta_star.read_group,
            ]
            return [meta_sort, bam]
        }

    SAMTOOLS_SORT(
        ch_sort_inputs,
    )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // Group sorted BAMs by sample key for merge decision
    //
    ch_bams_grouped = SAMTOOLS_SORT.out.bam
        .map { meta, bam -> [[key: meta.key], bam] }
        .groupTuple()

    ch_bais_grouped = SAMTOOLS_SORT.out.bai
        .map { meta, bai -> [[key: meta.key], bai] }
        .groupTuple()

    // Branch: multi-lane needs merge, single-lane passes through
    ch_bams_branched = ch_bams_grouped
        .branch { meta_group, bams ->
            runnable: bams.size() > 1
            skip: true
                return [meta_group, bams[0]]
        }

    ch_bais_branched = ch_bais_grouped
        .branch { meta_group, bais ->
            runnable: bais.size() > 1
            skip: true
                return [meta_group, bais[0]]
        }

    //
    // MODULE: Sambamba merge (multi-lane only)
    //
    ch_merge_inputs = WorkflowOncoanalyser.restoreMeta(ch_bams_branched.runnable, ch_inputs)
        .map { meta, bams ->
            def meta_sample = Utils.getTumorRnaSample(meta)
            def meta_merge = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: Utils.getTumorRnaSampleName(meta),
                subject_id: meta.subject_id,
                group_id: meta.group_id,
            ]
            return [meta_merge, bams]
        }

    SAMBAMBA_MERGE(
        ch_merge_inputs,
    )

    ch_versions = ch_versions.mix(SAMBAMBA_MERGE.out.versions)

    //
    // Collect BAMs ready for downstream
    //
    ch_bams_ready = channel.empty()
        .mix(
            // Multi-lane: merged BAMs
            WorkflowOncoanalyser.restoreMeta(SAMBAMBA_MERGE.out.bam, ch_inputs)
                .map { meta, bam -> [meta, bam, []] },
            // Single-lane: sorted BAM + BAI
            WorkflowOncoanalyser.restoreMeta(ch_bams_branched.skip, ch_inputs)
                .join(WorkflowOncoanalyser.restoreMeta(ch_bais_branched.skip, ch_inputs)),
        )

    // Set outputs
    ch_bam_out = channel.empty()
        .mix(
            ch_bams_ready,
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    rna_tumor = ch_bam_out   // channel: [ meta, bam, bai ]
    versions  = ch_versions  // channel: [ versions.yml ]
}
