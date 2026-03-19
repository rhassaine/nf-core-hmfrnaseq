//
// Align RNA reads with STAR (coordinate sort + merge)
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
    // Sample data
    ch_inputs         // channel: [mandatory] [ meta ]
    ch_fastq_inputs   // channel: [mandatory] [ meta_fastq, fwd, rev ] (from SORTMERNA_FILTER)

    // Reference data
    genome_star_index // channel: [mandatory] /path/to/genome_star_index/

    main:
    // channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = channel.empty()

    // Sort inputs
    // channel: [ meta ]
    ch_inputs_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_RNA_TUMOR)
            runnable: Utils.hasTumorRnaFastq(meta) && !has_existing
            skip: true
        }

    //
    // MODULE: STAR alignment
    //
    // Create process input channel
    // channel: [ meta_star, fastq_fwd, fastq_rev ]
    ch_star_inputs = ch_fastq_inputs
        .map { meta_fastq, fastq_fwd, fastq_rev ->
            def meta_star = [
                *:meta_fastq,
                read_group: "${meta_fastq.sample_id}.${meta_fastq.library_id}.${meta_fastq.lane}",
            ]

            return [meta_star, fastq_fwd, fastq_rev]
        }

    // Run process
    STAR_ALIGN(
        ch_star_inputs,
        genome_star_index,
    )

    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    //
    // Route BAMs: multi-lane samples need coord-sort + merge, single-lane skip both
    //
    // Count expected BAMs per sample for grouping
    // channel: [ [key: X], group_size ]
    ch_sample_fastq_counts = ch_star_inputs
        .map { meta_star, reads_fwd, reads_rev ->
            return [[key: meta_star.key], meta_star]
        }
        .groupTuple()
        .map { meta_count, meta_stars -> return [meta_count, meta_stars.size()] }

    // Tag each STAR BAM with its sample's lane count, then branch
    // channel: multi_lane:  [ meta_star, bam ]
    // channel: single_lane: [ meta_star, bam ]
    ch_star_bams_routed = STAR_ALIGN.out.bam
        .map { meta, bam -> [[key: meta.key], meta, bam] }
        .combine(ch_sample_fastq_counts, by: 0)
        .map { meta_key, meta_star, bam, group_size -> [meta_star, bam, group_size] }
        .branch { meta_star, bam, group_size ->
            multi_lane: group_size > 1
                return [meta_star, bam]
            single_lane: true
                return [meta_star, bam]
        }

    //
    // MODULE: SAMtools coordinate sort (multi-lane only — needed for SAMBAMBA_MERGE)
    //
    // channel: [ meta_sort, bam ]
    ch_sort_inputs = ch_star_bams_routed.multi_lane
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
    // MODULE: Sambamba merge (multi-lane only)
    //
    // Group sorted BAMs by sample key for merging
    // channel: [ meta_group, [bam, ...] ]
    ch_bams_for_merge = SAMTOOLS_SORT.out.bam
        .map { meta, bam -> [[key: meta.key], bam] }
        .combine(ch_sample_fastq_counts, by: 0)
        .map { meta_key, bam, group_size ->
            return tuple(groupKey(meta_key, group_size), bam)
        }
        .groupTuple()

    // Create process input channel
    // channel: [ meta_merge, [bams, ...] ]
    ch_merge_inputs = WorkflowOncoanalyser.restoreMeta(ch_bams_for_merge, ch_inputs)
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

    // Run process
    SAMBAMBA_MERGE(
        ch_merge_inputs,
    )

    ch_versions = ch_versions.mix(SAMBAMBA_MERGE.out.versions)

    //
    // Collect BAMs (restored meta, ready for downstream REDUX_PROCESSING)
    //
    // channel: [ meta, bam, bai ]
    ch_bams_ready = channel.empty()
        .mix(
            // Multi-lane: merged BAMs (Sambamba merge doesn't produce BAI)
            WorkflowOncoanalyser.restoreMeta(SAMBAMBA_MERGE.out.bam, ch_inputs)
                .map { meta, bam -> [meta, bam, []] },
            // Single-lane: unsorted STAR BAMs (no BAI needed, FIXMATE_SORT re-sorts)
            WorkflowOncoanalyser.restoreMeta(ch_star_bams_routed.single_lane, ch_inputs)
                .map { meta, bam -> [meta, bam, []] },
        )

    // Set outputs
    // channel: [ meta, bam, bai ]
    ch_bam_out = channel.empty()
        .mix(
            ch_bams_ready,
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    rna_tumor = ch_bam_out   // channel: [ meta, bam, bai ]
    versions  = ch_versions  // channel: [ versions.yml ]
}
