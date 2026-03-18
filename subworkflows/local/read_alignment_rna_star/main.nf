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
    // MODULE: SAMtools coordinate sort
    //
    // Create process input channel
    // channel: [ meta_sort, bam ]
    ch_sort_inputs = STAR_ALIGN.out.bam
        .map { meta_star, bam ->
            def meta_sort = [
                *:meta_star,
                prefix: meta_star.read_group,
            ]

            return [meta_sort, bam]
        }

    // Run process
    SAMTOOLS_SORT(
        ch_sort_inputs,
    )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // MODULE: Sambamba merge
    //
    // Reunite BAMs
    // First, count expected BAMs per sample for non-blocking groupTuple op
    // channel: [ meta_count, group_size ]
    ch_sample_fastq_counts = ch_star_inputs
        .map { meta_star, reads_fwd, reads_rev ->
            def meta_count = [key: meta_star.key]
            return [meta_count, meta_star]
        }
        .groupTuple()
        .map { meta_count, meta_stars -> return [meta_count, meta_stars.size()] }

    // Now, group with expected size then sort into tumor and normal channels
    // channel: [ meta_group, [bam, ...] ]
    // NOTE: using combine(by:0) instead of cross — cross uses queue semantics that
    // can silently drop samples on -resume when cached results emit in a different
    // order or only partially. combine(by:0) does key-based matching reliably and
    // handles multi-lane samples (multiple BAMs per key) correctly via groupTuple.
    ch_bams_united = SAMTOOLS_SORT.out.bam
        .map { meta, bam -> [[key: meta.key], bam] }
        .combine(ch_sample_fastq_counts, by: 0)
        .map { meta_key, bam, group_size ->
            return tuple(groupKey(meta_key, group_size), bam)
        }
        .groupTuple()

    // Also group BAI files in the same way
    // channel: [ meta_group, [bai, ...] ]
    ch_bais_united = SAMTOOLS_SORT.out.bai
        .map { meta, bai -> [[key: meta.key], bai] }
        .combine(ch_sample_fastq_counts, by: 0)
        .map { meta_key, bai, group_size ->
            return tuple(groupKey(meta_key, group_size), bai)
        }
        .groupTuple()

    // Sort into merge-eligible BAMs (at least two BAMs required)
    // channel: runnable: [ meta_group, [bam, ...] ]
    // channel: skip: [ meta_group, bam ]
    ch_bams_united_sorted = ch_bams_united
        .branch { meta_group, bams ->
            runnable: bams.size() > 1
            skip: true
                return [meta_group, bams[0]]
        }

    // Also sort BAIs to match
    // channel: skip: [ meta_group, bai ]
    ch_bais_united_sorted = ch_bais_united
        .branch { meta_group, bais ->
            runnable: bais.size() > 1
            skip: true
                return [meta_group, bais[0]]
        }

    // Create process input channel
    // channel: [ meta_merge, [bams, ...] ]
    ch_merge_inputs = WorkflowOncoanalyser.restoreMeta(ch_bams_united_sorted.runnable, ch_inputs)
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
    // Collect merged BAMs (restored meta, ready for downstream REDUX_PROCESSING)
    //
    // channel: [ meta, bam, bai ]
    ch_bams_ready = channel.empty()
        .mix(
            // Merged BAMs - Sambamba merge doesn't produce BAI
            WorkflowOncoanalyser.restoreMeta(SAMBAMBA_MERGE.out.bam, ch_inputs)
                .map { meta, bam -> [meta, bam, []] },
            // Single BAMs already have BAI from SAMTOOLS_SORT
            WorkflowOncoanalyser.restoreMeta(ch_bams_united_sorted.skip, ch_inputs)
                .join(WorkflowOncoanalyser.restoreMeta(ch_bais_united_sorted.skip, ch_inputs)),
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
