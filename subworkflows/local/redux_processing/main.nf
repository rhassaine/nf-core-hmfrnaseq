//
// REDUX processing: fixmate + duplicate marking/unmapping
// Runs samtools fixmate (MC tags) then REDUX on all BAM/CRAM inputs.
//

import Constants
import Utils
import WorkflowOncoanalyser

include { SAMTOOLS_FIXMATE_SORT } from '../../../modules/local/samtools/fixmate_sort/main'
include { REDUX                 } from '../../../modules/local/redux/main'

workflow REDUX_PROCESSING {
    take:
    // Sample data
    ch_inputs      // channel: [mandatory] [ meta ]
    ch_bam         // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta   // channel: [mandatory] /path/to/genome.fa
    genome_ver     // channel: [mandatory] genome version string (e.g., '38')
    genome_fai     // channel: [mandatory] /path/to/genome.fa.fai
    genome_dict    // channel: [mandatory] /path/to/genome.dict
    unmap_regions  // channel: [mandatory] /path/to/unmap_regions.tsv

    main:
    // channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = channel.empty()

    //
    // MODULE: SAMtools fixmate + sort
    // Sets MC (mate CIGAR) tag via streaming: name-sort → fixmate → coord-sort
    //
    // Build slim meta for process execution
    // channel: [ meta_process, bam ]
    ch_fixmate_inputs = ch_bam
        .map { meta, bam, bai ->
            def meta_sample = Utils.getTumorRnaSample(meta)
            def meta_process = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: Utils.getTumorRnaSampleName(meta),
                subject_id: meta.subject_id,
                group_id: meta.group_id,
                prefix: meta_sample.sample_id,
            ]
            return [meta_process, bam]
        }

    SAMTOOLS_FIXMATE_SORT(
        ch_fixmate_inputs,
    )

    ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE_SORT.out.versions)

    //
    // MODULE: REDUX duplicate marking
    //
    // channel: [ meta_process, [bam], [bai] ]
    ch_redux_inputs = SAMTOOLS_FIXMATE_SORT.out.bam
        .join(SAMTOOLS_FIXMATE_SORT.out.bai)
        .map { meta, bam, bai ->
            return [meta, [bam], bai instanceof List ? bai : [bai]]
        }

    // Run REDUX (MSI jitter skipped for RNA)
    REDUX(
        ch_redux_inputs,
        genome_fasta,
        genome_ver,
        genome_fai,
        genome_dict,
        unmap_regions,
        [],     // msi_jitter_sites: skip jitter analysis for RNA
        false,  // umi_enable
        '',     // umi_duplex_delim
    )

    ch_versions = ch_versions.mix(REDUX.out.versions)

    // Restore full meta from ch_inputs
    // REDUX.out.bam emits: [ meta, bam, bai ]
    // channel: [ meta, bam, bai ]
    ch_bams_ready = WorkflowOncoanalyser.groupByMeta(
        WorkflowOncoanalyser.restoreMeta(REDUX.out.bam.map { meta, bam, bai -> [meta, bam] }, ch_inputs),
        WorkflowOncoanalyser.restoreMeta(REDUX.out.bam.map { meta, bam, bai -> [meta, bai] }, ch_inputs),
    )

    emit:
    rna_tumor = ch_bams_ready  // channel: [ meta, bam, bai ]
    versions  = ch_versions    // channel: [ versions.yml ]
}
