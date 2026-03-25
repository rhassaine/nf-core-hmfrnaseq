import Constants
import Processes
import Utils

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse input samplesheet
// NOTE(SW): this is done early and outside of gpars so that we can access synchronously and prior to pipeline execution
inputs = Utils.parseInput(params.input, workflow.stubRun, log)

// Get run config
run_config = WorkflowMain.getRunConfig(params, inputs, log)

// Validate inputs
Utils.validateInput(inputs, run_config, params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.isofox_counts,
    params.isofox_gc_ratios,
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Used in Isofox subworkflow
isofox_read_length = params.isofox_read_length !== null ? params.isofox_read_length : Constants.DEFAULT_ISOFOX_READ_LENGTH

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { ISOFOX_QUANTIFICATION      } from '../subworkflows/local/isofox_quantification'
include { PREPARE_REFERENCE          } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_RNA_STAR    } from '../subworkflows/local/read_alignment_rna_star'
include { REDUX_PROCESSING           } from '../subworkflows/local/redux_processing'
include { RRNA_QC_GATE              } from '../subworkflows/local/rrna_qc_gate'
include { RSEQC_ANALYSIS             } from '../subworkflows/local/rseqc_analysis'
include { SORTMERNA_FILTER           } from '../subworkflows/local/sortmerna_filter'

include { AMBER                              } from '../modules/local/amber/main'
include { MULTIQC                            } from '../modules/local/multiqc/main'
include { MULTIQC as MULTIQC_AGGREGATED     } from '../modules/local/multiqc/main'
include { FASTQC                } from '../modules/nf-core/fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
samplesheet = Utils.getFileObject(params.input)

workflow RNA_REDUX_WORKFLOW {
    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = channel.fromList(inputs)

    // Validate BED file requirement for RSeQC
    if (run_config.stages.rseqc && !params.rseqc_bed_file) {
        error "RSeQC is enabled but --rseqc_bed_file is not set. Provide an rRNA BED file or use --processes_exclude rseqc"
    }

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    // Create channel for BED file input (staged via PREPARE_REFERENCE)
    // channel: [ meta2, bed ]
    ch_bed = ref_data.rseqc_bed
        .map { bed -> [ [ key: 'bedfile', id: 'bedfile' ], bed ] }

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // TASK: FastQC on raw reads
    //
    ch_fastqc_out = channel.empty()
    if (run_config.stages.alignment) {
        // Build FastQC channel from its own independent channel.fromList
        ch_fastq_for_qc = channel.fromList(inputs)
            .filter { meta -> Utils.hasTumorRnaFastq(meta) }
            .flatMap { meta ->
                def meta_sample = Utils.getTumorRnaSample(meta)
                meta_sample
                    .getAt(Constants.FileType.FASTQ)
                    .collect { key, fps ->
                        def (library_id, lane) = key
                        def meta_fastqc = [
                            key: meta.group_id,
                            id: "${meta.group_id}_${meta_sample.sample_id}_${library_id}_${lane}",
                            sample_id: meta_sample.sample_id,
                            library_id: library_id,
                            lane: lane,
                        ]
                        return [meta_fastqc, [fps['fwd'], fps['rev']]]
                    }
            }

        FASTQC(ch_fastq_for_qc)
        ch_fastqc_out = FASTQC.out.zip
    }

    //
    // TASK: Alignment (FASTQ → STAR → sort → merge)
    //
    ch_aligned = channel.empty()
    ch_sortmerna_log = channel.empty()

    if (run_config.stages.alignment) {

        // Build alignment channel from its own independent channel.fromList
        ch_raw_fastq_inputs = channel.fromList(inputs)
            .filter { meta ->
                Utils.hasTumorRnaFastq(meta) && !Utils.hasExistingInput(meta, Constants.INPUT.BAM_RNA_TUMOR) && !Utils.hasExistingInput(meta, Constants.INPUT.CRAM_RNA_TUMOR)
            }
            .flatMap { meta ->
                def meta_sample = Utils.getTumorRnaSample(meta)
                meta_sample
                    .getAt(Constants.FileType.FASTQ)
                    .collect { key, fps ->
                        def (library_id, lane) = key
                        def meta_fastq = [
                            key: meta.group_id,
                            id: "${meta.group_id}_${meta_sample.sample_id}",
                            sample_id: meta_sample.sample_id,
                            library_id: library_id,
                            lane: lane,
                        ]
                        [meta_fastq, fps['fwd'], fps['rev']]
                    }
            }

        // Optionally filter rRNA reads with SortMeRNA before alignment
        ch_reads_for_alignment = ch_raw_fastq_inputs
        if (!params.skip_sortmerna) {
            SORTMERNA_FILTER(
                ch_inputs,
                ref_data.sortmerna_db,
            )
            ch_reads_for_alignment = SORTMERNA_FILTER.out.reads
            ch_sortmerna_log = SORTMERNA_FILTER.out.sort_log
        }

        READ_ALIGNMENT_RNA_STAR(
            ch_inputs,
            ch_reads_for_alignment,
            ref_data.genome_star_index,
        )

        ch_versions = ch_versions.mix(READ_ALIGNMENT_RNA_STAR.out.versions)

        // Filter out skip entries (samples without FASTQ input)
        ch_aligned = READ_ALIGNMENT_RNA_STAR.out.rna_tumor
            .filter { meta, bam, bai -> bam }

    }

    //
    // TASK: Pre-aligned BAM/CRAM samples (bypass alignment)
    //
    // channel: [ meta, bam, bai ]
    ch_prealigned = ch_inputs
        .filter { meta ->
            Utils.hasExistingInput(meta, Constants.INPUT.BAM_RNA_TUMOR) || Utils.hasExistingInput(meta, Constants.INPUT.CRAM_RNA_TUMOR)
        }
        .map { meta ->
            def bam = Utils.getTumorRnaAlignment(meta) ?: []
            def bai = Utils.getTumorRnaAlignmentIndex(meta) ?: []
            [meta, bam, bai]
        }

    //
    // TASK: REDUX processing (fixmate + duplicate marking/unmapping for all BAMs)
    //
    REDUX_PROCESSING(
        ch_inputs,
        ch_aligned.mix(ch_prealigned),
        ref_data.genome_fasta,
        ref_data.genome_version,
        ref_data.genome_fai,
        ref_data.genome_dict,
        hmf_data.unmap_regions,
    )

    ch_versions = ch_versions.mix(REDUX_PROCESSING.out.versions)

    ch_align_rna_tumor_out = REDUX_PROCESSING.out.rna_tumor

    //
    // TASK: AMBER BAF profiling (tumor-only, RNA BAMs)
    //
    if (run_config.stages.amber) {
        ch_amber_inputs = ch_align_rna_tumor_out
            .map { meta, bam, bai ->
                def meta_sample = Utils.getTumorRnaSample(meta)
                def meta_amber = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: meta_sample.sample_id,
                ]
                [meta_amber, bam, bai]
            }

        AMBER(
            ch_amber_inputs,
            ref_data.genome_version,
            hmf_data.heterozygous_sites,
        )

        ch_versions = ch_versions.mix(AMBER.out.versions)
    }

    //
    // TASK: RSeQC QC analysis (must run before Isofox for rRNA contamination check)
    //
    // Split channel for multiple consumers
    ch_align_rna_tumor_out
        .multiMap { meta, bam, bai ->
            rseqc: [meta, bam, bai]
            isofox: [meta, bam, bai]
        }
        .set { ch_bam_split }

    ch_rseqc_out = channel.empty()
    ch_splitbam_stats = channel.empty()
    if (run_config.stages.rseqc) {
        // Run RSeQC QC on aligned BAMs
        // Note: RSeQC versions are collected via topics
        RSEQC_ANALYSIS(ch_inputs, ch_bam_split.rseqc, ch_bed)

        ch_rseqc_out = RSEQC_ANALYSIS.out.qc_reports
        ch_splitbam_stats = RSEQC_ANALYSIS.out.splitbam_stats

    } else {

        ch_rseqc_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // TASK: rRNA QC gate - parse splitbam stats and gate samples for Isofox
    //
    ch_samples_for_isofox = ch_bam_split.isofox
    ch_inputs_for_isofox = ch_inputs

    if (run_config.stages.rseqc) {
        RRNA_QC_GATE(
            ch_inputs,
            ch_splitbam_stats,
            ch_sortmerna_log,
            ch_bam_split.isofox,
            params.rrna_threshold_count,
            params.rrna_threshold_percent,
        )

        ch_samples_for_isofox = RRNA_QC_GATE.out.bam_pass
        ch_inputs_for_isofox = RRNA_QC_GATE.out.inputs_pass
    }

    //
    // MODULE: Run Isofox to analyse RNA data
    //
    // channel: [ meta, isofox_dir ]
    ch_isofox_out = channel.empty()
    if (run_config.stages.isofox) {

        isofox_counts = params.isofox_counts ? file(params.isofox_counts) : hmf_data.isofox_counts
        isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios

        ISOFOX_QUANTIFICATION(
            ch_inputs_for_isofox,
            ch_samples_for_isofox,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            isofox_counts,
            isofox_gc_ratios,
            [],  // isofox_gene_ids
            [],  // isofox_tpm_norm
            params.isofox_functions,
            isofox_read_length,
        )

        ch_versions = ch_versions.mix(ISOFOX_QUANTIFICATION.out.versions)

        ch_isofox_out = ch_isofox_out.mix(ISOFOX_QUANTIFICATION.out.isofox_dir)

    } else {

        ch_isofox_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // TASK: MultiQC (per-sample reports)
    //
    // Note: REDUX metrics are not MultiQC-compatible, so only FastQC and RSeQC are included
    if (run_config.stages.multiqc) {
        // Group QC files by sample (group_id) for per-sample reports
        ch_multiqc_per_sample = channel.empty()
            .mix(ch_fastqc_out.map { meta, files -> [meta.key, files] })
            .mix(ch_rseqc_out.map { meta, files -> [meta.group_id ?: meta.key, files] })
            .mix(ch_sortmerna_log.map { meta, log_file -> [meta.key, log_file] })
            .filter { group_id, files -> files }
            .groupTuple(by: 0)
            .map { group_id, file_lists ->
                def all_files = file_lists.flatten().findAll { it }
                def meta = [id: group_id, key: group_id]
                [meta, all_files]
            }
            .filter { meta, files -> files.size() > 0 }

        MULTIQC(
            ch_multiqc_per_sample,
            [],
            [],
            [],
            [],
            []
        )

        // Aggregated MultiQC report (all samples, no FastQC - one row per sample)
        ch_multiqc_aggregated = channel.empty()
            .mix(ch_rseqc_out.map { meta, files -> files })
            .mix(ch_sortmerna_log.map { meta, log_file -> log_file })
            .flatten()
            .filter { it }
            .collect()
            .map { files ->
                def meta = [id: 'aggregated', key: 'aggregated']
                [meta, files]
            }

        MULTIQC_AGGREGATED(
            ch_multiqc_aggregated,
            [],
            [],
            [],
            [],
            []
        )
    }

    // Note: MultiQC versions not collected - module emits tuples, not YAML files

    //
    // TASK: Aggregate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true,
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
