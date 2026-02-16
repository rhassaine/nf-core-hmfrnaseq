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

include { ISOFOX_QUANTIFICATION } from '../subworkflows/local/isofox_quantification'
include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_RNA    } from '../subworkflows/local/read_alignment_rna'
include { RSEQC_ANALYSIS        } from '../subworkflows/local/rseqc_analysis'
include { SORTMERNA_FILTER      } from '../subworkflows/local/sortmerna_filter'

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

workflow RNA_WORKFLOW {
    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = channel.fromList(inputs)

    // ch_inputs.view()

    // Create channel for BED file input
    // channel: [ meta2, bed ]
    ch_bed = params.rseqc_bed_file
        ? channel.of([ [ key: 'bedfile', id: 'bedfile', ], file(params.rseqc_bed_file) ])
        : channel.empty()

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // TASK: FastQC on raw reads
    //
    ch_fastqc_out = channel.empty()
    if (run_config.stages.alignment) {
        // Create FASTQ input channel for FastQC
        ch_fastq_for_qc = ch_inputs
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

        // Note: FastQC versions are collected via topics
        FASTQC(ch_fastq_for_qc)
        ch_fastqc_out = FASTQC.out.zip
    }

    ch_align_rna_tumor_out = channel.empty()
    ch_markdups_metrics = channel.empty()
    ch_sortmerna_log = channel.empty()

    if (run_config.stages.alignment) {

        SORTMERNA_FILTER(
            ch_inputs,
            hmf_data.sortmerna_db,
        )

        READ_ALIGNMENT_RNA(
            ch_inputs,
            SORTMERNA_FILTER.out.reads,
            ref_data.genome_star_index,
        )

        ch_versions = ch_versions.mix(READ_ALIGNMENT_RNA.out.versions)

        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.rna_tumor)
        ch_markdups_metrics = ch_markdups_metrics.mix(READ_ALIGNMENT_RNA.out.markdups_metrics)

        ch_sortmerna_log = SORTMERNA_FILTER.out.log

    } else {

    ch_align_rna_tumor_out = ch_inputs.map { meta ->
        // enrich meta like alignment would do
        def sample = Utils.getTumorRnaSample(meta)
        meta.id = sample.sample_id ?: meta.subject_id ?: meta.group_id ?: 'unknown'
        [meta, [], []]
        }
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
    // TASK: rRNA QC gate - parse splitbam stats and branch samples
    //
    // Parse stats and add rRNA metrics to meta, then branch into pass/fail
    ch_rrna_qc_result = ch_splitbam_stats
        .map { meta, stats_file ->
            def rrna_stats = Utils.parseRrnaStats(stats_file)
            def qc_result = Utils.checkRrnaQc(
                rrna_stats,
                params.rrna_threshold_count ?: 0,
                params.rrna_threshold_percent ?: 0
            )
            // Log the rRNA stats for visibility
            log.info "rRNA QC [${meta.group_id}]: ${rrna_stats.rrna_reads}/${rrna_stats.total_reads} reads (${String.format('%.2f', rrna_stats.rrna_percent)}%) - ${qc_result.pass ? 'PASS' : 'FAIL'}"
            if (!qc_result.pass) {
                log.warn "Sample ${meta.group_id} FAILED rRNA QC: ${qc_result.fail_reason}"
            }
            [meta, qc_result.pass, rrna_stats]
        }
        .branch { meta, pass, rrna_stats ->
            pass: pass
                return meta
            fail: true
                return meta
        }

    //
    // TASK: SortMeRNA QC gate - parse logs and branch samples
    //
    // Aggregate SortMeRNA stats per sample (one log per lane) and branch pass/fail
    ch_sortmerna_qc_result = ch_sortmerna_log
        .map { meta, log_file ->
            def stats = Utils.parseSortmernaLog(log_file)
            [meta.key, stats.total_reads, stats.rrna_reads]
        }
        .groupTuple()
        .map { key, total_list, rrna_list ->
            def total = total_list.sum()
            def rrna = rrna_list.sum()
            def pct = total > 0 ? (rrna / (double)total) * 100.0 : 0.0
            [key, [total_reads: total, rrna_reads: rrna, rrna_percent: pct]]
        }

    ch_sortmerna_qc = ch_inputs
        .map { meta -> [meta.group_id, meta] }
        .join(ch_sortmerna_qc_result)
        .map { group_id, meta, rrna_stats ->
            def qc_result = Utils.checkRrnaQc(
                rrna_stats,
                params.sortmerna_threshold_count ?: 0,
                params.sortmerna_threshold_percent ?: 0
            )
            log.info "SortMeRNA QC [${meta.group_id}]: ${rrna_stats.rrna_reads}/${rrna_stats.total_reads} reads (${String.format('%.2f', rrna_stats.rrna_percent)}%) - ${qc_result.pass ? 'PASS' : 'FAIL'}"
            if (!qc_result.pass) {
                log.warn "Sample ${meta.group_id} FAILED SortMeRNA QC: ${qc_result.fail_reason}"
            }
            [meta, qc_result.pass, rrna_stats]
        }
        .branch { meta, pass, rrna_stats ->
            pass: pass
                return meta
            fail: true
                return meta
        }

    // Samples must pass both SortMeRNA gate (if alignment ran) AND RSeQC gate (if RSeQC ran)
    ch_samples_for_isofox = ch_bam_split.isofox

    if (run_config.stages.alignment) {
        ch_samples_for_isofox = ch_samples_for_isofox
            .map { meta, bam, bai -> [meta.group_id, meta, bam, bai] }
            .join(ch_sortmerna_qc.pass.map { meta -> [meta.group_id, meta] }, by: 0)
            .map { group_id, meta_bam, bam, bai, meta_qc -> [meta_bam, bam, bai] }
    }

    if (run_config.stages.rseqc) {
        ch_samples_for_isofox = ch_samples_for_isofox
            .map { meta, bam, bai -> [meta.group_id, meta, bam, bai] }
            .join(ch_rrna_qc_result.pass.map { meta -> [meta.group_id, meta] }, by: 0)
            .map { group_id, meta_bam, bam, bai, meta_qc -> [meta_bam, bam, bai] }
    }

    //
    // MODULE: Run Isofox to analyse RNA data
    //
    // channel: [ meta, isofox_dir ]
    ch_isofox_out = channel.empty()
    if (run_config.stages.isofox) {

        isofox_counts = params.isofox_counts ? file(params.isofox_counts) : hmf_data.isofox_counts
        isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios

        // Get inputs for samples that passed rRNA QC (both gates)
        ch_inputs_for_isofox = ch_inputs

        if (run_config.stages.alignment) {
            ch_inputs_for_isofox = ch_inputs_for_isofox
                .map { meta -> [meta.group_id, meta] }
                .join(ch_sortmerna_qc.pass.map { meta -> [meta.group_id] }, by: 0)
                .map { group_id, meta -> meta }
        }

        if (run_config.stages.rseqc) {
            ch_inputs_for_isofox = ch_inputs_for_isofox
                .map { meta -> [meta.group_id, meta] }
                .join(ch_rrna_qc_result.pass.map { meta -> [meta.group_id] }, by: 0)
                .map { group_id, meta -> meta }
        }

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
    if (run_config.stages.multiqc) {
        // Group QC files by sample (group_id) for per-sample reports
        ch_multiqc_per_sample = channel.empty()
            .mix(ch_fastqc_out.map { meta, files -> [meta.key, files] })
            .mix(ch_rseqc_out.map { meta, files -> [meta.group_id ?: meta.key, files] })
            .mix(ch_markdups_metrics.map { meta, files -> [meta.group_id ?: meta.key, files] })
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
            .mix(ch_markdups_metrics.map { meta, files -> files })
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
