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
include { RUSTQC_ANALYSIS       } from '../subworkflows/local/rustqc_analysis'

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

workflow RNA_STANDARD_WORKFLOW {
    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = channel.fromList(inputs)

    // The `rseqc` process stage is the master switch for post-alignment rRNA QC; the two counters are then
    // independently controlled with --skip_rustqc / --skip_rseqc (disable one, both, or neither).
    def run_rustqc = run_config.stages.rseqc && !params.skip_rustqc
    def run_rseqc  = run_config.stages.rseqc && !params.skip_rseqc

    // Validate reference requirements only for the counter(s) that will actually run
    if (run_rustqc && !params.ref_data_genome_gtf) {
        error "RustQC is enabled but --ref_data_genome_gtf is not set. Provide a GTF file or use --skip_rustqc (or --processes_exclude rseqc)"
    }
    if (run_rseqc && !params.rseqc_bed_file) {
        error "RSeQC is enabled but --rseqc_bed_file is not set. Provide an rRNA BED file or use --skip_rseqc (or --processes_exclude rseqc)"
    }

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    // Create channel for GTF input (RustQC biotype counting; staged via PREPARE_REFERENCE)
    // channel: [ meta2, gtf ]
    ch_gtf = ref_data.genome_gtf
        .map { gtf -> [ [ key: 'gtf', id: 'gtf' ], gtf ] }

    // Create channel for rRNA BED input (RSeQC split_bam; staged via PREPARE_REFERENCE)
    // channel: [ meta2, bed ]
    ch_bed = ref_data.rseqc_bed
        .map { bed -> [ [ key: 'bedfile', id: 'bedfile' ], bed ] }

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // TASK: FastQC on raw reads
    //
    ch_fastqc_out = channel.empty()
    if (run_config.stages.alignment) {
        // Build FastQC channel from independent channel.fromList to avoid
        // consuming ch_inputs (queue channels can only be read once per consumer)
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

        // Note: FastQC versions are collected via topics
        FASTQC(ch_fastq_for_qc)
        ch_fastqc_out = FASTQC.out.zip
    }

    ch_align_rna_tumor_out = channel.empty()
    ch_markdups_metrics = channel.empty()

    // Samples with pre-aligned input (BAM/CRAM) bypass alignment
    ch_prealigned = ch_inputs
        .filter { meta ->
            Utils.hasExistingInput(meta, Constants.INPUT.BAM_RNA_TUMOR) || Utils.hasExistingInput(meta, Constants.INPUT.CRAM_RNA_TUMOR)
        }
        .map { meta ->
            def sample = Utils.getTumorRnaSample(meta)
            meta.id = sample.sample_id ?: meta.subject_id ?: meta.group_id ?: 'unknown'
            def bam = Utils.getTumorRnaAlignment(meta) ?: []
            def bai = Utils.getTumorRnaAlignmentIndex(meta) ?: []
            [meta, bam, bai]
        }

    ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(ch_prealigned)

    if (run_config.stages.alignment) {

        // No rRNA pre-filtering in this workflow — feed raw FASTQs straight to STAR.
        // This builds the same [meta_fastq, fwd, rev] shape READ_ALIGNMENT_RNA expects
        // (previously supplied by SORTMERNA_FILTER.out.reads).
        ch_fastq_inputs = channel.fromList(inputs)
            .filter { meta ->
                Utils.hasTumorRnaFastq(meta) && !Utils.hasExistingInput(meta, Constants.INPUT.BAM_RNA_TUMOR)
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

        READ_ALIGNMENT_RNA(
            ch_inputs,
            ch_fastq_inputs,
            ref_data.genome_star_index,
        )

        ch_versions = ch_versions.mix(READ_ALIGNMENT_RNA.out.versions)

        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.rna_tumor)
        ch_markdups_metrics = ch_markdups_metrics.mix(READ_ALIGNMENT_RNA.out.markdups_metrics)

    }

    //
    // TASK: rRNA counting + quantification — NO GATE.
    //
    // The aligned-BAM channel fans out into three independent consumers so RustQC, RSeQC
    // and Isofox all run in parallel. Isofox is NOT filtered by rRNA content; every aligned
    // (and pre-aligned) sample is quantified.
    ch_align_rna_tumor_out
        .multiMap { meta, bam, bai ->
            rustqc: [meta, bam, bai]
            rseqc:  [meta, bam, bai]
            isofox: [meta, bam, bai]
        }
        .set { ch_bam_split }

    //
    // TASK: RustQC analysis (rRNA biotype counts + full QC)
    //
    ch_rustqc_out = channel.empty()
    if (run_rustqc) {
        // Note: versions are collected via topics
        RUSTQC_ANALYSIS(ch_inputs, ch_bam_split.rustqc, ch_gtf)
        ch_rustqc_out = RUSTQC_ANALYSIS.out.qc_reports
    } else {
        ch_rustqc_out = ch_inputs.map { meta -> [meta, []] }
    }

    //
    // TASK: RSeQC analysis (split_bam rRNA count + QC) — cross-check for RustQC
    //
    ch_rseqc_out = channel.empty()
    if (run_rseqc) {
        // Note: RSeQC versions are collected via topics
        RSEQC_ANALYSIS(ch_inputs, ch_bam_split.rseqc, ch_bed)
        ch_rseqc_out = RSEQC_ANALYSIS.out.qc_reports
    } else {
        ch_rseqc_out = ch_inputs.map { meta -> [meta, []] }
    }

    //
    // MODULE: Run Isofox to analyse RNA data (runs on all samples — no rRNA gate)
    //
    // channel: [ meta, isofox_dir ]
    ch_isofox_out = channel.empty()
    if (run_config.stages.isofox) {

        isofox_counts = params.isofox_counts ? file(params.isofox_counts) : hmf_data.isofox_counts
        isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios
        isofox_excluded_regions = params.isofox_excluded_regions ? file(params.isofox_excluded_regions) : hmf_data.isofox_excluded_regions

        ISOFOX_QUANTIFICATION(
            ch_inputs,
            ch_bam_split.isofox,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            isofox_counts,
            isofox_gc_ratios,
            [],  // isofox_gene_ids
            [],  // isofox_tpm_norm
            isofox_excluded_regions,
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
        // Load MultiQC config
        ch_multiqc_config = channel.fromPath(
            "$projectDir/assets/multiqc_config.yml", checkIfExists: true)

        // Group QC files by sample (group_id) for per-sample reports.
        // Both rRNA counters (whichever ran) are included so their numbers sit side by side.
        ch_multiqc_per_sample = channel.empty()
            .mix(ch_fastqc_out.map { meta, files -> [meta.key, files] })
            .mix(ch_rustqc_out.map { meta, files -> [meta.group_id ?: meta.key, files] })
            .mix(ch_rseqc_out.map { meta, files -> [meta.group_id ?: meta.key, files] })
            .mix(ch_markdups_metrics.map { meta, files -> [meta.group_id ?: meta.key, files] })
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
            ch_multiqc_config.toList(),
            [],
            [],
            [],
            []
        )

        // Aggregated MultiQC report (all samples, no FastQC - one row per sample)
        ch_multiqc_aggregated = channel.empty()
            .mix(ch_rustqc_out.map { meta, files -> files })
            .mix(ch_rseqc_out.map { meta, files -> files })
            .mix(ch_markdups_metrics.map { meta, files -> files })
            .flatten()
            .filter { it }
            .collect()
            .map { files ->
                def meta = [id: 'Aggregated', key: 'Aggregated']
                [meta, files]
            }

        MULTIQC_AGGREGATED(
            ch_multiqc_aggregated,
            ch_multiqc_config.toList(),
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
