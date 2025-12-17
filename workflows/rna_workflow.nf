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

include { MULTIQC               } from '../modules/nf-core/multiqc/main'
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
    ch_bed = channel.of([ [ key: 'bedfile', id: 'bedfile', ], file(params.rseqc_bed_file) ])

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

    if (run_config.stages.alignment) {

        READ_ALIGNMENT_RNA(
            ch_inputs,
            ref_data.genome_star_index,
        )

        ch_versions = ch_versions.mix(READ_ALIGNMENT_RNA.out.versions)

        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.rna_tumor)

    } else {

    ch_align_rna_tumor_out = ch_inputs.map { meta ->
        // enrich meta like alignment would do
        def sample = Utils.getTumorRnaSample(meta)
        meta.id = sample.sample_id ?: meta.subject_id ?: meta.group_id ?: 'unknown'
        [meta, [], []]
        }
    }

    ch_align_rna_tumor_out.view { it[0].id }

    //
    // MODULE: Run Isofox to analyse RNA data
    //
   // channel: [ meta, isofox_dir ]
    ch_isofox_out = channel.empty()
    if (run_config.stages.isofox) {

        isofox_counts = params.isofox_counts ? file(params.isofox_counts) : hmf_data.isofox_counts
        isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios

        ISOFOX_QUANTIFICATION(
            ch_inputs,
            ch_align_rna_tumor_out,
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


    ch_rseqc_out = channel.empty()
    if (run_config.stages.rseqc) {
        // Run RSeQC QC on aligned BAMs
        // Note: RSeQC versions are collected via topics
        RSEQC_ANALYSIS(ch_inputs, ch_align_rna_tumor_out, ch_bed)

        ch_rseqc_out = RSEQC_ANALYSIS.out.qc_reports

    } else {

        ch_rseqc_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // TASK: MultiQC
    //
    ch_multiqc_files = channel.empty()
        .mix(ch_fastqc_out.map { meta, files -> files }.flatten().filter { it })
        .mix(ch_rseqc_out.map { meta, files -> files }.flatten().filter { it })
        .collect()

    MULTIQC(
        ch_multiqc_files,
        [],
        [],
        [],
        [],
        []
    )

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
