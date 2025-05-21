import Constants
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

// TODO(SW): consider whether we should check for null entries here for errors to be more informative
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Used in Isofox and Neo subworkflows
isofox_read_length = params.isofox_read_length !== null ? params.isofox_read_length : Constants.DEFAULT_ISOFOX_READ_LENGTH_WTS

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { ISOFOX_QUANTIFICATION } from '../subworkflows/local/isofox_quantification'
include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_RNA    } from '../subworkflows/local/read_alignment_rna'

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
    ch_versions = Channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = Channel.fromList(inputs)

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //
    READ_ALIGNMENT_RNA(
        ch_inputs,
        ref_data.genome_star_index,
        ref_data.genome_fasta,
        ref_data.genome_fai,
        params.max_fastq_records,
        params.fastp_umi,
        params.fastp_umi_location,
        params.fastp_umi_length,
        params.fastp_umi_skip
    )

    ch_versions = ch_versions.mix(READ_ALIGNMENT_RNA.out.versions)

    //
    // MODULE: Run Isofox to analyse RNA data
    //
    isofox_counts = params.isofox_counts ? file(params.isofox_counts) : hmf_data.isofox_counts
    isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios

    ISOFOX_QUANTIFICATION(
        ch_inputs,
        READ_ALIGNMENT_RNA.out.rna_tumor,
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