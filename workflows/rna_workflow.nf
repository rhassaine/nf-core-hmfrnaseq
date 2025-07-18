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

if (run_config.stages.lilac) {
    if (params.genome_version.toString() == '38' && params.genome_type == 'alt' && params.containsKey('ref_data_hla_slice_bed')) {
        checkPathParamList.add(params.ref_data_hla_slice_bed)
    }
}

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
include { RSEQC_ANALYSIS        } from '../subworkflows/local/rseqc_analysis'

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

    // ch_inputs.view()

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    // Set GRIDSS config
    gridss_config = params.gridss_config !== null ? file(params.gridss_config) : hmf_data.gridss_config

    ch_align_rna_tumor_out = Channel.empty()

    if (run_config.stages.alignment) {
        
        READ_ALIGNMENT_RNA(
            ch_inputs,
            ref_data.genome_star_index,
        )
        
        ch_versions = ch_versions.mix(READ_ALIGNMENT_RNA.out.versions)

        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.rna_tumor)
    
    // } else {

    //     ch_align_rna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
    // }
    } else {
        
    ch_align_rna_tumor_out = ch_inputs.map { meta ->
        // enrich meta like alignment would do
        def sample = Utils.getTumorRnaSample(meta)
        meta.id = sample.sample_id ?: meta.subject_id ?: meta.group_id ?: 'unknown'
        [meta, [], []]
        }
    }

    ch_align_rna_tumor_out.view { it[0].id }

    // ch_align_rna_tumor_out.view()
    // println "ISOFOX stage enabled: ${run_config.stages.isofox}"

    //
    // MODULE: Run Isofox to analyse RNA data
    //
   // channel: [ meta, isofox_dir ]
    ch_isofox_out = Channel.empty()
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
    

    ch_rsqec_out = Channel.empty()
    if (run_config.stages.rseqc) {
        // Run RSeQC QC on aligned BAMs
        RSEQC_ANALYSIS(ch_inputs, ch_align_rna_tumor_out)

        ch_versions = ch_versions.mix(RSEQC_ANALYSIS.out.versions)

    // ch_rseqc_out = ch_rseqc_out.mix(RSEQC_ANALYSIS.out.txt)
    // ch_rseqc_out = ch_rseqc_out.mix(RSEQC_ANALYSIS.out.bamstat)
    
        ch_rseqc_out = RSEQC_ANALYSIS.out.qc_reports


    } else {
    
        ch_rseqc_out_bamstat_out = ch_inputs.map { meta -> [meta, [], []] }
    
    }

    // Additional development logging 

    // ch_isofox_out.view()
    // ch_versions.view()
    // println "HMF Data: ${hmf_data}"
    // println "ISOFOX stage enabled: ${run_config.stages.isofox}"
    // println "ISOFOX functions: ${params.isofox_functions}"
    // println "ISOFOX read length: ${isofox_read_length}"
    // println "ISOFOX counts: ${isofox_counts}"
    // println "ISOFOX GC ratios: ${isofox_gc_ratios}"

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