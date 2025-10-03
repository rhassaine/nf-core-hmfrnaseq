#!/usr/bin/env nextflow
import Constants
import Utils

nextflow.enable.dsl=2

// log.info "Constants class: ${Constants.class}"
// log.info "Utils class: ${Utils.class}"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/hmfrnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/hmfrnaseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { HMFRNASEQ  } from './workflows/hmfrnaseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hmfrnaseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hmfrnaseq_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_hmfrnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`

params.ref_data_genome_fasta         = getGenomeAttribute('fasta')
params.ref_data_genome_fai           = getGenomeAttribute('fai')
params.ref_data_genome_dict          = getGenomeAttribute('dict')
params.ref_data_genome_img           = getGenomeAttribute('img')
params.ref_data_genome_bwamem2_index = getGenomeAttribute('bwamem2_index')
params.ref_data_genome_gridss_index  = getGenomeAttribute('gridss_index')
params.ref_data_genome_star_index    = getGenomeAttribute('star_index')

WorkflowMain.setParamsDefaults(params, log)
WorkflowMain.validateParams(params, log)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE PLACEHOLDER FILES FOR STUB RUNS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// NOTE(SW): required prior to workflow import

if (workflow.stubRun && params.create_stub_placeholders) {
    Utils.createStubPlaceholders(params)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNA_WORKFLOW }    from './workflows/rna_workflow'
include { FASTQC_WORKFLOW } from './workflows/fastqc_workflow'

run_mode = Utils.getRunMode(params.mode, log)

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_HMFRNASEQ {

    if (run_mode === Constants.RunMode.RNA_WORKFLOW){
        RNA_WORKFLOW()
    } else if (run_mode === Constants.RunMode.FASTQC_WORKFLOW) {
        FASTQC_WORKFLOW()
    } else {
        log.error("received bad run mode: ${run_mode}")
        Nextflow.exit(1)
    }
}

//     take:
//     samplesheet // channel: samplesheet read in from --input

//     main:

//     //
//     // WORKFLOW: Run pipeline
//     //
//     HMFRNASEQ (
//         samplesheet
//     )
//     emit:
//     multiqc_report = HMFRNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
// }
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        // params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_HMFRNASEQ (
        // PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        // NFCORE_HMFRNASEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
