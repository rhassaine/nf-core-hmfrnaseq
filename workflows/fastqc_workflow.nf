#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Workflow: FastQC + MultiQC
 * Location: workflows/fastqc_multiqc.nf
 * Description: Runs FastQC on all FASTQ files in the samplesheet, then generates a MultiQC report.
 */

// Import pipeline utilities and modules
import Utils
import Constants

include { FASTQC }                        from '../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_AGGREGATED } from '../modules/local/multiqc/main'

// Main workflow block
workflow FASTQC_WORKFLOW {

    /*
     * Input: params.input (CSV samplesheet path)
     * Output: FastQC HTML/ZIP, MultiQC report, versions
     */

    main:
        // Parse the samplesheet using pipeline logic
        inputs = Utils.parseInput(params.input, workflow.stubRun, log)
        ch_inputs = channel.fromList(inputs)

        // Extract FASTQ file tuples for FastQC
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

        // Run FastQC on all FASTQ pairs
        FASTQC(ch_fastq_for_qc)

        // Collect FastQC outputs for MultiQC
        // MULTIQC_AGGREGATED is the local per-sample-aware module (modules/local/multiqc)
        // aliased the same way rna_workflow.nf aliases its own single combined report —
        // this publishes flat to "${params.outdir}/multiqc" (conf/modules.config's
        // withName: 'MULTIQC_AGGREGATED' rule), not nested under a per-sample folder.
        // The fixed meta only feeds the report's filename prefix (FastQC_multiqc_report.*),
        // it no longer drives the publish path. Using the stock modules/nf-core/multiqc
        // module (no meta at all) here instead used to leave `meta` unbound in the
        // shared plain-MULTIQC publishDir directive, publishing into a literal "[:]" folder.
        ch_multiqc_files = channel.empty()
            .mix(FASTQC.out.zip.map { meta, file -> file })
            .collect()
            .map { files -> [[id: 'FastQC', key: 'FastQC'], files] }

        // Run MultiQC on FastQC outputs
        MULTIQC_AGGREGATED(
            ch_multiqc_files,
            [], // multiqc_config
            [], // extra_multiqc_config
            [], // multiqc_logo
            [], // replace_names
            []  // sample_names
        )

        // Collect version info
        // FASTQC now publishes via the `versions` topic channel (see
        // modules/nf-core/fastqc/main.nf), not a regular `versions` emit.
        ch_versions = channel.empty()
            .mix(MULTIQC_AGGREGATED.out.versions)

    emit:
        fastqc_html     = FASTQC.out.html
        fastqc_zip      = FASTQC.out.zip
        multiqc_report  = MULTIQC_AGGREGATED.out.report
        multiqc_data    = MULTIQC_AGGREGATED.out.data
        multiqc_plots   = MULTIQC_AGGREGATED.out.plots
        versions        = ch_versions
}
