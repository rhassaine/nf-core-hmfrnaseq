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

include { FASTQC }   from '../modules/nf-core/fastqc/main'
include { MULTIQC }  from '../modules/nf-core/multiqc/main'

// Main workflow block
workflow FASTQC_WORKFLOW {

    /*
     * Input: params.input (CSV samplesheet path)
     * Output: FastQC HTML/ZIP, MultiQC report, versions
     */

    main:
        // Parse the samplesheet using pipeline logic
        inputs = Utils.parseInput(params.input, workflow.stubRun, log)
        ch_inputs = Channel.fromList(inputs)

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
        ch_multiqc_files = Channel.empty()
            .mix(FASTQC.out.html.map { meta, file -> file })
            .mix(FASTQC.out.zip.map { meta, file -> file })
            .collect()

        // Run MultiQC on FastQC outputs
        MULTIQC(
            ch_multiqc_files,
            [], // multiqc_config
            [], // extra_multiqc_config
            [], // multiqc_logo
            [], // replace_names
            []  // sample_names
        )

        // Collect version info
        ch_versions = Channel.empty()
            .mix(FASTQC.out.versions)
            .mix(MULTIQC.out.versions)

    emit:
        fastqc_html     = FASTQC.out.html
        fastqc_zip      = FASTQC.out.zip
        multiqc_report  = MULTIQC.out.report
        multiqc_data    = MULTIQC.out.data
        multiqc_plots   = MULTIQC.out.plots
        versions        = ch_versions
}