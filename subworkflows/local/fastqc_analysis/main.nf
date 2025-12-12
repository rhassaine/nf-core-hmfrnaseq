#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import your existing Utils class for samplesheet parsing
import Utils

include { FASTQC }                from '../../../modules/nf-core/fastqc/main'

workflow FASTQC_WORKFLOW {
    /*
     * Input: params.input (CSV samplesheet path)
     * Output: FastQC HTML, ZIP, and versions channels
     */

    main:
        // Parse the samplesheet using your pipeline's logic
        inputs = Utils.parseInput(params.input, workflow.stubRun, log)

        // Create a channel from parsed input
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

        // Collect version info
        ch_versions = FASTQC.out.versions

    emit:
        html     = FASTQC.out.html
        zip      = FASTQC.out.zip
        versions = ch_versions
}

workflow {
    main:
        FASTQC_WORKFLOW()
}
