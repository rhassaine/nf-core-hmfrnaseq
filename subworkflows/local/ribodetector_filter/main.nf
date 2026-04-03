//
// Pre-alignment rRNA filtering using Ribodetector
// ML-based rRNA detection — no database required, just a read length parameter.
//

import Constants
import Utils

include { RIBODETECTOR } from '../../../modules/nf-core/ribodetector/main'

workflow RIBODETECTOR_FILTER {
    take:
    ch_inputs       // channel: [mandatory] [ meta ]
    read_length     // value:   [mandatory] read length for Ribodetector (e.g. 151)

    main:
    // Extract FASTQs from runnable samples
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = ch_inputs
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

    // nf-core ribodetector module expects [meta, [R1, R2]]
    ch_ribodetector_in = ch_fastq_inputs.map { meta, fwd, rev -> [meta, [fwd, rev]] }

    RIBODETECTOR(ch_ribodetector_in, read_length)

    // Split output back to [meta, R1, R2]
    ch_reads_filtered = RIBODETECTOR.out.fastq
        .map { meta, reads ->
            def (fwd, rev) = reads.sort()
            [meta, fwd, rev]
        }

    emit:
    reads    = ch_reads_filtered   // channel: [ meta_fastq, fwd, rev ]
    ribo_log = RIBODETECTOR.out.log // channel: [ meta, log_file ]
}
