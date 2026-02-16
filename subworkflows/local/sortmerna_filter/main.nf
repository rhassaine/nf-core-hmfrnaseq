//
// Pre-alignment rRNA filtering using SortMeRNA
//

import Constants
import Utils

include { SORTMERNA } from '../../../modules/nf-core/sortmerna/main'

workflow SORTMERNA_FILTER {
    take:
    ch_inputs       // channel: [mandatory] [ meta ]
    sortmerna_db    // channel: [mandatory] /path/to/sortmerna_rRNA_db.fasta

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

    // nf-core sortmerna module expects [meta, [R1, R2]]
    ch_sortmerna_in = ch_fastq_inputs.map { meta, fwd, rev -> [meta, [fwd, rev]] }
    ch_sortmerna_fastas = sortmerna_db.map { db -> [[:], db] }

    SORTMERNA(ch_sortmerna_in, ch_sortmerna_fastas, [[:], []])

    // Note: SORTMERNA versions are collected via topics

    // Split output back to [meta, R1, R2]
    ch_reads_filtered = SORTMERNA.out.reads
        .map { meta, reads ->
            def (fwd, rev) = reads.sort()
            [meta, fwd, rev]
        }

    emit:
    reads    = ch_reads_filtered  // channel: [ meta_fastq, fwd, rev ]
    sort_log = SORTMERNA.out.log  // channel: [ meta, log_file ]
}
