//
// rRNA QC gate: parse RSeQC splitbam stats, log SortMeRNA stats,
// and gate BAM/input channels for downstream Isofox processing
//

import Utils

workflow RRNA_QC_GATE {
    take:
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_splitbam_stats      // channel: [mandatory] [ meta, stats_file ]
    ch_sortmerna_log       // channel: [optional]  [ meta, log_file ]
    ch_bam                 // channel: [mandatory] [ meta, bam, bai ] (BAMs to gate)
    rrna_threshold_count   //   value: [mandatory] max rRNA read count (0 = disabled)
    rrna_threshold_percent //   value: [mandatory] max rRNA percentage (0 = disabled)

    main:
    //
    // STEP 1: Parse splitbam stats and branch samples into pass/fail
    //
    ch_rrna_qc_result = ch_splitbam_stats
        .map { meta, stats_file ->
            def rrna_stats = Utils.parseRrnaStats(stats_file)
            def qc_result = Utils.checkRrnaQc(
                rrna_stats,
                rrna_threshold_count ?: 0,
                rrna_threshold_percent ?: 0
            )
            // Log the rRNA stats for visibility
            log.info "rRNA QC [${meta.group_id}]: ${rrna_stats.rrna_reads}/${rrna_stats.total_reads} reads (${String.format('%.2f', rrna_stats.rrna_percent)}%) - ${qc_result.pass ? 'PASS' : 'FAIL'}"
            if (!qc_result.pass) {
                log.warn "Sample ${meta.group_id} FAILED rRNA QC: ${qc_result.fail_reason}"
            }
            [meta, qc_result.pass, rrna_stats]
        }
        .branch { meta, pass, rrna_stats ->
            pass: pass
                return meta
            fail: true
                return meta
        }

    //
    // STEP 2: SortMeRNA informational logging (aggregate per sample, log only)
    //
    ch_sortmerna_qc_result = ch_sortmerna_log
        .map { meta, log_file ->
            def stats = Utils.parseSortmernaLog(log_file)
            [meta.key, stats.total_reads, stats.rrna_reads]
        }
        .groupTuple()
        .map { key, total_list, rrna_list ->
            def total = total_list.sum()
            def rrna = rrna_list.sum()
            def pct = total > 0 ? (rrna / (double)total) * 100.0 : 0.0
            [key, [total_reads: total, rrna_reads: rrna, rrna_percent: pct]]
        }

    ch_inputs
        .map { meta -> [meta.group_id, meta] }
        .join(ch_sortmerna_qc_result)
        .map { group_id, meta, rrna_stats ->
            log.info "SortMeRNA [${meta.group_id}]: ${rrna_stats.rrna_reads}/${rrna_stats.total_reads} rRNA reads removed (${String.format('%.2f', rrna_stats.rrna_percent)}%)"
        }

    //
    // STEP 3: Duplicate the pass channel for two consumers (fixes queue channel bug)
    //
    ch_rrna_qc_result.pass
        .multiMap { meta ->
            for_bam: meta
            for_inputs: meta
        }
        .set { ch_pass_split }

    //
    // STEP 4: Gate BAMs - only pass samples proceed to Isofox
    //
    // channel: [ meta, bam, bai ]
    ch_bam_pass = ch_bam
        .map { meta, bam, bai -> [meta.group_id, meta, bam, bai] }
        .join(ch_pass_split.for_bam.map { meta -> [meta.group_id, meta] }, by: 0)
        .map { group_id, meta_bam, bam, bai, meta_qc -> [meta_bam, bam, bai] }

    //
    // STEP 5: Gate inputs - only pass sample metas proceed to Isofox
    //
    // channel: [ meta ]
    ch_inputs_pass = ch_inputs
        .map { meta -> [meta.group_id, meta] }
        .join(ch_pass_split.for_inputs.map { meta -> [meta.group_id] }, by: 0)
        .map { group_id, meta -> meta }

    emit:
    bam_pass    = ch_bam_pass    // channel: [ meta, bam, bai ] - BAMs that passed rRNA QC
    inputs_pass = ch_inputs_pass // channel: [ meta ] - input metas that passed rRNA QC
}
