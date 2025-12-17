//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import nextflow.Nextflow

class Utils {

    public static parseInput(input_fp_str, stub_run, log) {

        // NOTE(SW): using NF .splitCsv channel operator, hence should be easily interchangable with NF syntax

        def input_fp = Utils.getFileObject(input_fp_str)
        def inputs = nextflow.splitter.SplitterEx.splitCsv(input_fp, [header: true])
            .groupBy { it['group_id'] }
            .collect { group_id, entries ->

                def meta = [group_id: group_id]
                def sample_keys = [] as Set

                // Process each entry
                entries.each {
                    // Add subject id if absent or check if current matches existing
                    if (meta.containsKey('subject_id') && meta.subject_id != it.subject_id) {
                        log.error "got unexpected subject name for ${group_id} ${meta.subject_id}: ${it.subject_id}"
                        Nextflow.exit(1)
                    } else {
                        meta.subject_id = it.subject_id
                    }

                    // Sample type
                    def sample_type_enum = Utils.getEnumFromString(it.sample_type, Constants.SampleType)
                    if (!sample_type_enum) {
                        def sample_type_str = Utils.getEnumNames(Constants.SampleType).join('\n  - ')
                        log.error "received invalid sample type: '${it.sample_type}'. Valid options are:\n  - ${sample_type_str}"
                        Nextflow.exit(1)
                    }

                    // Sequence type
                    def sequence_type_enum = Utils.getEnumFromString(it.sequence_type, Constants.SequenceType)
                    if (!sequence_type_enum) {
                        def sequence_type_str = Utils.getEnumNames(Constants.SequenceType).join('\n  - ')
                        log.error "received invalid sequence type: '${it.sequence_type}'. Valid options are:\n  - ${sequence_type_str}"
                        Nextflow.exit(1)
                    }

                    // Filetype
                    def filetype_enum = Utils.getEnumFromString(it.filetype, Constants.FileType)
                    if (!filetype_enum) {
                        def filetype_str = Utils.getEnumNames(Constants.FileType).join('\n  - ')
                        log.error "received invalid file type: '${it.filetype}'. Valid options are:\n  - ${filetype_str}"
                        Nextflow.exit(1)
                    }

                    def sample_key = [sample_type_enum, sequence_type_enum]
                    def meta_sample = meta.get(sample_key, [sample_id: it.sample_id])

                    if (meta_sample.sample_id != it.sample_id) {
                        log.error "got unexpected sample name for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${it.sample_id}"
                        Nextflow.exit(1)
                    }

                    if (meta_sample.containsKey(filetype_enum) & filetype_enum != Constants.FileType.FASTQ) {
                        log.error "got duplicate file for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${filetype_enum}"
                        Nextflow.exit(1)
                    }

                    // Info data
                    def info_data = [:]
                    if (it.containsKey('info')) {
                        // Parse
                        it.info
                            .tokenize(';')
                            .each { e ->
                                def (k, v) = e.tokenize(':')
                                def info_field_enum = Utils.getEnumFromString(k, Constants.InfoField)

                                if (!info_field_enum) {
                                    def info_field_str = Utils.getEnumNames(Constants.InfoField).join('\n  - ')
                                    log.error "received invalid info field: '${k}'. Valid options are:\n  - ${info_field_str}"
                                    Nextflow.exit(1)
                                }

                                if (info_data.containsKey(info_field_enum)) {
                                    log.error "got duplicate info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${info_field_enum}"
                                    Nextflow.exit(1)
                                }

                                info_data[info_field_enum] = v
                            }

                    }


                    // Handle inputs appropriately
                    if (filetype_enum === Constants.FileType.FASTQ) {

                        if (!info_data.containsKey(Constants.InfoField.LIBRARY_ID)) {
                            log.error "missing 'library_id' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }

                        if (!info_data.containsKey(Constants.InfoField.LANE)) {
                            log.error "missing 'lane' info field for ${group_id} ${sample_type_enum}/${sequence_type_enum}"
                            Nextflow.exit(1)
                        }

                        def (fwd, rev) = it.filepath.tokenize(';')
                        def fastq_key = [info_data[Constants.InfoField.LIBRARY_ID], info_data[Constants.InfoField.LANE]]

                        if (meta_sample.containsKey(fastq_key)) {
                            log.error "got duplicate lane + library_id data for ${group_id} ${sample_type_enum}/${sequence_type_enum}: ${fastq_key}"
                            Nextflow.exit(1)
                        }

                        if (!meta_sample.containsKey(filetype_enum)) {
                            meta_sample[filetype_enum] = [:]
                        }

                        meta_sample[filetype_enum][fastq_key] = ['fwd': fwd, 'rev': rev]

                    } else {

                        meta_sample[filetype_enum] = Utils.getFileObject(it.filepath)

                    }

                    // Record sample key to simplify iteration later on
                    sample_keys << sample_key
                }

                // Check that required indexes are provided or are accessible
                sample_keys.each { sample_key ->

                    meta[sample_key]*.key.each { key ->

                        def index_enum
                        def index_str

                        if (key === Constants.FileType.BAM) {
                            index_enum = Constants.FileType.BAI
                            index_str = 'bai'
                        } else {
                            return
                        }

                        if (meta[sample_key].containsKey(index_enum)) {
                            return
                        }

                        def fp = meta[sample_key][key].toUriString()
                        def index_fp = nextflow.Nextflow.file("${fp}.${index_str}")

                        if (!index_fp.exists() && !stub_run) {
                            def (sample_type, sequence_type) = sample_key
                            log.error "no index provided or found for ${meta.group_id} ${sample_type}/${sequence_type}: ${key}: ${fp}"
                            Nextflow.exit(1)
                        }

                        meta[sample_key][index_enum] = index_fp

                    }
                }

                return meta
            }

        return inputs
    }

    public static void createStubPlaceholders(params) {

        def fps = [
            params.ref_data_genome_dict,
            params.ref_data_genome_fai,
            params.ref_data_genome_fasta,
            params.ref_data_genome_gtf,
            params.ref_data_genome_star_index,
        ]

        params.hmf_data_paths[params.genome_version.toString()]
            .each { k, v ->
                fps << "${params.ref_data_hmf_data_path.replaceAll('/$', '')}/${v}"
            }

        fps.each { fp_str ->
            if (fp_str === null) return

            def fp = Utils.getFileObject(fp_str)

            if (!fp_str || fp.exists()) return

            if (fp_str.endsWith('/')) {
                fp.mkdirs()
            } else {
                fp.getParent().mkdirs()
                fp.toFile().createNewFile()
            }
        }

    }

    public static void validateInput(inputs, run_config, params, log) {

        def sample_keys = [
            [Constants.SampleType.TUMOR, Constants.SequenceType.RNA],
        ]

        inputs.each { meta ->

            // Require BAMs or FASTQs for each defined sample type
            sample_keys.each { key ->

                if (!meta.containsKey(key)) {
                    return
                }

                def (sample_type, sequence_type) = key

                if (!meta[key].containsKey(Constants.FileType.BAM) &&
                    !meta[key].containsKey(Constants.FileType.FASTQ)) {

                    log.error "no BAMs nor FASTQs provided for ${meta.group_id} ${sample_type}/${sequence_type}\n\n" +
                        "NB: BAMs or FASTQs are always required as they are the basis to determine input sample type."
                    Nextflow.exit(1)
                }

            }

            // Do not allow CRAM RNA input
            if (Utils.hasTumorRnaBam(meta) && Utils.getTumorRnaBam(meta).toString().endsWith('cram')) {
                log.error "found tumor RNA CRAM input for ${meta.group_id} but RNA CRAM input is not supported"
                Nextflow.exit(1)
            }

        }

        // Require that an input GTF file is provided when creating STAR index
        def run_star_index = run_config.stages.alignment && run_config.has_rna_fastq && !params.ref_data_genome_star_index

        if (run_star_index && !params.ref_data_genome_gtf) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Creating a STAR index requires the appropriate genome transcript annotations\n" +
                "  as a GTF file. Please contact us on Slack for further information."
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

    }

    static public getEnumFromString(s, e) {
        try {
            return e.valueOf(s.toUpperCase())
        } catch(java.lang.IllegalArgumentException err) {
            return null
        }
    }

    public static getEnumNames(e) {
        e
            .values()
            *.name()
            *.toLowerCase()
    }


    static public getFileObject(path) {
        return path ? nextflow.Nextflow.file(path) : []
    }

    static public getRunMode(run_mode, log) {
        def run_mode_enum = Utils.getEnumFromString(run_mode, Constants.RunMode)
        if (!run_mode_enum) {
            def run_modes_str = Utils.getEnumNames(Constants.RunMode).join('\n  - ')
            log.error "received an invalid run mode: '${run_mode}'. Valid options are:\n  - ${run_modes_str}"
            Nextflow.exit(1)
        }
        return run_mode_enum
    }


    // Sample records
    static public getTumorRnaSample(meta) {
        return meta.getOrDefault([Constants.SampleType.TUMOR, Constants.SequenceType.RNA], [:])
    }

    // Sample names
    static public getTumorRnaSampleName(meta) {
        return getTumorRnaSample(meta)['sample_id']
    }

    // Files - Tumor RNA
    static public getTumorRnaFastq(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.FASTQ, null)
    }

    static public getTumorRnaBam(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.BAM, null)
    }

    static public getTumorRnaBai(meta) {
        return getTumorRnaSample(meta).getOrDefault(Constants.FileType.BAI, null)
    }


    static public hasTumorRnaFastq(meta) {
        return getTumorRnaFastq(meta) !== null
    }

    static public hasTumorRnaBam(meta) {
        return getTumorRnaBam(meta) !== null
    }


    // Status
    static public hasTumorRna(meta) {
        return hasTumorRnaBam(meta) || hasTumorRnaFastq(meta)
    }


    // Misc
    public static getInput(meta, key) {

        def result = []
        def (key_filetype, key_filetypes, key_sequencetypes) = key

        for (key_sample in [key_filetypes, key_sequencetypes].combinations()) {
            if (meta.containsKey(key_sample) && meta[key_sample].containsKey(key_filetype)) {
                result = meta[key_sample].get(key_filetype)
                break
            }
        }
        return result
    }

    public static hasExistingInput(meta, key) {
        return getInput(meta, key) != []
    }

    public static selectCurrentOrExisting(val, meta, key) {
        if (hasExistingInput(meta, key)) {
            return getInput(meta, key)
        } else {
            return val
        }
    }

    // Parse splitbam stats file and return rRNA metrics
    public static parseRrnaStats(stats_file) {
        def total_reads = 0
        def rrna_reads = 0

        stats_file.eachLine { line ->
            if (line.contains('Total records:')) {
                total_reads = line.split(':')[-1].trim() as Long
            } else if (line.contains('consumed by input gene list')) {
                rrna_reads = line.split(':')[-1].trim() as Long
            }
        }

        def rrna_percent = total_reads > 0 ? (rrna_reads / total_reads) * 100.0 : 0.0

        return [
            total_reads: total_reads,
            rrna_reads: rrna_reads,
            rrna_percent: rrna_percent
        ]
    }

    // Check if sample passes rRNA QC thresholds
    public static checkRrnaQc(rrna_stats, threshold_count, threshold_percent) {
        def pass = true
        def fail_reasons = []

        // Check absolute count threshold (0 = disabled)
        if (threshold_count > 0 && rrna_stats.rrna_reads > threshold_count) {
            pass = false
            fail_reasons << "rRNA reads (${rrna_stats.rrna_reads}) exceeds threshold (${threshold_count})"
        }

        // Check percentage threshold (0 = disabled)
        if (threshold_percent > 0 && rrna_stats.rrna_percent > threshold_percent) {
            pass = false
            fail_reasons << "rRNA percentage (${String.format('%.2f', rrna_stats.rrna_percent)}%) exceeds threshold (${threshold_percent}%)"
        }

        return [
            pass: pass,
            fail_reason: fail_reasons.join('; ') ?: 'none'
        ]
    }

}
