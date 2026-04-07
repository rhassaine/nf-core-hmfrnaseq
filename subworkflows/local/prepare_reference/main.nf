//
// Prepare reference data as required
//

import Constants

include { SAMTOOLS_DICT         } from '../../../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE   } from '../../../modules/nf-core/star/genomegenerate/main'

include { CUSTOM_EXTRACTTARBALL as DECOMP_HMF_DATA         } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_SORTMERNA_DB     } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_STAR_INDEX       } from '../../../modules/local/custom/extract_tarball/main'

workflow PREPARE_REFERENCE {
    take:
    run_config // channel: [mandatory] run configuration

    main:
    // channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = channel.empty()

    //
    // Set genome reference channels
    //
    ch_genome_fasta = channel.fromPath(params.ref_data_genome_fasta).first()
    ch_genome_version = channel.value(params.genome_version)

    //
    // Set .fai and .dict indexes, create if required
    //
    ch_genome_fai = getRefFilechannel('ref_data_genome_fai')
    if (!params.ref_data_genome_fai) {
        SAMTOOLS_FAIDX(ch_genome_fasta)
        ch_genome_fai = SAMTOOLS_FAIDX.out.fai.first()
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    ch_genome_dict = getRefFilechannel('ref_data_genome_dict')
    if (!params.ref_data_genome_dict) {
        SAMTOOLS_DICT(ch_genome_fasta)
        ch_genome_dict = SAMTOOLS_DICT.out.dict.first()
        ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)
    }

    //
    // Set STAR index, unpack or create if required
    //
    ch_genome_star_index = channel.empty()
    if (run_config.has_rna_fastq && run_config.stages.alignment) {
        if (!params.ref_data_genome_star_index) {

            STAR_GENOMEGENERATE(
                ch_genome_fasta,
                file(params.ref_data_genome_gtf),
            )
            ch_genome_star_index = STAR_GENOMEGENERATE.out.index.first()
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

        } else if (params.ref_data_genome_star_index.endsWith('.tar.gz')) {

            ch_genome_star_index_inputs = channel.fromPath(params.ref_data_genome_star_index)
                .map { [[id: "star_index_${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_STAR_INDEX(ch_genome_star_index_inputs)
            ch_genome_star_index = DECOMP_STAR_INDEX.out.extracted_dir.first()

        } else {

            ch_genome_star_index = getRefFilechannel('ref_data_genome_star_index')

        }
    }

    //
    // Set SortMeRNA rRNA database, unpack if required
    //
    ch_sortmerna_db = channel.empty()
    if (params.sortmerna_fastas && params.mode != 'rna_redux_workflow' && run_config.has_rna_fastq && run_config.stages.alignment) {
        if (params.sortmerna_fastas.endsWith('.tar.gz')) {

            ch_sortmerna_db_inputs = channel.fromPath(params.sortmerna_fastas)
                .map { [[id: 'sortmerna_db'], it] }

            DECOMP_SORTMERNA_DB(ch_sortmerna_db_inputs)

            ch_sortmerna_db = DECOMP_SORTMERNA_DB.out.extracted_dir
                .map { dir -> file("${dir}/smr_v4.3_default_db.fasta") }
                .first()

        } else {

            ch_sortmerna_db = channel.fromPath(params.sortmerna_fastas).first()

        }
    }

    //
    // Set genome GTF for RustQC biotype counting
    //
    ch_genome_gtf = channel.empty()
    if (params.ref_data_genome_gtf && run_config.stages.rseqc) {
        ch_genome_gtf = channel.fromPath(params.ref_data_genome_gtf).first()
    }

    //
    // Set HMF reference data, unpack if required
    //
    ch_hmf_data = channel.empty()
    hmf_data_paths = params.hmf_data_paths[params.genome_version.toString()]
    if (params.ref_data_hmf_data_path.endsWith('tar.gz')) {

        ch_hmf_data_inputs = channel.fromPath(params.ref_data_hmf_data_path)
            .map { [[id: "hmf_data_${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

        DECOMP_HMF_DATA(ch_hmf_data_inputs)

        ch_hmf_data = DECOMP_HMF_DATA.out.extracted_dir
            .collect()
            .map { dir_list ->
                assert dir_list.size() == 1
                def dirpath = dir_list[0].toUriString()
                return createDataMap(hmf_data_paths, dirpath)
            }

    } else {

        ch_hmf_data = channel.value(createDataMap(hmf_data_paths, params.ref_data_hmf_data_path))

    }

    emit:
    genome_fasta         = ch_genome_fasta                  // path: genome_fasta
    genome_fai           = ch_genome_fai                   // path: genome_fai
    genome_dict          = ch_genome_dict                  // path: genome_dict
    genome_star_index    = ch_genome_star_index            // path: genome_star_index
    genome_version       = ch_genome_version               // val:  genome_version

    hmf_data             = ch_hmf_data                     // map:  HMF data paths
    sortmerna_db         = ch_sortmerna_db                 // path: sortmerna rRNA database fasta
    genome_gtf           = ch_genome_gtf                  // path: genome GTF for RustQC biotype counting

    versions             = ch_versions                     // channel: [ versions.yml ]
}

// IMPORTANT: must use channel.value() not channel.of() here.
//
// Nextflow has two channel types (https://www.nextflow.io/docs/latest/channel.html):
//
//   Queue channel (channel.of, channel.fromPath, channel.fromList, process outputs):
//     - Each item can only be consumed ONCE across all readers
//     - When a process reads an item, it's gone — other tasks/processes never see it
//     - Designed for per-sample data that flows through the pipeline
//
//   Value channel (channel.value, channel.of with single item in some contexts):
//     - Can be read unlimited times by any number of processes/tasks
//     - The same value is available to every task invocation
//     - Designed for shared reference data (genomes, indices, BED files)
//
// Using channel.of() here created a queue channel for reference files like
// genome_star_index. The first STAR task consumed it, and all subsequent
// samples silently received nothing — causing only 1/N samples to align.
//
// See: https://www.nextflow.io/docs/latest/channel.html#queue-channel
//      https://www.nextflow.io/docs/latest/channel.html#value-channel
def getRefFilechannel(key) {
    def fp = params.get(key) ? file(params.getAt(key)) : []
    return channel.value(fp) // was channel.of(fp) — the bug that caused only 1/N samples to align when using a STAR index tarball
}

def createDataMap(entries, ref_data_path) {
    return entries
        .collectEntries { name, path ->
            def ref_data_file = path == [] ? [] : getRefdataFile(path, ref_data_path)
            return [name, ref_data_file]
        }
}

def getRefdataFile(filepath, ref_data_path) {
    // If filepath is already absolute (starts with /), use it as-is
    // Otherwise, concatenate with base path
    def full_path = filepath.startsWith('/') ? filepath : "${ref_data_path.toString().replaceAll('/$', '')}/${filepath}"
    return file(full_path, checkIfExists: true)
}
