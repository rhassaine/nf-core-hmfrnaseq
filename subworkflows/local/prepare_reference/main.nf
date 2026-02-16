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
    ch_genome_fasta = channel.fromPath(params.ref_data_genome_fasta)
    ch_genome_version = channel.value(params.genome_version)

    //
    // Set .fai and .dict indexes, create if required
    //
    ch_genome_fai = getRefFilechannel('ref_data_genome_fai')
    if (!params.ref_data_genome_fai) {
        SAMTOOLS_FAIDX(ch_genome_fasta)
        ch_genome_fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    ch_genome_dict = getRefFilechannel('ref_data_genome_dict')
    if (!params.ref_data_genome_dict) {
        SAMTOOLS_DICT(ch_genome_fasta)
        ch_genome_dict = SAMTOOLS_DICT.out.dict
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
            ch_genome_star_index = STAR_GENOMEGENERATE.out.index
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

        } else if (params.ref_data_genome_star_index.endsWith('.tar.gz')) {

            ch_genome_star_index_inputs = channel.fromPath(params.ref_data_genome_star_index)
                .map { [[id: "star_index_${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_STAR_INDEX(ch_genome_star_index_inputs)
            ch_genome_star_index = DECOMP_STAR_INDEX.out.extracted_dir

        } else {

            ch_genome_star_index = getRefFilechannel('ref_data_genome_star_index')

        }
    }

    //
    // Set SortMeRNA rRNA database, unpack if required
    //
    ch_sortmerna_db = channel.empty()
    if (params.sortmerna_fastas && run_config.has_rna_fastq && run_config.stages.alignment) {
        if (params.sortmerna_fastas.endsWith('.tar.gz')) {

            ch_sortmerna_db_inputs = channel.fromPath(params.sortmerna_fastas)
                .map { [[id: 'sortmerna_db'], it] }

            DECOMP_SORTMERNA_DB(ch_sortmerna_db_inputs)

            ch_sortmerna_db = DECOMP_SORTMERNA_DB.out.extracted_dir
                .map { dir -> file("${dir}/smr_v4.3_default_db.fasta") }

        } else {

            ch_sortmerna_db = channel.fromPath(params.sortmerna_fastas)

        }
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
    genome_fasta         = ch_genome_fasta.first()         // path: genome_fasta
    genome_fai           = ch_genome_fai.first()           // path: genome_fai
    genome_dict          = ch_genome_dict.first()          // path: genome_dict
    genome_star_index    = ch_genome_star_index.first()    // path: genome_star_index
    genome_version       = ch_genome_version               // val:  genome_version

    hmf_data             = ch_hmf_data                     // map:  HMF data paths
    sortmerna_db         = ch_sortmerna_db                 // path: sortmerna rRNA database fasta

    versions             = ch_versions                     // channel: [ versions.yml ]
}

def getRefFilechannel(key) {
    def fp = params.get(key) ? file(params.getAt(key)) : []
    return channel.of(fp)
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
