import Constants
import Utils

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REDUX TEST WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Minimal test: run REDUX on a BAM file
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { REDUX                  } from '../modules/local/redux/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REDUX_TEST {
    // Create channel for versions
    ch_versions = channel.empty()

    //
    // Reference data
    //
    genome_fasta = file(params.ref_data_genome_fasta)
    genome_fai   = file(params.ref_data_genome_fai)
    genome_dict  = file(params.ref_data_genome_dict)
    genome_ver   = params.genome_version

    // HMF data paths
    hmf_data_path = params.ref_data_hmf_data_path
    hmf_paths = params.hmf_data_paths[genome_ver.toString()]
    unmap_regions    = file("${hmf_data_path}/${hmf_paths.unmap_regions}")
    msi_jitter_sites = file("${hmf_data_path}/${hmf_paths.msi_jitter_sites}")

    //
    // Input BAM
    //
    bam_file = file(params.redux_test_bam)
    bai_file = file("${params.redux_test_bam}.bai")

    ch_redux_input = channel.of([
        [key: 'test', id: 'redux_test', sample_id: 'test_sample'],
        [bam_file],
        [bai_file]
    ])

    //
    // Run REDUX
    //
    REDUX(
        ch_redux_input,
        genome_fasta,
        genome_ver,
        genome_fai,
        genome_dict,
        unmap_regions,
        msi_jitter_sites,
        false,  // umi_enable
        '',     // umi_duplex_delim
    )

    ch_versions = ch_versions.mix(REDUX.out.versions)

    //
    // Aggregate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true,
        )
}
