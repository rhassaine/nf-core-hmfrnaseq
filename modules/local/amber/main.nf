process AMBER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-amber:4.2--hdfd78af_0' :
        'biocontainers/hmftools-amber:4.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val genome_ver
    path heterozygous_sites

    output:
    tuple val(meta), path('amber/'), emit: amber_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def xmx_mod = task.ext.xmx_mod ?: 0.75

    // Germline-only mode: uses -reference/-reference_bam instead of -tumor/-tumor_bam
    // This skips contamination checks and PCF fitting, and names output with ref sample
    // See: https://github.com/hartwigmedical/hmftools/tree/master/amber#germline-only-mode
    """
    amber \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -reference ${meta.sample_id} \\
        -reference_bam ${bam} \\
        -ref_genome_version ${genome_ver} \\
        -loci ${heterozygous_sites} \\
        -threads ${task.cpus} \\
        -output_dir amber/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: \$(amber -version | sed -n '/^Amber version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p amber/
    touch amber/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
