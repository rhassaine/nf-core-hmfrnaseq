process BAMCHECKER {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/rhassaine/bam-tools-non-arm:1.6' :
        'ghcr.io/rhassaine/bam-tools-non-arm:1.6' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path genome_fasta
    path genome_fai

    output:
    tuple val(meta), path("${meta.sample_id}.checked.bam"), emit: bam
    tuple val(meta), path("${meta.sample_id}.checked.bam.bai"), emit: bai
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def xmx_mod = task.ext.xmx_mod ?: 0.95

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        -cp /usr/share/java/bam-tools.jar \\
        com.hartwig.hmftools.bamtools.checker.BamChecker \\
        ${args} \\
        -bam_file ${bam} \\
        -ref_genome ${genome_fasta} \\
        -output_bam ${meta.sample_id}.checked.bam \\
        -bamtool \$(which samtools) \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamchecker: \$(java -cp /usr/share/java/bam-tools.jar com.hartwig.hmftools.bamtools.checker.BamChecker -version 2>&1 | grep -oP 'version \\K[0-9.]+' || echo 'unknown')
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.checked.bam
    touch ${meta.sample_id}.checked.bam.bai

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
