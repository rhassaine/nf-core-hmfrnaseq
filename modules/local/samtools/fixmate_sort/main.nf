process SAMTOOLS_FIXMATE_SORT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::samtools=1.22.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"),     emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai
    path 'versions.yml',                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools sort \\
        -n \\
        -@ ${task.cpus} \\
        ${bam} \\
    | samtools fixmate \\
        -m \\
        -u \\
        - - \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.bam##idx##${prefix}.bam.bai \\
        --write-index \\
        -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | sed -n '1s/^.*samtools //p')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
