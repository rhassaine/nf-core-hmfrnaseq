process WRITE_REFERENCE_DATA {
    tag "${fp.name}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
    path fp
    val workflow_version

    output:
    path fp, includeInputs: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """
}