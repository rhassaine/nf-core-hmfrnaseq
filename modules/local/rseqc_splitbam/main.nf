process RSEQC_SPLITBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_1' :
        'biocontainers/rseqc:5.0.4--pyhdfd78af_1' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.in.bam")  , emit: in_bam
    tuple val(meta), path("*.ex.bam")  , emit: ex_bam
    tuple val(meta), path("*.junk.bam"), emit: junk_bam
    tuple val(meta), path("*.stats.txt"), emit: stats
    tuple val(meta), path("*_mqc.yaml"), emit: multiqc
    tuple val("${task.process}"), val('rseqc'), eval('split_bam.py --version | sed "s/split_bam.py //"'), emit: versions_rseqc, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_unzip = bed.toString().endsWith('.gz') ? "gunzip -c ${bed} > bed_file.bed" : ""
    def bed_file = bed.toString().endsWith('.gz') ? "bed_file.bed" : bed
    """
    ${bed_unzip}
    split_bam.py \\
        -i $bam \\
        -r ${bed_file} \\
        -o $prefix \\
        $args > ${prefix}.stats.txt

    # Parse stats and create MultiQC-compatible YAML
    # Output format is "%-55s%d" - number follows immediately after text with no space
    # Use grep -oE to extract just the trailing digits
    total_reads=\$(grep "Total records:" ${prefix}.stats.txt | grep -oE '[0-9]+\$')
    rrna_reads=\$(grep "\\.in\\.bam" ${prefix}.stats.txt | grep -oE '[0-9]+\$')
    non_rrna_reads=\$(grep "\\.ex\\.bam" ${prefix}.stats.txt | grep -oE '[0-9]+\$')
    junk_reads=\$(grep "\\.junk\\.bam" ${prefix}.stats.txt | grep -oE '[0-9]+\$')

    # Default to 0 if empty
    total_reads=\${total_reads:-0}
    rrna_reads=\${rrna_reads:-0}
    non_rrna_reads=\${non_rrna_reads:-0}
    junk_reads=\${junk_reads:-0}

    # Calculate percentage (handle division by zero)
    if [ "\$total_reads" -gt 0 ]; then
        rrna_pct=\$(awk "BEGIN {printf \\"%.2f\\", (\$rrna_reads / \$total_reads) * 100}")
    else
        rrna_pct="0.00"
    fi

    cat > ${prefix}.splitbam_mqc.yaml <<EOF
id: 'rseqc_splitbam'
section_name: 'RSeQC: Split BAM (rRNA)'
description: 'rRNA contamination metrics from RSeQC split_bam.py using the provided BED file regions.'
plot_type: 'generalstats'
pconfig:
  - rrna_reads:
      title: 'rRNA Reads'
      description: 'Number of reads overlapping BED regions (rRNA)'
      min: 0
      format: '{:,.0f}'
      scale: 'RdYlGn-rev'
  - rrna_pct:
      title: 'rRNA %'
      description: 'Percentage of reads overlapping BED regions (rRNA)'
      min: 0
      max: 100
      suffix: '%'
      format: '{:,.2f}'
      scale: 'RdYlGn-rev'
  - non_rrna_reads:
      title: 'Non-rRNA Reads'
      description: 'Number of reads not overlapping BED regions'
      min: 0
      format: '{:,.0f}'
      hidden: true
  - junk_reads:
      title: 'Unmapped/QC-failed'
      description: 'Number of unmapped or QC-failed reads'
      min: 0
      format: '{:,.0f}'
      hidden: true
data:
  ${prefix}:
    rrna_reads: \$rrna_reads
    rrna_pct: \$rrna_pct
    non_rrna_reads: \$non_rrna_reads
    junk_reads: \$junk_reads
EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.in.bam
    touch ${prefix}.ex.bam
    touch ${prefix}.junk.bam
    touch ${prefix}.stats.txt

    cat > ${prefix}.splitbam_mqc.yaml <<EOF
id: 'rseqc_splitbam'
section_name: 'RSeQC: Split BAM (rRNA)'
description: 'rRNA contamination metrics from RSeQC split_bam.py using the provided BED file regions.'
plot_type: 'generalstats'
pconfig:
  - rrna_reads:
      title: 'rRNA Reads'
      description: 'Number of reads overlapping BED regions (rRNA)'
      min: 0
      format: '{:,.0f}'
      scale: 'RdYlGn-rev'
  - rrna_pct:
      title: 'rRNA %'
      description: 'Percentage of reads overlapping BED regions (rRNA)'
      min: 0
      max: 100
      suffix: '%'
      format: '{:,.2f}'
      scale: 'RdYlGn-rev'
data:
  ${prefix}:
    rrna_reads: 0
    rrna_pct: 0.0
EOF
    """
}
