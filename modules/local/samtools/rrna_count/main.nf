process SAMTOOLS_RRNA_COUNT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.rrna.stats.txt"), emit: stats
    tuple val(meta), path("*_mqc.yaml")      , emit: multiqc
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | sed -n "1s/samtools //p"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_unzip = bed.toString().endsWith('.gz') ? "gunzip -c ${bed} > rrna_regions.bed" : ""
    def bed_file  = bed.toString().endsWith('.gz') ? "rrna_regions.bed" : bed
    """
    ${bed_unzip}

    # Count primary, mapped, non-supplementary reads (-F 0x904 = exclude unmapped 0x4 +
    # secondary 0x100 + supplementary 0x800) overlapping the rRNA regions vs the total.
    # Same interval-overlap method as RSeQC split_bam, but a single samtools pass and no
    # in/ex/junk BAMs written -> much faster.
    rrna_reads=\$(samtools view -c -F 0x904 -L ${bed_file} --threads ${task.cpus} ${bam})
    total_reads=\$(samtools view -c -F 0x904 --threads ${task.cpus} ${bam})
    non_rrna_reads=\$(( total_reads - rrna_reads ))

    if [ "\$total_reads" -gt 0 ]; then
        rrna_pct=\$(awk "BEGIN {printf \\"%.2f\\", (\$rrna_reads / \$total_reads) * 100}")
    else
        rrna_pct="0.00"
    fi

    printf 'Total records: %s\\nrRNA reads: %s\\nnon-rRNA reads: %s\\nrRNA pct: %s\\n' \\
        "\$total_reads" "\$rrna_reads" "\$non_rrna_reads" "\$rrna_pct" > ${prefix}.rrna.stats.txt

# MultiQC custom-content YAML (column-0 so the emitted YAML is not indented).
cat > ${prefix}.rrna_mqc.yaml <<EOF
id: 'rseqc_splitbam'
section_name: 'rRNA content (samtools)'
description: 'rRNA fraction from primary mapped reads overlapping the rRNA BED regions (samtools view -L); same interval method as RSeQC split_bam.'
plot_type: 'generalstats'
pconfig:
  - rrna_reads:
      title: 'rRNA Reads'
      description: 'Primary mapped reads overlapping rRNA regions'
      min: 0
      format: '{:,.0f}'
      scale: 'RdYlGn-rev'
  - rrna_pct:
      title: 'rRNA %'
      description: 'Percentage of primary mapped reads overlapping rRNA regions'
      min: 0
      max: 100
      suffix: '%'
      format: '{:,.2f}'
      scale: 'RdYlGn-rev'
  - non_rrna_reads:
      title: 'Non-rRNA Reads'
      description: 'Primary mapped reads not overlapping rRNA regions'
      min: 0
      format: '{:,.0f}'
      hidden: true
data:
  ${prefix}:
    rrna_reads: \$rrna_reads
    rrna_pct: \$rrna_pct
    non_rrna_reads: \$non_rrna_reads
EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rrna.stats.txt
cat > ${prefix}.rrna_mqc.yaml <<EOF
id: 'rseqc_splitbam'
section_name: 'rRNA content (samtools)'
plot_type: 'generalstats'
data:
  ${prefix}:
    rrna_reads: 0
    rrna_pct: 0.0
    non_rrna_reads: 0
EOF
    """
}
