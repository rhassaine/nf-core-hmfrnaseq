# nf-core/hmfrnaseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [Alignment](#alignment) - STAR alignment and BAM processing
- [RSeQC](#rseqc) - BAM quality metrics and rRNA contamination check
- [Isofox](#isofox) - Transcript quantification and fusion detection
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Output directory structure

Results are organized by sample:

```
results/
├── <group_id>/
│   ├── alignments/rna/      # Aligned BAM files
│   ├── fastqc/              # FastQC reports
│   ├── isofox/              # Isofox output
│   ├── rseqc/               # RSeQC metrics
│   └── multiqc/             # Per-sample MultiQC report
├── multiqc/                 # Aggregated MultiQC report (all samples)
└── pipeline_info/           # Execution reports
```

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/alignments/rna/`
  - `*.markdup.bam`: Coordinate-sorted BAM file with duplicates marked.
  - `*.markdup.bam.bai`: BAM index file.
  - `star_logs/`: STAR alignment log files.

</details>

The pipeline aligns reads using [STAR](https://github.com/alexdobin/STAR) and processes the output through SAMtools sorting, Sambamba merging (for multi-lane samples), and GATK MarkDuplicates for duplicate marking.

### RSeQC

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/rseqc/bamstat/`
  - `*.bam_stat.txt`: BAM statistics including mapping rates and proper pairs.

- `<group_id>/rseqc/readduplication/`
  - `*.pos.DupRate.xls`: Position-based read duplication rate.
  - `*.seq.DupRate.xls`: Sequence-based read duplication rate.
  - `*.DupRate_plot.pdf`: Duplication rate plot.

- `<group_id>/rseqc/splitbam/`
  - `*.stat.txt`: rRNA contamination statistics showing reads overlapping rRNA regions.

</details>

[RSeQC](http://rseqc.sourceforge.net/) provides RNA-seq quality metrics. The splitbam module checks for rRNA contamination by counting reads overlapping rRNA regions defined in a BED file. Samples exceeding rRNA thresholds are excluded from Isofox analysis.

### Isofox

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/isofox/`
  - `*.isf.summary.csv`: Summary statistics including QC status and fragment counts.
  - `*.isf.gene_data.csv`: Gene-level expression data.
  - `*.isf.transcript_data.csv`: Transcript-level expression data.
  - `*.isf.fusions.csv`: Detected gene fusions.
  - `*.isf.alt_splice_junc.csv`: Alternative splice junctions.
  - `*.isf.retained_intron.csv`: Retained introns.

</details>

[Isofox](https://github.com/hartwigmedical/hmftools/tree/master/isofox) is a dedicated tool for RNA transcript analysis, part of the WiGiTS suite developed by the Hartwig Medical Foundation. It performs transcript quantification, alternative splicing detection, and fusion calling. Only samples passing rRNA QC are processed by Isofox.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `<group_id>/multiqc/`
  - `<group_id>_multiqc_report.html`: Per-sample report including FastQC, RSeQC, and GATK MarkDuplicates metrics.
  - `<group_id>_multiqc_report_data/`: Parsed statistics for the sample.

- `multiqc/`
  - `Aggregated_multiqc_report.html`: Aggregated report across all samples (RSeQC and MarkDuplicates metrics only, one row per sample).
  - `Aggregated_multiqc_report_data/`: Parsed statistics for all samples.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates HTML reports summarising QC metrics. The pipeline generates two types of reports:

1. **Per-sample reports** - Include all QC metrics (FastQC, RSeQC, MarkDuplicates) for each sample
2. **Aggregated report** - Combines RSeQC and MarkDuplicates metrics across all samples for easy comparison

For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
