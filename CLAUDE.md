# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

nf-core/hmfrnaseq is a Nextflow DSL2 bioinformatics pipeline for RNA-seq analysis, built on the nf-core framework. It integrates Hartwig Medical Foundation (HMF) tools for RNA analysis including transcript quantification and fusion detection.

## Common Commands

### Running the Pipeline

```bash
# Run with test profile (docker or singularity)
nextflow run . -profile test,docker --outdir results
nextflow run . -profile test,singularity --outdir results

# Full run with custom input
nextflow run . -profile docker --input samplesheet.csv --outdir results --genome GRCh38_hmf

# Stub run (dry-run with placeholders, no actual processing)
nextflow run . -profile test_stub,docker --outdir results -stub

# Resume a failed/interrupted run
nextflow run . -profile test,docker --outdir results -resume
```

### Linting and Validation

```bash
# Run nf-core linting (requires nf-core tools)
nf-core lint

# Validate parameters
nextflow run . --help
nextflow run . --help_full  # Show all parameters including hidden
```

### Testing

The pipeline uses nf-test for module and subworkflow testing. Test files (`.nf.test`) are located in `tests/` subdirectories alongside modules.

```bash
# Run a specific module test (must cd to module directory)
cd modules/local/rseqc_splitbam && nf-test test tests/main.nf.test
cd modules/nf-core/fastqc && nf-test test tests/main.nf.test
```

## Architecture

### Run Modes

The pipeline supports multiple run modes controlled by `--mode` parameter (default: `rna_workflow`):
- `rna_workflow`: Full RNA analysis (alignment, QC, Isofox quantification)
- `fastqc_workflow`: FastQC only

### Entry Points

- `main.nf` - Pipeline entry point, dispatches to workflows based on mode
- `workflows/rna_workflow.nf` - Main RNA analysis workflow
- `workflows/fastqc_workflow.nf` - FastQC-only workflow

### Pipeline Stages (rna_workflow)

1. **PREPARE_REFERENCE** - Downloads/stages HMF reference data and genome indices
2. **FASTQC** - Quality control on raw reads (runs if alignment stage enabled)
3. **READ_ALIGNMENT_RNA** - STAR alignment → SAMtools sort → Sambamba merge → GATK MarkDuplicates
4. **RSEQC_ANALYSIS** - RSeQC quality metrics including rRNA contamination check via split_bam
5. **rRNA QC Gate** - Filters samples exceeding rRNA thresholds (configurable via `--rrna_threshold_count`, `--rrna_threshold_percent`)
6. **ISOFOX_QUANTIFICATION** - Transcript quantification and fusion detection (only for samples passing rRNA QC)
7. **MULTIQC** - Aggregate QC reports

The pipeline uses Nextflow channels to pass data between stages. Key channel pattern: `ch_inputs` is parsed from the samplesheet in `Utils.parseInput()`, then flows through alignment → RSeQC → Isofox with BAM/BAI tuples.

### Key Subworkflows (subworkflows/local/)

- `prepare_reference/` - Prepares/downloads reference data (genome indices, HMF data)
- `read_alignment_rna/` - STAR alignment → SAMtools sort → Sambamba merge → GATK MarkDuplicates
- `isofox_quantification/` - Isofox transcript quantification and fusion calling
- `rseqc_analysis/` - RSeQC quality control metrics

### Library Classes (lib/)

Groovy classes providing utilities:
- `Constants.groovy` - Enums (RunMode, Process, FileType, SampleType), HMF reference data URLs, input type mappings
- `Utils.groovy` - Input CSV parsing, sample extraction (`getTumorRnaSample`, `getTumorRnaBam`, etc.), rRNA QC logic
- `WorkflowMain.groovy` - Parameter validation, run config generation
- `Processes.groovy` - Process stage configuration, determines which stages run based on inputs

### Configuration Files (conf/)

- `base.config` - Default resource allocations
- `hmf_genomes.config` - HMF genome reference paths (STAR index, BED files, etc.)
- `hmf_data.config` - HMF tool reference data paths (Isofox counts, known fusions, etc.)
- `modules.config` - Module-specific ext.args and publish settings
- `test.config` / `test_full.config` / `test_stub.config` - Test profile settings

### Module Organization

- `modules/nf-core/` - Standard nf-core modules (FastQC, MultiQC, STAR, GATK4, SAMtools, RSeQC, etc.)
- `modules/local/` - Custom modules (fastp, isofox, gridss_index, star_align, rseqc_splitbam)

## Input Samplesheet Format

CSV with columns: `group_id`, `subject_id`, `sample_id`, `sample_type`, `sequence_type`, `filetype`, `filepath`, `info`

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath,info
SAMPLE1,SUBJECT1,SAMPLE1_T,tumor,rna,fastq,/path/to/R1.fq.gz;/path/to/R2.fq.gz,library_id:LIB001;lane:L001
```

- `sample_type`: `tumor`
- `sequence_type`: `rna`
- `filetype`: `fastq`, `bam`, `bai`, `isofox_dir`, `rseqc_dir`
- `info`: Semicolon-separated key:value pairs (required for FASTQ: `library_id`, `lane`)
- For FASTQs, `filepath` contains forward and reverse reads separated by `;`

## Key Parameters

- `--genome`: Reference genome (`GRCh38_hmf` supported)
- `--input`: Samplesheet CSV path
- `--mode`: Workflow mode (`rna_workflow`, `fastqc_workflow`)
- `--isofox_functions`: Semicolon-separated Isofox analysis functions (default: `TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS;RETAINED_INTRONS`)
- `--max_fastq_records`: FASTQ splitting threshold (0 = no split)
- `--rrna_threshold_count`: Maximum rRNA read count before failing QC (default: 143303744)
- `--rrna_threshold_percent`: Maximum rRNA percentage before failing QC (0 = disabled)
- `--prepare_reference_only`: Only stage reference data, don't run analysis

### Process Control

Control which pipeline stages run using `--processes_include` and `--processes_exclude` (comma-separated):
- `alignment` - STAR alignment and BAM processing
- `isofox` - Isofox transcript quantification
- `rseqc` - RSeQC quality metrics

```bash
# Run only alignment and RSeQC (skip Isofox)
nextflow run . -profile docker --input samplesheet.csv --outdir results --processes_exclude isofox

# Run only Isofox on pre-aligned BAMs
nextflow run . -profile docker --input samplesheet_with_bams.csv --outdir results --processes_include isofox
```

## Reference Data

HMF reference data is automatically downloaded/extracted from public URLs defined in `Constants.groovy`. The pipeline uses HMF pipeline resources v2.3.0 for GRCh38. Use `--prepare_reference_only true` to pre-stage reference data.

## Debugging

```bash
# Enable debug profile for verbose output
nextflow run . -profile test,docker,debug --outdir results

# Check Nextflow logs
cat .nextflow.log

# View execution reports (generated in outdir/pipeline_info/)
# - execution_timeline_*.html - Timeline visualization
# - execution_report_*.html - Resource usage report
# - execution_trace_*.txt - Detailed process trace
# - pipeline_dag_*.html - DAG visualization
```
