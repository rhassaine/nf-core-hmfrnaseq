# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

nf-core/hmfrnaseq is a Nextflow DSL2 bioinformatics pipeline for RNA-seq analysis, built on the nf-core framework. It integrates Hartwig Medical Foundation (HMF) tools for RNA analysis including transcript quantification and fusion detection.

## Common Commands

### Running the Pipeline

```bash
# Run with test profile
nextflow run . -profile test,docker --outdir results

# Run with singularity
nextflow run . -profile test,singularity --outdir results

# Full run with custom input
nextflow run . -profile docker --input samplesheet.csv --outdir results --genome GRCh38_hmf
```

### Linting and Validation

```bash
# Validate pipeline schema
nextflow run . --help

# Run nf-core linting (requires nf-core tools)
nf-core lint
```

### Testing

The pipeline uses nf-test for module and subworkflow testing. Test files are located alongside modules in `tests/` subdirectories with `.nf.test` extension.

```bash
# Run nf-test (from module/subworkflow directory)
nf-test test
```

## Architecture

### Run Modes

The pipeline supports multiple run modes controlled by `--mode` parameter (default: `rna_workflow`):
- `rna_workflow`: Full RNA analysis (alignment, QC, Isofox quantification)
- `fastqc_workflow`: FastQC only

### Entry Points

- [main.nf](main.nf) - Pipeline entry point, dispatches to workflows based on mode
- [workflows/rna_workflow.nf](workflows/rna_workflow.nf) - Main RNA analysis workflow
- [workflows/fastqc_workflow.nf](workflows/fastqc_workflow.nf) - FastQC-only workflow

### Key Subworkflows (subworkflows/local/)

- `prepare_reference/` - Prepares/downloads reference data (genome indices, HMF data)
- `read_alignment_rna/` - STAR alignment → SAMtools sort → Sambamba merge → GATK MarkDuplicates
- `isofox_quantification/` - Isofox transcript quantification and fusion calling
- `rseqc_analysis/` - RSeQC quality control metrics

### Library Classes (lib/)

Groovy classes providing utilities:
- [Constants.groovy](lib/Constants.groovy) - Enums (RunMode, Process, FileType, SampleType), reference data URLs, input type mappings
- [Utils.groovy](lib/Utils.groovy) - Input parsing, sample extraction utilities
- [WorkflowMain.groovy](lib/WorkflowMain.groovy) - Parameter validation, run config generation
- [Processes.groovy](lib/Processes.groovy) - Process stage configuration

### Configuration Files (conf/)

- `base.config` - Default resource allocations
- `hmf_genomes.config` - HMF genome reference paths
- `hmf_data.config` - HMF tool reference data paths
- `modules.config` - Module-specific ext.args and publish settings
- `test.config` - Test profile settings

### Module Organization

- `modules/nf-core/` - Standard nf-core modules (FastQC, MultiQC, STAR, GATK4, etc.)
- `modules/local/` - Custom modules (fastp, isofox, gridss index, star align wrapper)

## Key Parameters

- `--genome`: Reference genome (GRCh37_hmf, GRCh38_hmf)
- `--input`: Samplesheet CSV path
- `--mode`: Workflow mode (rna_workflow, fastqc_workflow)
- `--isofox_functions`: Comma-separated Isofox analysis functions
- `--max_fastq_records`: FASTQ splitting threshold (0 = no split)

## Input Samplesheet Format

CSV with columns defining samples, FASTQs, and optional pre-computed results. Sample types include tumor RNA with paired-end FASTQs.

## Reference Data

HMF reference data is automatically downloaded/extracted from public URLs defined in Constants.groovy. Supports both GRCh37 and GRCh38. Use `--prepare_reference_only true` to pre-stage reference data.
