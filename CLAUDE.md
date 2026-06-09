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
# Newer QC/filtering modules also have nf-test coverage:
cd modules/nf-core/rustqc && nf-test test tests/main.nf.test
cd modules/nf-core/ribodetector && nf-test test tests/main.nf.test
cd modules/nf-core/sortmerna && nf-test test tests/main.nf.test
```

## Architecture

### Run Modes

The pipeline supports multiple run modes controlled by `--mode` parameter (default: `rna_workflow`).

**Current state (2026-06-04):** `Constants.RunMode` defines five modes, dispatched in `main.nf` → `NFCORE_HMFRNASEQ`:
- `rna_standard` (**DEFAULT**): No rRNA pre-filter → STAR (2.7.3a) → GATK MarkDuplicates → then **RustQC, RSeQC, and Isofox fan out in parallel** (no rRNA gate; Isofox runs on all samples). RustQC and RSeQC both count rRNA post-alignment and are independently disablable via `--skip_rustqc` / `--skip_rseqc`. STAR is pinned to 2.7.3a; the `star_index` default is the HMF 2.7.3a tarball, overridable at runtime with `--ref_data_genome_star_index` (a tarball, an uncompressed index directory, or unset to build from fasta+GTF). Isofox uses the coordinate-compatible HMF fasta.
- `rna_workflow`: Full RNA analysis with **SortMeRNA** rRNA pre-filtering → STAR → GATK MarkDuplicates → RustQC → rRNA gate → Isofox
- `rna_redux_workflow`: Alternative pipeline using **Ribodetector** rRNA pre-filtering → STAR → **REDUX** (fixmate + duplicate marking) → optional **AMBER** BAF profiling → RustQC → Isofox
- `fastqc_workflow`: FastQC only
- `redux_test`: REDUX processing test harness

> Historically (earlier doc): only `rna_workflow` and `fastqc_workflow` existed; `rna_workflow` was the default. As of 2026-06-04 the default is `rna_standard`.

### Entry Points

- `main.nf` - Pipeline entry point, dispatches to workflows based on mode
- `workflows/rna_standard_workflow.nf` - **Default** workflow (no pre-filter; STAR 2.7.3a + GATK MarkDuplicates; RustQC/RSeQC/Isofox in parallel, no gate)
- `workflows/rna_workflow.nf` - Main RNA analysis workflow (SortMeRNA + GATK MarkDuplicates path)
- `workflows/rna_redux_workflow.nf` - Alternative RNA workflow (Ribodetector + REDUX + AMBER path)
- `workflows/fastqc_workflow.nf` - FastQC-only workflow
- `workflows/redux_test.nf` - REDUX test workflow
- `workflows/hmfrnaseq.nf` - Legacy/scaffold workflow (currently not wired into `main.nf`)

### Pipeline Stages (rna_workflow)

**Current state (2026-06-04)** — `rna_workflow.nf`:

1. **PREPARE_REFERENCE** - Downloads/stages HMF reference data and genome indices (including the SortMeRNA rRNA database)
2. **FASTQC** - Quality control on raw reads (runs if alignment stage enabled)
3. **SORTMERNA_FILTER** - rRNA read filtering on FASTQs before alignment (skippable via `--skip_sortmerna`)
4. **READ_ALIGNMENT_RNA** - STAR alignment → SAMtools sort → Sambamba merge → GATK MarkDuplicates
5. **RUSTQC_ANALYSIS** - Single-pass RNA QC via RustQC, with GTF-based **biotype counts** used for rRNA detection (replaces the old RSeQC `split_bam` approach). Requires `--ref_data_genome_gtf`.
6. **RRNA_QC_GATE** - Parses RustQC biotype counts and gates samples exceeding rRNA thresholds (`--rrna_threshold_count`, `--rrna_threshold_percent`); only passing samples proceed to Isofox
7. **ISOFOX_QUANTIFICATION** - Transcript quantification and fusion detection (only for samples passing the rRNA gate)
8. **MULTIQC** - Two reports are produced: a **per-sample** report (`MULTIQC`) grouped by `group_id`, and an **aggregated** report (`MULTIQC_AGGREGATED`, one row per sample, no FastQC)

> Historically (earlier doc): rRNA contamination was checked by RSeQC `split_bam` in an `RSEQC_ANALYSIS` stage, and only a single MultiQC report was produced. RSeQC has been replaced by RustQC; the `rseqc` process flag (see Process Control) is retained as the gate that now enables `RUSTQC_ANALYSIS`.

**Important gotcha:** the process/stage flag is still named `rseqc` (in `Constants.Process` and `run_config.stages.rseqc`), but it now controls **RustQC**, not RSeQC. Use `--processes_exclude rseqc` to skip RustQC.

### Pipeline Stages (rna_redux_workflow)

`rna_redux_workflow.nf` is the alternative path. Differences from `rna_workflow`:
- rRNA pre-filtering uses **Ribodetector** (`RIBODETECTOR_FILTER`, skippable via `--skip_rrna_filter`) instead of SortMeRNA
- Alignment uses `READ_ALIGNMENT_RNA_STAR` (STAR + sort + merge, **no** GATK MarkDuplicates)
- **REDUX_PROCESSING** handles fixmate + duplicate marking/unmapping for all BAMs (including pre-aligned input)
- Optional **AMBER** stage runs germline-only BAF profiling on the RNA BAMs/CRAMs
- MultiQC uses a dynamically generated `replace_names.tsv` to normalize sample names across tools (Qualimap embeds the BAM filename; FastQC uses per-lane names) to the canonical `<group_id>_<sample_id>`

The pipeline uses Nextflow channels to pass data between stages. Key channel pattern: `ch_inputs` is parsed from the samplesheet in `Utils.parseInput()`, then flows through alignment → RustQC → rRNA gate → Isofox with `[meta, bam, bai]` tuples. Pre-aligned BAM/CRAM inputs bypass alignment via `Utils.hasExistingInput(...)`. Because queue channels can only be consumed once, FastQC and alignment each build their own `channel.fromList(inputs)` rather than sharing `ch_inputs`.

### MultiQC reporting — sample naming & metric selection

MultiQC merges metrics by sample name, so tools that emit a non-canonical name show up as **extra rows** in the general-stats table. Two deliberate choices keep the report to **one row per sample** (canonical `<group_id>_<sample_id>`):

1. **`replace_names.tsv` (the `--replace-names` slot of the MULTIQC module) + `sample_names_replace_complete: true`** in `assets/multiqc_config.yml`. Each workflow generates this TSV from the samplesheet to remap tool-specific names → canonical:
   - `rna_redux_workflow`: Qualimap embeds the REDUX BAM filename (`<sample_id>.redux`); FastQC uses per-lane names.
   - `rna_standard_workflow`: only **FastQC** per-lane/per-read names (`<group_id>_<sample_id>_<lib>_<lane>_{1,2}`) need remapping (RustQC/RSeQC/Qualimap/Samtools already report canonical).
2. **MarkDuplicates (GATK4/Picard) metrics are excluded from MultiQC in `rna_standard`.** MultiQC's Picard module names samples from the metrics' LIBRARY field (the read-group `sample_id.library.lane`), which created a duplicate per-sample row. The metric is also **redundant** for an RNA report — duplication is covered by **dupRadar** (RNA dup-vs-expression) + **Samtools** flagstat, and library complexity by **Preseq**. So `ch_markdups_metrics` is not mixed into either MULTIQC channel (the `.markdup.bam.metrics` file is still published to `outdir`; duplicates are still marked in the BAM). Renaming via `replace_names` would also work, but dropping the redundant metric is cleaner.

> Note: this Picard-naming quirk has existed in every `rna_workflow` variant (they all mix `markdups_metrics` into MultiQC without `replace_names`) but was rarely seen because real `rna_workflow` + aggregated-MultiQC runs were uncommon. `rna_standard` addresses it by exclusion.

### Key Subworkflows (subworkflows/local/)

**Current state (2026-06-04):**
- `prepare_reference/` - Prepares/downloads reference data (genome indices, HMF data, SortMeRNA DB)
- `prepare_inputs/` - Input preparation helpers
- `read_alignment_rna/` - STAR alignment → SAMtools sort → Sambamba merge → GATK MarkDuplicates (used by `rna_workflow`)
- `read_alignment_rna_star/` - STAR alignment → sort → merge, no MarkDuplicates (used by `rna_redux_workflow`)
- `redux_processing/` - REDUX fixmate + duplicate marking/unmapping
- `sortmerna_filter/` - SortMeRNA rRNA read filtering (rna_workflow)
- `ribodetector_filter/` - Ribodetector rRNA read filtering (rna_redux_workflow)
- `rustqc_analysis/` - RustQC single-pass RNA QC + biotype counts (current QC stage)
- `rseqc_analysis/` - RSeQC quality control metrics (legacy; superseded by `rustqc_analysis`)
- `rrna_qc_gate/` - Parses biotype counts and gates samples for Isofox
- `isofox_quantification/` - Isofox transcript quantification and fusion calling

> Historically (earlier doc): only `prepare_reference`, `read_alignment_rna`, `isofox_quantification`, and `rseqc_analysis` were listed.

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
- **Current state (2026-06-04):** additional test profiles `test_full_stub.config`, `test_multi_stub.config` (multi-sample stub), `test_redux.config` (REDUX workflow), plus `panel_data.config`, `targeted_parameters.config`, and `igenomes.config` / `igenomes_ignored.config`

### Module Organization

- `modules/nf-core/` - Standard nf-core modules (FastQC, MultiQC, STAR, GATK4, SAMtools, RSeQC, SortMeRNA, Ribodetector, RustQC, bwamem2, etc.)
- `modules/local/` - Custom modules. **Current state (2026-06-04):** `amber`, `bamchecker`, `bwa-mem2`, `custom` (write_reference_data, extract_tarball), `fastp`, `gatk4`, `gridss`, `isofox`, `multiqc`, `redux`, `rseqc_splitbam`, `sambamba`, `samtools` (sort, fixmate_sort), `star`.

> Historically (earlier doc): local modules were listed as fastp, isofox, gridss_index, star_align, rseqc_splitbam.

## Input Samplesheet Format

CSV with columns: `group_id`, `subject_id`, `sample_id`, `sample_type`, `sequence_type`, `filetype`, `filepath`, `info`

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath,info
SAMPLE1,SUBJECT1,SAMPLE1_T,tumor,rna,fastq,/path/to/R1.fq.gz;/path/to/R2.fq.gz,library_id:LIB001;lane:L001
```

- `sample_type`: `tumor`
- `sequence_type`: `rna`
- `filetype`: `fastq`, `bam`, `bai`, `cram`, `crai`, `isofox_dir`, `rseqc_dir` (pre-aligned BAM/CRAM input bypasses the alignment stage)
- `info`: Semicolon-separated key:value pairs (required for FASTQ: `library_id`, `lane`)
- For FASTQs, `filepath` contains forward and reverse reads separated by `;`

## Key Parameters

- `--genome`: Reference genome (`GRCh38_hmf` supported)
- `--input`: Samplesheet CSV path
- `--mode`: Workflow mode (`rna_workflow`, `rna_redux_workflow`, `fastqc_workflow`, `redux_test`)
- `--isofox_functions`: Semicolon-separated Isofox analysis functions (default: `TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS;RETAINED_INTRONS`)
- `--max_fastq_records`: FASTQ splitting threshold (0 = no split)
- `--rrna_threshold_count`: Maximum rRNA read count before failing QC (default: 143303744)
- `--rrna_threshold_percent`: Maximum rRNA percentage before failing QC (0 = disabled)
- `--prepare_reference_only`: Only stage reference data, don't run analysis
- `--ref_data_genome_gtf`: GTF used by RustQC for biotype-based rRNA detection — required when the `rseqc`/RustQC stage runs

**Current state (2026-06-04)** — rRNA filtering and QC parameters added since the earlier doc:
- `--strandedness`: Library strandedness for RustQC (default: `forward`)
- `--skip_sortmerna`: Skip SortMeRNA rRNA pre-filtering in `rna_workflow` (default: false)
- `--sortmerna_fastas`: SortMeRNA rRNA database URL/path
- `--skip_rrna_filter`: Skip Ribodetector rRNA pre-filtering in `rna_redux_workflow` (default: false)
- `--ribodetector_read_length`: Read length passed to Ribodetector (default: 151)

### Process Control

Control which pipeline stages run using `--processes_include` and `--processes_exclude` (comma-separated). Valid stage names come from `Constants.Process`.

**Current state (2026-06-04):** `alignment`, `amber`, `isofox`, `rseqc`, `multiqc`.
- `alignment` - STAR alignment and BAM processing
- `amber` - AMBER germline-only BAF profiling (`rna_redux_workflow` only)
- `isofox` - Isofox transcript quantification
- `rseqc` - **now gates RustQC** (the QC + rRNA-biotype stage), despite the legacy name
- `multiqc` - MultiQC report generation

> Historically (earlier doc): only `alignment`, `isofox`, and `rseqc` were documented, and `rseqc` referred to actual RSeQC.

```bash
# Run only alignment and QC (skip Isofox)
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
