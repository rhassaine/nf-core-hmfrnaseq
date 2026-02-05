# nf-core/hmfrnaseq

[![GitHub Actions CI Status](https://github.com/nf-core/hmfrnaseq/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/hmfrnaseq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/hmfrnaseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/hmfrnaseq/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/hmfrnaseq)

## Introduction

**nf-core/hmfrnaseq** is a bioinformatics pipeline for RNA-seq analysis integrating WiGiTS tools developed by the Hartwig Medical Foundation. It accepts paired-end Illumina FASTQ files or pre-aligned BAM files as input. The pipeline performs read alignment, quality control (including rRNA contamination checks), and uses Isofox for transcript quantification, alternative splicing detection, and fusion calling. Results are summarised in per-sample and aggregated MultiQC reports.

```mermaid
flowchart LR
    subgraph Input
        FASTQ[FASTQ files]
    end

    subgraph QC
        FASTQC[FastQC]
    end

    subgraph Alignment
        STAR[STAR] --> SORT[SAMtools Sort] --> MERGE[Sambamba Merge] --> MARKDUP[MarkDuplicates]
    end

    subgraph "RNA QC"
        BAMSTAT[RSeQC BamStat]
        READDUP[RSeQC ReadDup]
        SPLITBAM[RSeQC SplitBAM]
    end

    subgraph "rRNA Gate"
        GATE{Pass/Fail}
    end

    subgraph Analysis
        ISOFOX[Isofox]
    end

    subgraph Reports
        MULTIQC[MultiQC]
    end

    FASTQ --> FASTQC
    FASTQ --> STAR
    MARKDUP --> BAMSTAT
    MARKDUP --> READDUP
    MARKDUP --> SPLITBAM
    SPLITBAM --> GATE
    GATE -->|Pass| ISOFOX
    BAMSTAT --> MULTIQC
    READDUP --> MULTIQC
    SPLITBAM --> MULTIQC
    FASTQC --> MULTIQC
    ISOFOX --> MULTIQC
```

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Alignment ([`STAR`](https://github.com/alexdobin/STAR))
3. BAM processing ([`SAMtools`](http://www.htslib.org/), [`Sambamba`](https://lomereiter.github.io/sambamba/))
4. Duplicate marking ([`GATK MarkDuplicates`](https://gatk.broadinstitute.org/))
5. RNA QC metrics ([`RSeQC`](http://rseqc.sourceforge.net/))
6. rRNA contamination check and filtering
7. Transcript quantification and fusion detection ([`Isofox`](https://github.com/hartwigmedical/hmftools/tree/master/isofox))
8. QC report ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
SAMPLE1,SUBJECT1,SAMPLE1_T,tumor,rna,fastq,library_id:LIB001;lane:L001,/path/to/R1.fastq.gz;/path/to/R2.fastq.gz
```

Each row represents a pair of FASTQ files for one lane. For samples sequenced across multiple lanes, add one row per lane with the same identifiers.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/hmfrnaseq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --genome GRCh38_hmf
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

nf-core/hmfrnaseq was originally written by Rayan Hassaine.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/hmfrnaseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
