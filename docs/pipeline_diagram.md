# nf-core/hmfrnaseq Pipeline Diagram

## Workflow Overview

```mermaid
flowchart TD
    subgraph Input["Input"]
        FASTQ[/"FASTQ files"/]
        BAM_IN[/"Pre-aligned BAM"/]
    end

    subgraph QC["Quality Control"]
        FASTQC["FastQC"]
    end

    subgraph Alignment["Read Alignment"]
        STAR["STAR Align"]
        SORT["SAMtools Sort"]
        MERGE["Sambamba Merge"]
        MARKDUP["GATK MarkDuplicates"]
    end

    subgraph RSeQC["RSeQC Analysis"]
        BAMSTAT["bam_stat.py"]
        READDUP["read_duplication.py"]
        SPLITBAM["split_bam.py"]
    end

    subgraph Gate["rRNA QC Gate"]
        CHECK{"rRNA Check"}
        PASS["Pass"]
        FAIL["Fail - Skip Isofox"]
    end

    subgraph Analysis["Transcript Analysis"]
        ISOFOX["Isofox"]
    end

    subgraph Reports["Reports"]
        MQC_SAMPLE["Per-sample MultiQC"]
        MQC_AGG["Aggregated MultiQC"]
    end

    %% Input paths
    FASTQ --> FASTQC
    FASTQ --> STAR
    BAM_IN --> BAMSTAT
    BAM_IN --> READDUP
    BAM_IN --> SPLITBAM

    %% Alignment path
    STAR --> SORT --> MERGE --> MARKDUP

    %% RSeQC from aligned BAM
    MARKDUP --> BAMSTAT
    MARKDUP --> READDUP
    MARKDUP --> SPLITBAM

    %% rRNA gate
    SPLITBAM --> CHECK
    CHECK --> PASS
    CHECK --> FAIL

    %% Analysis
    PASS --> ISOFOX
    MARKDUP --> ISOFOX

    %% Reports
    FASTQC --> MQC_SAMPLE
    BAMSTAT --> MQC_SAMPLE
    BAMSTAT --> MQC_AGG
    READDUP --> MQC_SAMPLE
    SPLITBAM --> MQC_SAMPLE
    SPLITBAM --> MQC_AGG
    MARKDUP --> MQC_SAMPLE
    MARKDUP --> MQC_AGG
```

## Simplified Linear View

```mermaid
flowchart LR
    A[FASTQ] --> B[FastQC]
    A --> C[STAR]
    C --> D[Sort/Merge]
    D --> E[MarkDuplicates]
    E --> F[RSeQC]
    F --> G{rRNA OK?}
    G -->|Yes| H[Isofox]
    G -->|No| I[Skip]
    H --> J[MultiQC]
    F --> J
    B --> J
```

## Process Labels

| Process | Label | Resources (Attempt 1) |
|---------|-------|----------------------|
| STAR_ALIGN | process_high | 12 CPU, 72GB |
| GATK4_MARKDUPLICATES | process_high | 12 CPU, 72GB |
| ISOFOX | process_high | 12 CPU, 72GB |
| RSEQC_* | process_high | 12 CPU, 72GB |
| FASTQC | process_medium | 6 CPU, 36GB |
| MULTIQC | process_single | 1 CPU, 6GB |
