# nf-core/hmfrnaseq: Citations

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

> Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

> Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

- [STAR](https://pubmed.ncbi.nlm.nih.gov/23104886/)

> Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. PubMed PMID: 23104886; PubMed Central PMCID: PMC3530905.

- [SAMtools](https://pubmed.ncbi.nlm.nih.gov/19505943/)

> Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.

- [Sambamba](https://pubmed.ncbi.nlm.nih.gov/25697820/)

> Tarasov A, Vilella AJ, Cuppen E, Nijman IJ, Prins P. Sambamba: fast processing of NGS alignment formats. Bioinformatics. 2015 Jun 15;31(12):2032-4. doi: 10.1093/bioinformatics/btv098. PubMed PMID: 25697820; PubMed Central PMCID: PMC4765878.

- [GATK4 / Picard MarkDuplicates](https://gatk.broadinstitute.org/) — duplicate marking on RNA alignments.

> Van der Auwera GA, O'Connor BD. Genomics in the Cloud: Using Docker, GATK, and WDL in Terra. 1st ed. O'Reilly Media; 2020. ("Picard Toolkit", Broad Institute, GitHub repository: https://broadinstitute.github.io/picard/)

- [RustQC](https://github.com/seqeralabs/rustqc) — single-pass post-alignment RNA-seq QC (replaces dupRadar, featureCounts biotype QC, RSeQC, Preseq, Qualimap, SAMtools stats); used here for biotype-based rRNA counting.

> Seqera Labs. RustQC: high-performance RNA-seq QC. https://seqeralabs.github.io/RustQC/

- [RSeQC](https://pubmed.ncbi.nlm.nih.gov/22743226/) — independent rRNA quantification (split_bam) used as a cross-check against RustQC.

> Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments. Bioinformatics. 2012 Aug 15;28(16):2184-5. doi: 10.1093/bioinformatics/bts356. PubMed PMID: 22743226.

- [Isofox](https://github.com/hartwigmedical/hmftools/tree/master/isofox) — RNA transcript quantification, novel splice junctions, and fusion calling (Hartwig Medical Foundation hmftools).

> Hartwig Medical Foundation. hmftools: Isofox. https://github.com/hartwigmedical/hmftools

- [SortMeRNA](https://pubmed.ncbi.nlm.nih.gov/23071270/) — pre-alignment rRNA filtering (`rna_workflow`).

> Kopylova E, Noé L, Touzet H. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics. 2012 Dec 15;28(24):3211-7. doi: 10.1093/bioinformatics/bts611. PubMed PMID: 23071270.

- [RiboDetector](https://pubmed.ncbi.nlm.nih.gov/35631930/) — pre-alignment rRNA filtering (`rna_redux_workflow`).

> Deng ZL, Münch PC, Mreches R, McHardy AC. Rapid and accurate identification of ribosomal RNA sequences via deep learning. Nucleic Acids Res. 2022 Jun 24;50(11):e60. doi: 10.1093/nar/gkac112. PubMed PMID: 35631930; PubMed Central PMCID: PMC9226531.

- [REDUX](https://github.com/hartwigmedical/hmftools/tree/master/redux) and [AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber) — duplicate marking/unmapping and BAF profiling (`rna_redux_workflow`; Hartwig Medical Foundation hmftools).

> Hartwig Medical Foundation. hmftools. https://github.com/hartwigmedical/hmftools

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
