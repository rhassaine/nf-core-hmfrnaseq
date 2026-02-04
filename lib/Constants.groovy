class Constants {

    static List GENOMES_VERSION_38 = ['GRCh38_hmf', 'GRCh38', 'hg38']
    static List GENOMES_SUPPORTED  = ['GRCh38_hmf']
    static List GENOMES_DEFINED    = Constants.GENOMES_VERSION_38

    // HMF reference data URL (from oncoanalyser v2.3.0)
    static String HMF_DATA_38_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.3.0--2.tar.gz'

    static Integer DEFAULT_ISOFOX_READ_LENGTH = 151

    static enum RunMode {
        RNA_WORKFLOW,
        RNA_REDUX_WORKFLOW,
        FASTQC_WORKFLOW,
        REDUX_TEST,
    }

    static enum Process {
        ALIGNMENT,
        ISOFOX,
        RSEQC,
        MULTIQC,
    }

    static enum FileType {
        // Generic
        BAM,
        BAI,
        FASTQ,
        // Process outputs
        ISOFOX_DIR,
        RSEQC_DIR,
    }

    static enum SampleType {
        TUMOR,
    }

    static enum SequenceType {
        RNA,
    }

    static enum InfoField {
        LANE,
        LIBRARY_ID,
    }

    static Map PLACEHOLDER_META = [meta_placeholder: null]
    static List PLACEHOLDER_OPTIONAL_CHANNEL = []

    static Map INPUT = [

        // RNA BAMs
        BAM_RNA_TUMOR: [
            FileType.BAM,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],

        BAI_RNA_TUMOR: [
            FileType.BAI,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],

        // Process outputs
        ISOFOX_DIR: [
            FileType.ISOFOX_DIR,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],

        RSEQC_DIR: [
            FileType.RSEQC_DIR,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],

    ]
}
