nextflow.enable.dsl=2

// Import caching system
try {
    evaluate(new File('../../scripts/cache_manager.py'))
    CACHE_ENABLED = params.cache?.enabled ?: true
    CACHE_DIR = params.cache?.dir ?: '.cache'
    cache_manager = get_cache_manager(CACHE_DIR, CACHE_ENABLED)
} catch (Exception e) {
    CACHE_ENABLED = false
    cache_manager = null
}

include { FASTQC; MULTIQC } from './modules/qc.nf'
include { SALMON_INDEX; SALMON_QUANT } from './modules/salmon.nf'
include { TXIMPORT_DESEQ2; FGSEA } from './modules/tximport_de.nf'
include { SINGLECELL_ANALYSIS } from './modules/singlecell.nf'
include { REPORT } from './modules/report.nf'

params.project = params.project ?: 'rnaseq-mini'
params.paths = params.paths ?: [:]
params.paths.samples = params.paths.samples ?: 'config/samples.tsv'
params.paths.salmon = params.paths.salmon ?: 'results/salmon'
params.paths.qc = params.paths.qc ?: 'results/qc'
params.paths.outdir = params.paths.outdir ?: 'results'
params.reference = params.reference ?: [:]
params.reference.transcripts_fa = params.reference.transcripts_fa ?: 'references/yeast/transcripts.fa.gz'
params.reference.annotation_gtf = params.reference.annotation_gtf ?: 'references/yeast/annotation.gtf.gz'
params.reference.decoy_fasta = params.reference.decoy_fasta ?: null
params.reference.salmon_index = params.reference.salmon_index ?: 'references/yeast/salmon_index'
params.r = params.r ?: [:]
params.r.contrasts_file = params.r.contrasts_file ?: 'config/contrasts.tsv'
params.se = params.se ?: false

// Single-cell parameters
params.singlecell = params.singlecell ?: [:]
params.singlecell.enabled = params.singlecell.enabled ?: false
params.singlecell.technology = params.singlecell.technology ?: 'auto'
params.singlecell.method = params.singlecell.method ?: 'auto'
params.singlecell.chemistry = params.singlecell.chemistry ?: 'auto'
params.singlecell.min_genes = params.singlecell.min_genes ?: 200
params.singlecell.max_genes = params.singlecell.max_genes ?: 6000
params.singlecell.mito_threshold = params.singlecell.mito_threshold ?: 10.0

workflow {
    samples_file = file(params.paths.samples)
    contrasts_file = file(params.r.contrasts_file)

    samples_ch = Channel.fromPath(samples_file)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def fastq1 = file(row.fastq_1)
            def fastq2 = row.fastq_2 ? file(row.fastq_2) : null
            def meta = row + [fastq_1: fastq1, fastq_2: fastq2]
            tuple(row.sample, meta)
        }
        .share()

    fastqc_input = samples_ch.flatMap { sample, meta ->
        def entries = [ tuple(sample, 'R1', meta.fastq_1) ]
        if (!params.se && meta.fastq_2) {
            entries << tuple(sample, 'R2', meta.fastq_2)
        }
        entries
    }

    fastqc_results = FASTQC(fastqc_input)
    multiqc_trigger = fastqc_results.collect()
    multiqc = MULTIQC(multiqc_trigger)

    index_input = Channel.of(tuple(file(params.reference.transcripts_fa), file(params.reference.annotation_gtf), params.reference.decoy_fasta ? file(params.reference.decoy_fasta) : null))
    salmon_index = SALMON_INDEX(index_input)

    quant_input = samples_ch.map { sample, meta -> tuple(sample, meta, meta.fastq_1, meta.fastq_2) }
    quant_results = SALMON_QUANT(quant_input, salmon_index.out.index)

    quant_trigger = quant_results.collect()
    sample_sheet_ch = Channel.value(samples_file)
    contrast_file_ch = Channel.value(contrasts_file)
    tximport = TXIMPORT_DESEQ2(quant_trigger, sample_sheet_ch, contrast_file_ch, salmon_index.out.tx2gene)

    fgsea = FGSEA(tximport.out.de_dir, contrast_file_ch)

    // Single-cell analysis (if enabled)
    if (params.singlecell.enabled) {
        // Prepare transcript-to-gene mapping for kallisto|bustools
        t2g_file = file("${params.reference.transcripts_fa}".replace('.fa.gz', '_t2g.txt'))

        singlecell_results = SINGLECELL_ANALYSIS(samples_ch, params.reference.salmon_index, t2g_file)

        // Include single-cell results in report
        report = REPORT(multiqc.out, tximport.out.counts_dir, tximport.out.de_dir,
                       fgsea.out.fgsea_dir, singlecell_results.qc_reports.collect(),
                       singlecell_results.clustering_results.collect())
    } else {
        report = REPORT(multiqc.out, tximport.out.counts_dir, tximport.out.de_dir, fgsea.out.fgsea_dir)
    }
}
