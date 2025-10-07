#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

option_list <- list(
  make_option(c("--params"), type = "character", default = "config/params.yaml"),
  make_option(c("--sample-sheet"), type = "character", help = "Sample metadata TSV"),
  make_option(c("--quant-dir"), type = "character", help = "Directory containing Salmon quant.sf results"),
  make_option(c("--tx2gene"), type = "character", help = "Transcript-to-gene mapping TSV"),
  make_option(c("--out-counts"), type = "character", help = "Gene-level counts TSV output"),
  make_option(c("--out-tpm"), type = "character", help = "Gene-level TPM TSV output"),
  make_option(c("--out-txi"), type = "character", help = "tximport RDS output path"),
  make_option(c("--out-dds"), type = "character", help = "DESeq2 dataset RDS output"),
  make_option(c("--out-summary"), type = "character", help = "Summary TSV output"),
  make_option(c("--contrast-file"), type = "character", help = "Contrasts TSV"),
  make_option(c("--figdir"), type = "character", help = "Directory for DE plots", default = "results/de/figures"),
  make_option(c("--se"), action = "store_true", default = FALSE, help = "Flag for single-end reads")
)

opt <- parse_args(OptionParser(option_list = option_list))
params <- read_yaml(opt$params)
set.seed(42)

required_args <- c("sample_sheet", "quant_dir", "tx2gene", "out_counts", "out_tpm", "out_txi", "out_dds", "out_summary", "contrast_file")
for (arg in required_args) {
  if (is.null(opt[[arg]]) || is.na(opt[[arg]])) {
    stop(sprintf("Missing required argument: %s", arg), call. = FALSE)
  }
}

for (path in c(opt$`out_counts`, opt$`out_tpm`, opt$`out_txi`, opt$`out_dds`, opt$`out_summary`)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

samples <- readr::read_tsv(opt$`sample_sheet`, show_col_types = FALSE)
if (!"sample" %in% names(samples)) {
  stop("Sample sheet must contain a 'sample' column", call. = FALSE)
}
if (!"condition" %in% names(samples)) {
  stop("Sample sheet must contain a 'condition' column", call. = FALSE)
}

samples <- samples %>% mutate(sample = as.character(sample))
rownames(samples) <- samples$sample

quant_paths <- file.path(opt$`quant_dir`, samples$sample, "quant.sf")
if (any(!file.exists(quant_paths))) {
  missing <- quant_paths[!file.exists(quant_paths)]
  stop(sprintf("Missing quant.sf files for samples: %s", paste(basename(dirname(missing)), collapse = ", ")), call. = FALSE)
}

files <- setNames(quant_paths, samples$sample)
tx2gene <- readr::read_tsv(opt$tx2gene, col_names = c("transcript", "gene"), show_col_types = FALSE)

message("Running tximport...")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
counts_tbl <- as.data.frame(txi$counts) %>% tibble::rownames_to_column(params$r$gene_id_column %||% "gene_id")
tpm_tbl <- as.data.frame(txi$abundance) %>% tibble::rownames_to_column(params$r$gene_id_column %||% "gene_id")
readr::write_tsv(counts_tbl, opt$`out_counts`)
readr::write_tsv(tpm_tbl, opt$`out_tpm`)

saveRDS(txi, file = opt$`out_txi`)

message("Building DESeq2 dataset...")
design_formula <- stats::as.formula(params$r$design)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design_formula)
if (any(sizeFactors(dds) == 0)) {
  dds <- estimateSizeFactors(dds)
}
dds <- DESeq(dds, quiet = TRUE)

saveRDS(dds, file = opt$`out_dds`)

contrasts <- readr::read_tsv(opt$`contrast_file`, show_col_types = FALSE)
if (!all(c("groupA", "groupB") %in% names(contrasts))) {
  stop("Contrast file must contain 'groupA' and 'groupB' columns", call. = FALSE)
}
contrast_variable <- params$r$contrast_variable %||% "condition"
alpha <- params$r$alpha %||% 0.05
lfc_shrink <- isTRUE(params$r$lfc_shrink)

summary_rows <- list()

for (i in seq_len(nrow(contrasts))) {
  groupA <- contrasts$groupA[i]
  groupB <- contrasts$groupB[i]
  label <- sprintf("%s_vs_%s", groupA, groupB)
  message(sprintf("Running DESeq2 contrast: %s", label))
  res <- results(dds, contrast = c(contrast_variable, groupA, groupB), alpha = alpha, tidy = TRUE)
  res <- res %>% rename(gene_id = row) %>% arrange(padj)
  coef_name <- sprintf("%s_%s_vs_%s", contrast_variable, groupA, groupB)
  if (lfc_shrink && coef_name %in% resultsNames(dds)) {
    suppressWarnings({
      shrink <- try(lfcShrink(dds, coef = coef_name, type = "apeglm"), silent = TRUE)
    })
    if (!inherits(shrink, "try-error")) {
      res$log2FoldChange <- shrink$log2FoldChange
      res$lfcSE <- shrink$lfcSE
    }
  }
  res <- res %>% mutate(comparison = label)
  out_file <- file.path(dirname(opt$`out_summary`), sprintf("DE_%s.tsv", label))
  readr::write_tsv(res, out_file)

  ma_plot <- ggplot(res, aes(x = baseMean, y = log2FoldChange, colour = padj < alpha)) +
    geom_point(alpha = 0.6, size = 1.1) +
    scale_x_log10() +
    scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
    labs(title = paste0("MA plot: ", groupA, " vs ", groupB), x = "Base mean", y = "log2FC", colour = paste0("padj < ", alpha)) +
    theme_minimal()
  ggsave(filename = file.path(opt$figdir, sprintf("MA_%s.png", label)), plot = ma_plot, width = 6, height = 4, dpi = 150)

  volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), colour = padj < alpha)) +
    geom_point(alpha = 0.6, size = 1.1) +
    scale_colour_manual(values = c(`TRUE` = "steelblue", `FALSE` = "grey70")) +
    labs(title = paste0("Volcano: ", groupA, " vs ", groupB), x = "log2FC", y = "-log10(pvalue)", colour = paste0("padj < ", alpha)) +
    theme_minimal()
  ggsave(filename = file.path(opt$figdir, sprintf("Volcano_%s.png", label)), plot = volcano, width = 6, height = 4, dpi = 150)

  up_n <- sum(res$padj < alpha & res$log2FoldChange > 0, na.rm = TRUE)
  down_n <- sum(res$padj < alpha & res$log2FoldChange < 0, na.rm = TRUE)
  summary_rows[[label]] <- tibble::tibble(
    comparison = label,
    groupA = groupA,
    groupB = groupB,
    detected_genes = sum(!is.na(res$padj)),
    significant = sum(res$padj < alpha, na.rm = TRUE),
    up = up_n,
    down = down_n
  )
}

summary_tbl <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else tibble::tibble()
readr::write_tsv(summary_tbl, opt$`out_summary`)

norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% tibble::rownames_to_column(params$r$gene_id_column %||% "gene_id")
readr::write_tsv(norm_counts, file.path(dirname(opt$`out_summary`), "normalized_counts.tsv"))
