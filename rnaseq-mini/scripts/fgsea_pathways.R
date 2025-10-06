#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(readr)
  library(dplyr)
  library(fgsea)
  library(ggplot2)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

option_list <- list(
  make_option(c("--params"), type = "character", default = "config/params.yaml"),
  make_option(c("--de-dir"), type = "character", help = "Directory containing DESeq2 results"),
  make_option(c("--contrast-file"), type = "character", help = "Contrasts TSV"),
  make_option(c("--outdir"), type = "character", help = "Output directory"),
  make_option(c("--figdir"), type = "character", help = "Directory for enrichment plots", default = "results/fgsea/figures")
)

opt <- parse_args(OptionParser(option_list = option_list))
params <- read_yaml(opt$params)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)

fgsea_params <- params$fgsea %||% list()
score_column <- fgsea_params$score_column %||% "log2FoldChange"
min_size <- fgsea_params$min_size %||% 15
max_size <- fgsea_params$max_size %||% 500
nperm <- fgsea_params$nperm %||% 1000
padj_cutoff <- fgsea_params$padj_cutoff %||% 0.1

load_pathways <- function(path) {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) {
    return(list(
      Example_Positive = c("YAL001C", "YBR002W"),
      Example_Negative = c("YAL001C")
    ))
  }
  ext <- tools::file_ext(path)
  if (tolower(ext) == "gmt") {
    return(gmtPathways(path))
  }
  if (tolower(ext) == "rds") {
    return(readRDS(path))
  }
  if (tolower(ext) == "tsv" || tolower(ext) == "csv") {
    tbl <- readr::read_delim(path, delim = ifelse(grepl("csv$", path), ",", "\t"), show_col_types = FALSE)
    stopifnot(all(c("pathway", "gene") %in% names(tbl)))
    return(split(tbl$gene, tbl$pathway))
  }
  stop(sprintf("Unsupported pathway file format: %s", path), call. = FALSE)
}

pathways <- load_pathways(fgsea_params$genesets)

contrasts <- readr::read_tsv(opt$`contrast_file`, show_col_types = FALSE)
if (!all(c("groupA", "groupB") %in% names(contrasts))) {
  stop("Contrast file must contain 'groupA' and 'groupB' columns", call. = FALSE)
}

results_summary <- list()

for (i in seq_len(nrow(contrasts))) {
  groupA <- contrasts$groupA[i]
  groupB <- contrasts$groupB[i]
  label <- sprintf("%s_vs_%s", groupA, groupB)
  de_file <- file.path(opt$`de_dir`, sprintf("DE_%s.tsv", label))
  if (!file.exists(de_file)) {
    warning(sprintf("DE file not found for contrast %s", label))
    next
  }
  de_res <- readr::read_tsv(de_file, show_col_types = FALSE)
  if (!score_column %in% names(de_res)) {
    stop(sprintf("Score column '%s' not found in DE table for %s", score_column, label), call. = FALSE)
  }
  if (!"gene_id" %in% names(de_res)) {
    stop("DE table must contain 'gene_id' column", call. = FALSE)
  }
  stats <- de_res %>% dplyr::filter(!is.na(.data[[score_column]])) %>% dplyr::arrange(desc(.data[[score_column]]))
  ranks <- stats %>% dplyr::pull(score_column)
  names(ranks) <- stats$gene_id
  ranks <- stats::na.omit(ranks)
  if (length(ranks) == 0) {
    warning(sprintf("No ranks available for %s", label))
    next
  }
  fgsea_res <- fgsea(pathways = pathways, stats = ranks, minSize = min_size, maxSize = max_size, nperm = nperm)
  fgsea_res <- fgsea_res %>% arrange(padj) %>% mutate(comparison = label)
  out_file <- file.path(opt$outdir, sprintf("fgsea_%s.tsv", label))
  readr::write_tsv(fgsea_res, out_file)
  results_summary[[label]] <- fgsea_res

  top_pathway <- fgsea_res %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(1)
  if (nrow(top_pathway) > 0 && is.finite(top_pathway$padj[1]) && top_pathway$padj[1] <= padj_cutoff) {
    pathway_name <- top_pathway$pathway[1]
    plot <- plotEnrichment(pathways[[pathway_name]], ranks) + ggtitle(sprintf("%s (%s)", pathway_name, label))
    ggsave(filename = file.path(opt$figdir, sprintf("FGSEA_%s_%s.png", label, gsub("[/ ]", "_", pathway_name))), plot = plot, width = 6, height = 4, dpi = 150)
  }
}

if (length(results_summary) > 0) {
  combined <- dplyr::bind_rows(results_summary)
  readr::write_tsv(combined, file.path(opt$outdir, "fgsea_summary.tsv"))
} else {
  readr::write_tsv(tibble::tibble(), file.path(opt$outdir, "fgsea_summary.tsv"))
}
