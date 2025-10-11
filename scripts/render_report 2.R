#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(rmarkdown)
})

option_list <- list(
  make_option(c("--params"), type = "character", default = "config/params.yaml"),
  make_option(c("--template"), type = "character", default = "report/rnaseq_report.Rmd"),
  make_option(c("--output"), type = "character", help = "Path to output HTML"),
  make_option(c("--results-dir"), type = "character", default = "results")
)

opt <- parse_args(OptionParser(option_list = option_list))
params <- yaml::read_yaml(opt$params)
if (is.null(opt$output) || is.na(opt$output)) {
  stop("--output must be provided", call. = FALSE)
}

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

render(
  input = opt$template,
  output_file = opt$output,
  params = list(
    params_yaml = params,
    results_dir = opt$`results_dir`
  ),
  envir = new.env(parent = globalenv())
)
