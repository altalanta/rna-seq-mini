#!/usr/bin/env Rscript
# NHANES BMI Body Fat Analysis Pipeline Runner
# This script helps run the targets-based pipeline

# Load required packages
library(targets)
library(future)

# Set up parallel processing
plan(multisession, workers = availableCores() - 1)

# Load the pipeline
source("_targets.R")

# Check if pipeline is up to date
outdated_targets <- tar_outdated()

if (length(outdated_targets) == 0) {
  cat("Pipeline is up to date! All targets are current.\n")
} else {
  cat("Running outdated targets:\n")
  print(outdated_targets)

  # Run the pipeline
  cat("\nRunning pipeline...\n")
  start_time <- Sys.time()
  tar_make()
  end_time <- Sys.time()

  cat(paste("\nPipeline completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n"))

  # Show results
  cat("\nPipeline status:\n")
  tar_progress()

  cat("\nLatest results available in outputs/ directory\n")
}

# Provide helpful information
cat("\nUseful commands:\n")
cat("  make pipeline        - Run full pipeline\n")
cat("  make pipeline-vis    - Visualize pipeline graph\n")
cat("  make pipeline-status - Show pipeline progress\n")
cat("  make clean           - Clean output files\n")
cat("  make help            - Show all available targets\n")



