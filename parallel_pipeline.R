#!/usr/bin/env Rscript
# NHANES BMI Body Fat Analysis - Parallel Processing Pipeline
# Alternative implementation using future/furrr for parallelization

# Load required packages
library(future)
library(furrr)
library(dplyr)
library(ggplot2)
library(survey)
library(foreign)
library(readr)
library(yaml)
library(digest) # For caching

# Source data versioning utilities
source("R/data_versioning.R")

# Set up parallel processing
plan(multisession, workers = availableCores() - 1)

# Load configuration
config <- read_yaml("config/config.yml")

# Source utility functions
source("R/load_config.R")
source("R/error_handling.R")
source("R/data_validation.R")

# Cache directory
cache_dir <- "cache"
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

# Cache utility functions
get_cache_path <- function(name) {
  file.path(cache_dir, paste0(digest(name), ".rds"))
}

load_from_cache <- function(name) {
  cache_file <- get_cache_path(name)
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  return(NULL)
}

save_to_cache <- function(name, data) {
  cache_file <- get_cache_path(name)
  saveRDS(data, cache_file)
  return(data)
}

# Main pipeline function
run_parallel_pipeline <- function() {
  start_time <- Sys.time()
  cat("Starting parallel NHANES BMI Body Fat analysis pipeline...\n")

  # Step 0: Data version management and integrity checks
  cat("Step 0: Data version management and integrity checks...\n")

  # Initialize data registry if needed
  if (!file.exists("data/registry/data_registry.json")) {
    cat("  Initializing data registry...\n")
    initialize_data_registry()
  }

  # Check for data updates
  cat("  Checking for data updates...\n")
  update_check <- check_for_updates()

  # Validate data integrity
  cat("  Validating data integrity...\n")
  integrity_check <- validate_data_integrity()

  if (!integrity_check$valid) {
    warning("⚠️ Data integrity issues detected. Analysis may produce unreliable results.")
    cat("💡 Run 'make data-health' to see detailed quality report.\n")
  }

  # Generate quality report
  quality_report <- generate_quality_report()

  cat("  Data quality report generated.\n")

  # Step 1: Data fetching and loading (with caching)
  cat("Step 1: Loading NHANES datasets...\n")
  datasets_key <- "nhanes_datasets"
  datasets <- load_from_cache(datasets_key)

  if (is.null(datasets)) {
    # Load datasets
    demo_data <- safe_read_xpt(file.path(config$data$raw_dir, config$nhanes$demo_file), "DEMO_J dataset")
    bmx_data <- safe_read_xpt(file.path(config$data$raw_dir, config$nhanes$bmx_file), "BMX_J dataset")
    dxx_data <- safe_read_xpt(file.path(config$data$raw_dir, config$nhanes$dxx_file), "DXX_J dataset")
    dxxag_data <- safe_read_xpt(file.path(config$data$raw_dir, config$nhanes$dxxag_file), "DXXAG_J dataset")

    datasets <- list(
      demo = demo_data,
      bmx = bmx_data,
      dxx = dxx_data,
      dxxag = dxxag_data
    )

    save_to_cache(datasets_key, datasets)
    cat("  Datasets loaded and cached.\n")
  } else {
    cat("  Datasets loaded from cache.\n")
  }

  # Step 2: Identify body fat variables (parallel)
  cat("Step 2: Identifying body fat variables...\n")
  variables_key <- "bodyfat_variables"
  variables <- load_from_cache(variables_key)

  if (is.null(variables)) {
    # Identify variables in parallel
    variables <- future_map(list(datasets$dxx, datasets$dxxag), function(dataset) {
      if (is.null(dataset) || nrow(dataset) == 0) return(NULL)

      # Identify body fat variable
      dxx_labels <- attr(dataset, "var.labels")
      if (is.null(dxx_labels)) {
        dxx_labels <- sapply(names(dataset), function(x) {
          label <- attr(dataset[[x]], "label")
          if (is.null(label)) return("")
          return(label)
        })
      }

      fat_patterns <- c("percent fat", "%fat", "total % fat", "% fat")
      bodyfat_vars <- c()

      for (pattern in fat_patterns) {
        matches <- grep(pattern, dxx_labels, ignore.case = TRUE)
        if (length(matches) > 0) {
          bodyfat_vars <- c(bodyfat_vars, names(dataset)[matches])
        }
      }

      common_names <- c("DXDTOFAT", "DXDTOPF", "DXDPFAT")
      for (name in common_names) {
        if (name %in% names(dataset)) {
          bodyfat_vars <- c(bodyfat_vars, name)
        }
      }

      if (length(bodyfat_vars) > 0) {
        return(list(bodyfat = bodyfat_vars[1]))
      }

      # Look for android/gynoid fat
      android_matches <- grep("android.*fat", dxx_labels, ignore.case = TRUE)
      gynoid_matches <- grep("gynoid.*fat", dxx_labels, ignore.case = TRUE)

      result <- list()
      if (length(android_matches) > 0) {
        result$android <- names(dataset)[android_matches[1]]
      }
      if (length(gynoid_matches) > 0) {
        result$gynoid <- names(dataset)[gynoid_matches[1]]
      }

      return(result)
    })

    # Extract variables
    bodyfat_var <- variables[[1]]$bodyfat
    android_var <- variables[[2]]$android
    gynoid_var <- variables[[2]]$gynoid

    variables <- list(
      bodyfat = bodyfat_var,
      android = android_var,
      gynoid = gynoid_var
    )

    save_to_cache(variables_key, variables)
    cat("  Variables identified and cached.\n")
  } else {
    cat("  Variables loaded from cache.\n")
  }

  # Step 3: Merge and clean data (with caching)
  cat("Step 3: Merging and cleaning datasets...\n")
  cleaned_key <- "cleaned_data"
  cleaned_data <- load_from_cache(cleaned_key)

  if (is.null(cleaned_data)) {
    # Merge datasets
    nhanes <- datasets$demo %>%
      select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, WTMEC2YR, SDMVSTRA, SDMVPSU) %>%
      left_join(datasets$bmx %>% select(SEQN, BMXBMI), by = "SEQN") %>%
      left_join(datasets$dxx %>% select_at(c("SEQN", variables$bodyfat)), by = "SEQN")

    if (!is.null(variables$android) && !is.null(variables$gynoid)) {
      nhanes <- nhanes %>%
        left_join(datasets$dxxag %>% select_at(c("SEQN", variables$android, variables$gynoid)), by = "SEQN")
    }

    # Apply inclusion criteria
    nhanes_adults <- nhanes %>%
      filter(RIDAGEYR >= config$analysis$age_range[1] &
             RIDAGEYR <= config$analysis$age_range[2])

    nhanes_complete <- nhanes_adults %>%
      filter(!is.na(BMXBMI) & !is.na(bodyfat_pct) &
             !is.na(WTMEC2YR) & !is.na(SDMVSTRA) & !is.na(SDMVPSU))

    # Create BMI categories and sex factor
    nhanes_complete <- nhanes_complete %>%
      mutate(
        bmi_cat = case_when(
          BMXBMI < 18.5 ~ "Underweight",
          BMXBMI >= 18.5 & BMXBMI < 25 ~ "Normal",
          BMXBMI >= 25 & BMXBMI < 30 ~ "Overweight",
          BMXBMI >= 30 & BMXBMI < 35 ~ "Obesity I",
          BMXBMI >= 35 & BMXBMI < 40 ~ "Obesity II",
          BMXBMI >= 40 ~ "Obesity III"
        ),
        bmi_cat = factor(bmi_cat, levels = c("Underweight", "Normal", "Overweight",
                                             "Obesity I", "Obesity II", "Obesity III")),
        sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female"))
      )

    # Rename body fat variable
    names(nhanes_complete)[names(nhanes_complete) == variables$bodyfat] <- "bodyfat_pct"

    cleaned_data <- nhanes_complete
    save_to_cache(cleaned_key, cleaned_data)
    cat("  Data cleaned and cached.\n")
  } else {
    cat("  Cleaned data loaded from cache.\n")
  }

  # Step 4: Survey design (with caching)
  cat("Step 4: Creating survey design...\n")
  design_key <- "survey_design"
  survey_design <- load_from_cache(design_key)

  if (is.null(survey_design)) {
    validate_survey_design(
      cleaned_data,
      config$analysis$survey_weights_col,
      config$analysis$strata_col,
      config$analysis$psu_col
    )

    survey_design <- svydesign(
      ids = as.formula(paste("~", config$analysis$psu_col)),
      strata = as.formula(paste("~", config$analysis$strata_col)),
      weights = as.formula(paste("~", config$analysis$survey_weights_col)),
      nest = TRUE,
      data = cleaned_data
    )

    save_to_cache(design_key, survey_design)
    cat("  Survey design created and cached.\n")
  } else {
    cat("  Survey design loaded from cache.\n")
  }

  # Step 5: Parallel correlation analysis
  cat("Step 5: Computing correlations in parallel...\n")
  corr_key <- "correlation_results"
  correlation_results <- load_from_cache(corr_key)

  if (is.null(correlation_results)) {
    correlation_results <- future_map(c("Overall", "Male", "Female"), function(group) {
      if (group == "Overall") {
        design_subset <- survey_design
      } else if (group == "Male") {
        design_subset <- subset(survey_design, sex == "Male")
      } else {
        design_subset <- subset(survey_design, sex == "Female")
      }

      corr_data <- svyvar(~BMXBMI + bodyfat_pct, design_subset)
      correlation <- corr_data[1,2] / sqrt(corr_data[1,1] * corr_data[2,2])

      # Standard error using delta method
      corr_se <- sqrt((1 - correlation^2)^2 / (4 * correlation^2) *
                      (corr_data[1,1]/corr_data[1,2]^2 + corr_data[2,2]/corr_data[1,2]^2 -
                       2/(corr_data[1,1] * corr_data[2,2])))

      data.frame(
        group = group,
        correlation = correlation,
        std_error = corr_se,
        ci_lower = correlation - 1.96 * corr_se,
        ci_upper = correlation + 1.96 * corr_se
      )
    }) %>% bind_rows()

    save_to_cache(corr_key, correlation_results)
    cat("  Correlations computed and cached.\n")
  } else {
    cat("  Correlations loaded from cache.\n")
  }

  # Step 6: Parallel BMI class analysis
  cat("Step 6: Computing BMI class statistics in parallel...\n")
  bmi_key <- "bmi_class_results"
  bmi_class_results <- load_from_cache(bmi_key)

  if (is.null(bmi_class_results)) {
    # Create BMI-sex combinations for parallel processing
    bmi_sex_grid <- expand.grid(
      bmi_cat = levels(cleaned_data$bmi_cat),
      sex = levels(cleaned_data$sex),
      stringsAsFactors = FALSE
    )

    bmi_class_results <- future_pmap(bmi_sex_grid, function(bmi_cat, sex) {
      subset_design <- subset(survey_design, bmi_cat == bmi_cat & sex == sex)

      if (nrow(subset_design$variables) > 0) {
        n_unweighted <- nrow(subset_design$variables)
        pop_total <- sum(weights(subset_design))

        mean_bf <- svymean(~bodyfat_pct, subset_design)
        mean_est <- as.numeric(mean_bf)
        mean_se <- as.numeric(SE(mean_bf))
        mean_ci <- confint(mean_bf)

        q05 <- as.numeric(svyquantile(~bodyfat_pct, subset_design, quantiles = 0.05)[[1]])
        q50 <- as.numeric(svyquantile(~bodyfat_pct, subset_design, quantiles = 0.50)[[1]])
        q95 <- as.numeric(svyquantile(~bodyfat_pct, subset_design, quantiles = 0.95)[[1]])

        data.frame(
          bmi_cat = bmi_cat,
          sex = sex,
          n_unweighted = n_unweighted,
          pop_total = pop_total,
          mean_bodyfat = mean_est,
          mean_se = mean_se,
          mean_ci_lower = as.numeric(mean_ci[1]),
          mean_ci_upper = as.numeric(mean_ci[2]),
          q05 = q05,
          q50 = q50,
          q95 = q95
        )
      } else {
        data.frame(
          bmi_cat = bmi_cat,
          sex = sex,
          n_unweighted = 0,
          pop_total = 0,
          mean_bodyfat = NA,
          mean_se = NA,
          mean_ci_lower = NA,
          mean_ci_upper = NA,
          q05 = NA,
          q50 = NA,
          q95 = NA
        )
      }
    }) %>% bind_rows()

    save_to_cache(bmi_key, bmi_class_results)
    cat("  BMI class statistics computed and cached.\n")
  } else {
    cat("  BMI class statistics loaded from cache.\n")
  }

  # Step 7: Linearity assessment
  cat("Step 7: Assessing linearity...\n")
  linear_key <- "linearity_results"
  linearity_results <- load_from_cache(linear_key)

  if (is.null(linearity_results)) {
    linear_model <- svyglm(bodyfat_pct ~ BMXBMI, design = survey_design)
    quad_model <- svyglm(bodyfat_pct ~ BMXBMI + I(BMXBMI^2), design = survey_design)

    quad_p <- summary(quad_model)$coefficients["I(BMXBMI^2)", "Pr(>|t|)"]

    linearity_results <- list(
      linear_aic = AIC(linear_model),
      quad_aic = AIC(quad_model),
      quad_p_value = quad_p,
      significant_nonlinearity = quad_p < 0.05
    )

    save_to_cache(linear_key, linearity_results)
    cat("  Linearity assessment completed and cached.\n")
  } else {
    cat("  Linearity results loaded from cache.\n")
  }

  # Step 8: Generate visualizations (with caching)
  cat("Step 8: Creating visualizations...\n")
  plots_key <- "visualization_files"
  plot_files <- load_from_cache(plots_key)

  if (is.null(plot_files)) {
    # Main scatter plot
    plot_data <- cleaned_data %>%
      mutate(plot_weight = WTMEC2YR / sum(WTMEC2YR) * nrow(cleaned_data))

    main_plot <- ggplot(plot_data, aes(x = BMXBMI, y = bodyfat_pct)) +
      geom_point(aes(size = plot_weight), alpha = 0.3) +
      geom_smooth(aes(weight = plot_weight), method = "loess", se = TRUE) +
      facet_wrap(~sex) +
      labs(
        title = "BMI vs Whole-Body % Body Fat by Sex",
        subtitle = "U.S. civilian non-institutionalized adults (20-59), NHANES 2017-2018",
        x = "Body Mass Index (kg/m²)",
        y = "Whole-Body % Body Fat (DXA)",
        size = "Survey Weight"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")

    main_plot_file <- file.path(config$outputs$figures_dir, "bmi_vs_bodyfat_plot_sex_facets.png")
    ggsave(main_plot_file, main_plot, width = 10, height = 6, dpi = 300)

    # BMI class plot
    bmi_plot <- ggplot(bmi_class_results, aes(x = bmi_cat, y = mean_bodyfat, fill = sex)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = mean_ci_lower, ymax = mean_ci_upper),
                    position = position_dodge(width = 0.9), width = 0.25) +
      labs(
        title = "Mean Body Fat by BMI Class and Sex",
        x = "BMI Category",
        y = "Mean % Body Fat",
        fill = "Sex"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    bmi_plot_file <- file.path(config$outputs$figures_dir, "bodyfat_by_bmi_class.png")
    ggsave(bmi_plot_file, bmi_plot, width = 10, height = 6, dpi = 300)

    plot_files <- c(main_plot_file, bmi_plot_file)
    save_to_cache(plots_key, plot_files)
    cat("  Visualizations created and cached.\n")
  } else {
    cat("  Visualizations loaded from cache.\n")
  }

  # Step 9: Export results
  cat("Step 9: Exporting results...\n")
  write.csv(correlation_results, file.path(config$outputs$tables_dir, "corr_bmi_bodyfat_overall_and_by_sex.csv"), row.names = FALSE)
  write.csv(bmi_class_results, file.path(config$outputs$tables_dir, "bodyfat_by_bmi_class_by_sex.csv"), row.names = FALSE)

  # Create methods documentation
  methods_text <- paste0(
    "NHANES 2017-2018 BMI vs % Body Fat Analysis (Parallel Pipeline)\n",
    "================================================================\n\n",
    "Pipeline: Custom parallel processing with caching and data versioning\n",
    "Workers: ", availableCores() - 1, "\n",
    "Analysis date: ", Sys.Date(), "\n",
    "Processing time: ", round(difftime(Sys.time(), start_time, units = "secs"), 2), " seconds\n\n",
    "Data Management:\n",
    "  Registry: ", DATA_REGISTRY_FILE, "\n",
    "  Manifest: ", DATA_MANIFEST_FILE, "\n",
    "  Integrity: ", ifelse(quality_report$integrity_check$valid, "VALID", "ISSUES DETECTED"), "\n",
    "  Total files: ", quality_report$registry_summary$total_files, "\n"
  )

  writeLines(methods_text, file.path(config$outputs$logs_dir, "parallel_pipeline_methods.txt"))
  cat("  Results exported.\n")

  # Final summary
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "secs")

  cat("\n" %>% paste(strrep("=", 60), "\n"))
  cat("Pipeline completed successfully!\n")
  cat(paste("Total processing time:", round(total_time, 2), "seconds\n"))
  cat(paste("Workers used:", availableCores() - 1, "\n"))
  cat("Results available in outputs/ directory\n")
  cat(strrep("=", 60) %>% paste("\n"))

  return(list(
    correlation_results = correlation_results,
    bmi_class_results = bmi_class_results,
    linearity_results = linearity_results,
    plot_files = plot_files,
    processing_time = total_time
  ))
}

# Run pipeline if script is executed directly
if (!interactive()) {
  results <- run_parallel_pipeline()
}
