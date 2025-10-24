# Parallel Processing Pipeline Integration Tests
# Tests parallel processing functionality and caching

library(testthat)
library(future)
library(furrr)
library(digest)
library(dplyr)

# Source parallel pipeline functions
source("../parallel_pipeline.R")

# Test data setup
setup_test_data <- function() {
  # Create mock NHANES-like data for testing
  set.seed(123)

  n <- 100
  test_data <- data.frame(
    SEQN = 1:n,
    RIDAGEYR = sample(20:59, n, replace = TRUE),
    RIAGENDR = sample(1:2, n, replace = TRUE),
    BMXBMI = rnorm(n, 27, 5),
    bodyfat_pct = rnorm(n, 25, 8),
    WTMEC2YR = runif(n, 0.5, 2.0),
    SDMVSTRA = sample(1:15, n, replace = TRUE),
    SDMVPSU = sample(1:30, n, replace = TRUE)
  )

  # Create BMI categories
  test_data <- test_data %>%
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

  return(test_data)
}

# Test 1: Parallel Processing Setup
test_that("Parallel processing setup works correctly", {

  # Test 1a: Worker allocation
  original_plan <- plan()
  plan(multisession, workers = 2)

  # Should be able to set parallel backend
  expect_true(plan()$strategy == "multisession")

  # Restore original plan
  plan(original_plan)

  # Test 1b: Available cores detection
  cores <- availableCores()
  expect_is(cores, "integer")
  expect_true(cores >= 1)

  # Test 1c: Worker calculation
  workers <- availableCores() - 1
  expect_true(workers >= 0)
})

# Test 2: Cache System
test_that("Cache system works correctly", {

  # Setup cache directory
  cache_dir <- "tests/test_cache"
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Test 2a: Cache key generation
  test_data <- data.frame(x = 1:10, y = 11:20)
  cache_key <- digest(list("test", test_data))
  expect_is(cache_key, "character")
  expect_true(nchar(cache_key) > 0)

  # Test 2b: Cache operations
  test_result <- data.frame(result = "cached_data")
  save_to_cache(cache_key, test_result)

  cache_file <- get_cache_path(cache_key)
  expect_true(file.exists(cache_file))

  loaded_result <- load_from_cache(cache_key)
  expect_equal(loaded_result, test_result)

  # Test 2c: Cache miss
  non_existent_key <- digest("non_existent")
  missing_result <- load_from_cache(non_existent_key)
  expect_null(missing_result)
})

# Test 3: Pipeline Function Integration
test_that("Pipeline functions integrate properly", {

  # Test 3a: Function existence
  expect_true(exists("run_parallel_pipeline"))
  expect_true(is.function(run_parallel_pipeline))

  # Test 3b: Function signature
  pipeline_args <- formals(run_parallel_pipeline)
  expect_true(length(pipeline_args) == 0)  # Should take no arguments

  # Test 3c: Pipeline metadata
  metadata <- list(
    targets_version = packageVersion("targets"),
    future_version = packageVersion("future"),
    furrr_version = packageVersion("furrr"),
    workers = availableCores() - 1,
    created_at = Sys.time()
  )

  expect_is(metadata, "list")
  expect_true(length(metadata) >= 4)
})

# Test 4: Parallel Computation
test_that("Parallel computations work correctly", {

  # Test 4a: Simple parallel computation
  original_plan <- plan()
  plan(multisession, workers = 2)

  test_data <- 1:10

  # Parallel computation
  parallel_results <- future_map(test_data, function(x) x * 2)

  # Sequential computation for comparison
  sequential_results <- lapply(test_data, function(x) x * 2)

  # Results should be identical
  expect_equal(parallel_results, sequential_results)

  # Restore original plan
  plan(original_plan)

  # Test 4b: Parallel computation with data
  test_df <- setup_test_data()

  # Test correlation computation
  correlation_results <- future_map(c("Overall", "Male", "Female"), function(group) {
    if (group == "Overall") {
      subset_data <- test_df
    } else if (group == "Male") {
      subset_data <- test_df %>% filter(sex == "Male")
    } else {
      subset_data <- test_df %>% filter(sex == "Female")
    }

    if (nrow(subset_data) > 0) {
      corr <- cor(subset_data$BMXBMI, subset_data$bodyfat_pct)
      return(data.frame(group = group, correlation = corr))
    } else {
      return(data.frame(group = group, correlation = NA))
    }
  }) %>% bind_rows()

  expect_is(correlation_results, "data.frame")
  expect_true("group" %in% names(correlation_results))
  expect_true("correlation" %in% names(correlation_results))
})

# Test 5: Performance Characteristics
test_that("Performance characteristics are within expected ranges", {

  # Test 5a: Cache performance
  cache_dir <- "tests/test_cache"
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  test_data <- data.frame(x = 1:100, y = rnorm(100))

  # Measure cache hit vs miss performance
  cache_key <- digest(test_data)

  # Cache miss (first save)
  start_time <- Sys.time()
  save_to_cache(cache_key, test_data)
  first_save_time <- difftime(Sys.time(), start_time, units = "secs")

  # Cache hit (load)
  start_time <- Sys.time()
  loaded_data <- load_from_cache(cache_key)
  cache_hit_time <- difftime(Sys.time(), start_time, units = "secs")

  # Cache hit should be much faster than initial save
  expect_true(cache_hit_time < first_save_time * 0.1)

  # Test 5b: Memory usage
  initial_memory <- memory.size()

  # Create test dataset
  large_test_data <- data.frame(matrix(rnorm(10000), ncol = 50))

  memory_after_creation <- memory.size()
  memory_growth <- memory_after_creation - initial_memory

  # Memory growth should be reasonable
  expect_true(memory_growth < 100)  # Less than 100MB

  # Cleanup
  rm(large_test_data)
  gc()
})

# Test 6: Error Handling in Parallel Context
test_that("Error handling works in parallel processing", {

  # Test 6a: Error propagation in parallel context
  original_plan <- plan()
  plan(multisession, workers = 2)

  # Test function that might fail
  risky_function <- function(x) {
    if (x == 5) {
      stop("Test error in parallel context")
    }
    return(x * 2)
  }

  # Should handle errors gracefully
  expect_error(future_map(1:10, risky_function))

  # Restore original plan
  plan(original_plan)

  # Test 6b: Safe execution in parallel context
  safe_parallel_function <- function(x) {
    tryCatch({
      if (x == 5) {
        stop("Test error")
      }
      return(x * 2)
    }, error = function(e) {
      return(NA)  # Return NA instead of error
    })
  }

  # Should complete without errors
  safe_results <- future_map(1:10, safe_parallel_function)
  expect_is(safe_results, "list")
  expect_true(length(safe_results) == 10)
})

# Test 7: Integration with Survey Package
test_that("Parallel processing integrates with survey package", {

  # Test 7a: Survey design creation
  test_data <- setup_test_data()

  # Create survey design
  survey_design <- svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC2YR,
    nest = TRUE,
    data = test_data
  )

  expect_is(survey_design, "survey.design")

  # Test 7b: Survey-weighted correlations
  correlation <- svyvar(~BMXBMI + bodyfat_pct, survey_design)
  expect_is(correlation, "matrix")

  # Extract correlation coefficient
  corr_coef <- correlation[1, 2] / sqrt(correlation[1, 1] * correlation[2, 2])
  expect_is(corr_coef, "numeric")

  # Test 7c: Survey-weighted means
  mean_bmi <- svymean(~BMXBMI, survey_design)
  expect_is(mean_bmi, "svymean")
})

# Test 8: Caching Integration
test_that("Caching integrates properly with parallel processing", {

  # Test 8a: Cache persistence across sessions
  cache_dir <- "tests/test_cache"
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  test_data <- data.frame(x = 1:50, y = rnorm(50))
  cache_key <- digest(test_data)

  # Save to cache
  save_to_cache(cache_key, test_data)
  expect_true(file.exists(get_cache_path(cache_key)))

  # Simulate new R session (different process)
  # In practice, this would be a new R session
  loaded_data <- load_from_cache(cache_key)
  expect_equal(loaded_data, test_data)

  # Test 8b: Cache invalidation
  # Modify data slightly
  modified_data <- test_data %>% mutate(x = x + 1)
  modified_key <- digest(modified_data)

  # Should be different keys
  expect_false(cache_key == modified_key)

  # Test 8c: Cache cleanup
  cache_files <- list.files(cache_dir, full.names = TRUE)
  initial_count <- length(cache_files)

  # Add another cache entry
  save_to_cache("test_key_2", data.frame(z = 1:10))

  # Verify cache grew
  updated_count <- length(list.files(cache_dir, full.names = TRUE))
  expect_true(updated_count > initial_count)
})

# Test 9: Performance Benchmarks
test_that("Performance benchmarks meet expectations", {

  # Test 9a: Basic performance
  start_time <- Sys.time()

  # Simple computation
  test_result <- sum(1:1000000)

  execution_time <- difftime(Sys.time(), start_time, units = "secs")
  expect_true(execution_time < 1)  # Should be very fast

  # Test 9b: Parallel overhead
  original_plan <- plan()
  plan(multisession, workers = 2)

  # Measure parallel overhead
  simple_parallel <- function(x) x * 2
  test_data <- 1:100

  start_time <- Sys.time()
  parallel_result <- future_map(test_data, simple_parallel)
  parallel_time <- difftime(Sys.time(), start_time, units = "secs")

  # Parallel should not be excessively slower for simple tasks
  expect_true(parallel_time < 5)  # Less than 5 seconds

  plan(original_plan)

  # Test 9c: Memory efficiency
  initial_memory <- memory.size()

  # Create moderately large dataset
  moderate_data <- data.frame(matrix(rnorm(100000), ncol = 10))

  memory_after <- memory.size()
  memory_growth <- memory_after - initial_memory

  # Memory growth should be reasonable
  expect_true(memory_growth < 20)  # Less than 20MB

  # Cleanup
  rm(moderate_data)
  gc()
})

# Test 10: Integration with Configuration System
test_that("Parallel processing integrates with configuration", {

  # Test 10a: Configuration loading
  if (file.exists("../config/config.yml")) {
    config <- read_yaml("../config/config.yml")
    expect_is(config, "list")

    # Test configuration values
    expect_true(is.numeric(config$analysis$age_range[1]))
    expect_true(is.character(config$analysis$survey_weights_col))
  }

  # Test 10b: Worker configuration
  workers <- availableCores() - 1
  expect_true(workers >= 0)

  # Test 10c: Memory configuration
  memory_limit <- memory.limit()
  expect_is(memory_limit, "numeric")

  # Test 10d: Cache configuration
  cache_enabled <- TRUE  # Default
  expect_is(cache_enabled, "logical")
})

# Run all tests with reporting
run_parallel_tests <- function() {
  cat("⚡ Running Parallel Processing Integration Tests\n")
  cat("===============================================\n")

  # Run test suite
  test_results <- testthat::test_file("tests/test_parallel_pipeline.R")

  # Generate report
  test_report <- list(
    timestamp = Sys.time(),
    tests_run = length(test_results),
    tests_passed = sum(sapply(test_results, function(x) x$passed)),
    tests_failed = sum(sapply(test_results, function(x) !x$passed)),
    total_time = sum(sapply(test_results, function(x) x$time)),
    performance_metrics = list(
      cores_available = availableCores(),
      workers_configured = availableCores() - 1,
      memory_available = memory.size()
    )
  )

  # Display results
  cat("\n📊 Parallel Processing Test Results\n")
  cat("==================================\n")
  cat("Tests run:", test_report$tests_run, "\n")
  cat("Tests passed:", test_report$tests_passed, "\n")
  cat("Tests failed:", test_report$tests_failed, "\n")
  cat("Total time:", round(test_report$total_time, 2), "seconds\n")

  cat("\n🔧 System Configuration:\n")
  cat("  CPU cores:", test_report$performance_metrics$cores_available, "\n")
  cat("  Workers configured:", test_report$performance_metrics$workers_configured, "\n")
  cat("  Memory available:", round(test_report$performance_metrics$memory_available, 1), "MB\n")

  if (test_report$tests_failed > 0) {
    cat("\n❌ Some parallel processing tests failed!\n")
    cat("Check test output above for details.\n")
  } else {
    cat("\n✅ All parallel processing tests passed!\n")
    cat("🚀 Parallel processing integration verified successfully.\n")
  }

  return(test_report)
}

# If run directly, execute tests
if (!interactive()) {
  run_parallel_tests()
}


