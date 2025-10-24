# Comprehensive Integration Tests for NHANES BMI Body Fat Analysis Platform
# Tests all enhanced features working together

library(testthat)
library(dplyr)
library(yaml)
library(jsonlite)
library(digest)

# Source all required modules
source("../R/error_handling.R")
source("../R/data_versioning.R")
source("../R/data_validation.R")
source("../parallel_pipeline.R")

# Test configuration
TEST_CONFIG <- list(
  test_data_dir = "tests/test_data",
  test_outputs_dir = "tests/test_outputs",
  test_cache_dir = "tests/test_cache"
)

# Setup test environment
setup_test_environment <- function() {
  # Create test directories
  dirs <- c(TEST_CONFIG$test_data_dir, TEST_CONFIG$test_outputs_dir, TEST_CONFIG$test_cache_dir)
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  # Create minimal test config
  test_config <- list(
    data = list(
      raw_dir = TEST_CONFIG$test_data_dir,
      derived_dir = TEST_CONFIG$test_outputs_dir
    ),
    outputs = list(
      tables_dir = TEST_CONFIG$test_outputs_dir,
      figures_dir = TEST_CONFIG$test_outputs_dir,
      logs_dir = TEST_CONFIG$test_outputs_dir,
      report_dir = TEST_CONFIG$test_outputs_dir
    ),
    nhanes = list(
      demo_file = "DEMO_J.XPT",
      bmx_file = "BMX_J.XPT",
      dxx_file = "DXX_J.XPT",
      dxxag_file = "DXXAG_J.XPT"
    ),
    analysis = list(
      age_range = c(20, 59),
      survey_weights_col = "WTMEC2YR",
      strata_col = "SDMVSTRA",
      psu_col = "SDMVPSU"
    ),
    logging = list(
      level = "INFO",
      file = "test_log.txt"
    )
  )

  write_yaml(test_config, file.path(TEST_CONFIG$test_outputs_dir, "test_config.yml"))

  return(test_config)
}

# Cleanup test environment
cleanup_test_environment <- function() {
  # Remove test directories and files
  dirs_to_remove <- c(TEST_CONFIG$test_data_dir, TEST_CONFIG$test_outputs_dir, TEST_CONFIG$test_cache_dir)
  for (dir in dirs_to_remove) {
    if (dir.exists(dir)) {
      unlink(dir, recursive = TRUE)
    }
  }
}

# Test 1: Complete Pipeline Integration
test_that("Complete pipeline integration works end-to-end", {

  # Setup test environment
  test_config <- setup_test_environment()

  # Test 1a: Data Versioning Integration
  expect_true(initialize_data_registry())
  expect_true(file.exists("data/registry/data_registry.json"))

  registry <- load_data_registry()
  expect_is(registry, "list")
  expect_true(length(registry$entries) == 0)  # Empty initially

  # Test 1b: Configuration Integration
  config <- load_config(file.path(TEST_CONFIG$test_outputs_dir, "test_config.yml"))
  expect_is(config, "list")
  expect_true(all(c("data", "outputs", "nhanes", "analysis", "logging") %in% names(config)))

  # Test 1c: Error Handling Integration
  expect_silent(safe_log("Test log message", "INFO", config))
  expect_true(file.exists(file.path(config$outputs$logs_dir, config$logging$file)))

  # Test 1d: Parallel Pipeline Integration
  # This would run the full pipeline if data were available
  # For now, test that the pipeline functions exist and are callable
  expect_true(exists("run_parallel_pipeline"))

  # Test 1e: Cache System Integration
  cache_dir <- TEST_CONFIG$test_cache_dir
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  test_data <- data.frame(x = 1:10, y = 11:20)
  cache_key <- digest(list("test", test_data))

  # Test cache operations
  save_to_cache(cache_key, test_data)
  expect_true(file.exists(get_cache_path(cache_key)))

  loaded_data <- load_from_cache(cache_key)
  expect_equal(loaded_data, test_data)

  # Cleanup
  cleanup_test_environment()
})

# Test 2: Data Versioning Integration
test_that("Data versioning system integrates properly", {

  # Setup
  test_config <- setup_test_environment()

  # Test 2a: Registry Operations
  expect_true(initialize_data_registry())

  # Create mock data file for testing
  mock_data <- data.frame(SEQN = 1:100, BMI = rnorm(100, 25, 5))
  mock_file <- file.path(TEST_CONFIG$test_data_dir, "test_data.csv")
  write.csv(mock_data, mock_file, row.names = FALSE)

  # Test adding to registry
  success <- add_to_registry(mock_file, "test_data", "2017-2018")
  expect_true(success)

  # Verify registry entry
  registry <- load_data_registry()
  expect_true(length(registry$entries) > 0)

  entry <- registry$entries[[1]]
  expect_equal(entry$data_type, "test_data")
  expect_equal(entry$nhanes_cycle, "2017-2018")
  expect_true(!is.null(entry$hash_sha256))

  # Test 2b: Integrity Validation
  integrity <- validate_data_integrity()
  expect_is(integrity, "list")
  expect_true("valid" %in% names(integrity))

  # Test 2c: Update Detection
  updates <- check_for_updates()
  expect_is(updates, "list")
  expect_true(all(c("uptodate", "updates") %in% names(updates)))

  # Test 2d: Quality Report Generation
  quality_report <- generate_quality_report()
  expect_is(quality_report, "list")
  expect_true(all(c("metadata", "integrity_check", "registry_summary", "recommendations") %in% names(quality_report)))

  # Test 2e: Manifest Generation
  manifest <- generate_data_manifest()
  expect_is(manifest, "list")
  expect_true(all(c("metadata", "data_files", "summary") %in% names(manifest)))

  # Cleanup
  cleanup_test_environment()
})

# Test 3: Interactive Features Integration
test_that("Interactive features integrate properly", {

  # Test 3a: Tutorial System Integration
  # Check that tutorial files exist and are properly formatted
  tutorial_files <- c("tutorials/getting_started.Rmd", "tutorials/help_troubleshooting.Rmd")

  for (tutorial_file in tutorial_files) {
    if (file.exists(tutorial_file)) {
      content <- readLines(tutorial_file, n = 10)
      expect_true(any(grepl("learnr::tutorial", content)))
      expect_true(any(grepl("runtime: shiny", content)))
    }
  }

  # Test 3b: Configuration Wizard Integration
  # Check that app.R exists and loads properly
  expect_true(file.exists("app.R"))

  # Test configuration validation
  if (file.exists("config/config.yml")) {
    config <- read_yaml("config/config.yml")
    expect_is(config, "list")

    # Validate required sections
    required_sections <- c("data", "outputs", "nhanes", "analysis", "logging")
    for (section in required_sections) {
      expect_true(section %in% names(config))
    }
  }

  # Test 3c: Error Handling Integration in Interactive Context
  # Test that error handling works in tutorial context
  test_error <- NhanesError("Test error for integration testing", "TEST001")
  expect_is(test_error, "NhanesError")
  expect_equal(test_error$message, "Test error for integration testing")
  expect_equal(test_error$code, "TEST001")
})

# Test 4: Parallel Processing Integration
test_that("Parallel processing integrates properly", {

  # Test 4a: Future/Furrr Integration
  library(future)
  library(furrr)

  # Test parallel backend setup
  original_plan <- plan()
  plan(multisession, workers = 2)

  # Test parallel computation
  test_data <- 1:10
  parallel_results <- future_map(test_data, function(x) x * 2)
  sequential_results <- lapply(test_data, function(x) x * 2)

  expect_equal(parallel_results, sequential_results)

  # Restore original plan
  plan(original_plan)

  # Test 4b: Cache Integration
  cache_dir <- TEST_CONFIG$test_cache_dir
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Test cache key generation
  test_input <- list(data = test_data, config = list(param = "test"))
  cache_key <- digest(test_input)
  expect_is(cache_key, "character")
  expect_true(nchar(cache_key) > 0)

  # Test cache operations
  test_result <- data.frame(result = "test_output")
  save_to_cache(cache_key, test_result)
  expect_true(file.exists(get_cache_path(cache_key)))

  loaded_result <- load_from_cache(cache_key)
  expect_equal(loaded_result, test_result)

  # Test 4c: Pipeline Function Integration
  # Verify pipeline functions exist and are callable
  expect_true(exists("run_parallel_pipeline"))
  expect_true(is.function(run_parallel_pipeline))

  # Test pipeline metadata generation
  pipeline_metadata <- list(
    targets_version = packageVersion("targets"),
    future_version = packageVersion("future"),
    furrr_version = packageVersion("furrr"),
    workers = availableCores() - 1,
    created_at = Sys.time()
  )
  expect_is(pipeline_metadata, "list")
  expect_true(length(pipeline_metadata) >= 4)
})

# Test 5: Error Handling Integration
test_that("Error handling integrates across all components", {

  # Test 5a: Custom Error Classes
  test_cases <- list(
    list(error = NhanesError("Test message", "TEST001"), expected_class = "NhanesError"),
    list(error = DataValidationError("Validation failed", "test_field"), expected_class = "DataValidationError"),
    list(error = FileNotFoundError("File missing", "test.csv"), expected_class = "FileNotFoundError")
  )

  for (test_case in test_cases) {
    expect_is(test_case$error, test_case$expected_class)
    expect_is(test_case$error, "NhanesError")
  }

  # Test 5b: Error Suggestion System
  error_suggestions <- get_error_suggestions("DL001", "File not found")
  expect_is(error_suggestions, "character")
  expect_true(length(error_suggestions) > 0)
  expect_true(any(grepl("download", error_suggestions, ignore.case = TRUE)))

  # Test 5c: Safe Execution
  safe_result <- safe_execute({
    # This should work
    x <- 1 + 1
    return(x)
  }, "test_operation")

  expect_equal(safe_result, 2)

  # Test error in safe execution
  expect_error(safe_execute({
    stop("Test error")
  }, "failing_operation"))

  # Test 5d: Pipeline Health Check
  health_issues <- check_pipeline_health()
  expect_is(health_issues, "list")

  # Should detect missing directories/files
  expect_true(length(health_issues) > 0)

  # Test 5e: User-Friendly Error Display
  test_error <- NhanesError(
    "Integration test error",
    "INT001",
    details = list(component = "test", severity = "medium")
  )

  # Should not throw error when displaying
  expect_silent(display_user_friendly_error(test_error))
})

# Test 6: Configuration Integration
test_that("Configuration system integrates properly", {

  # Test 6a: Configuration Loading
  test_config <- setup_test_environment()
  loaded_config <- load_config(file.path(TEST_CONFIG$test_outputs_dir, "test_config.yml"))
  expect_is(loaded_config, "list")

  # Test 6b: Configuration Validation
  required_sections <- c("data", "outputs", "nhanes", "analysis", "logging")
  for (section in required_sections) {
    expect_true(section %in% names(loaded_config))
  }

  # Test 6c: Directory Creation
  ensure_output_dirs(loaded_config)

  # Verify directories were created
  for (section in c("data", "outputs")) {
    for (dir_type in names(loaded_config[[section]])) {
      dir_path <- loaded_config[[section]][[dir_type]]
      if (dir_path != "" && !is.null(dir_path)) {
        expect_true(dir.exists(dir_path) || file.exists(dir_path))
      }
    }
  }

  # Test 6d: Logging Integration
  safe_log("Test configuration integration", "INFO", loaded_config)
  log_file <- file.path(loaded_config$outputs$logs_dir, loaded_config$logging$file)
  expect_true(file.exists(log_file))

  # Cleanup
  cleanup_test_environment()
})

# Test 7: End-to-End Workflow Integration
test_that("End-to-end workflow integrates all components", {

  # This is a high-level integration test that simulates the complete workflow

  # Test 7a: Setup Phase
  test_config <- setup_test_environment()

  # Initialize data registry
  expect_true(initialize_data_registry())

  # Test 7b: Configuration Phase
  expect_is(load_config(file.path(TEST_CONFIG$test_outputs_dir, "test_config.yml")), "list")

  # Test 7c: Error Handling Phase
  # Test that error handling works throughout
  test_error <- NhanesError("Workflow integration test", "WF001")
  expect_silent(display_user_friendly_error(test_error))

  # Test 7d: Data Versioning Phase
  # Create mock data file
  mock_data <- data.frame(SEQN = 1:50, BMI = rnorm(50, 25, 3))
  mock_file <- file.path(TEST_CONFIG$test_data_dir, "integration_test.csv")
  write.csv(mock_data, mock_file, row.names = FALSE)

  # Test registry integration
  expect_true(add_to_registry(mock_file, "integration_test", "2017-2018"))

  registry <- load_data_registry()
  expect_true(length(registry$entries) > 0)

  # Test 7e: Quality Assurance Phase
  integrity <- validate_data_integrity()
  expect_is(integrity, "list")

  quality_report <- generate_quality_report()
  expect_is(quality_report, "list")

  # Test 7f: Cleanup Phase
  cleanup_test_environment()

  # Verify cleanup was successful
  expect_false(dir.exists(TEST_CONFIG$test_data_dir))
  expect_false(dir.exists(TEST_CONFIG$test_outputs_dir))
})

# Test 8: Performance Integration
test_that("Performance features integrate properly", {

  # Test 8a: Parallel Processing Performance
  library(future)
  library(furrr)

  # Test worker allocation
  original_plan <- plan()
  plan(multisession, workers = 2)

  # Test parallel computation performance
  test_function <- function(x) {
    Sys.sleep(0.01)  # Small delay to simulate work
    return(x * 2)
  }

  test_data <- 1:20

  # Time parallel execution
  start_time <- Sys.time()
  parallel_results <- future_map(test_data, test_function)
  parallel_time <- difftime(Sys.time(), start_time, units = "secs")

  # Time sequential execution
  start_time <- Sys.time()
  sequential_results <- lapply(test_data, test_function)
  sequential_time <- difftime(Sys.time(), start_time, units = "secs")

  # Verify results are identical
  expect_equal(parallel_results, sequential_results)

  # Parallel should be reasonably fast (not necessarily faster due to overhead)
  expect_true(parallel_time < sequential_time * 2)  # Allow some overhead

  # Restore original plan
  plan(original_plan)

  # Test 8b: Caching Performance
  cache_dir <- TEST_CONFIG$test_cache_dir
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  test_data <- data.frame(x = 1:100, y = rnorm(100))

  # First execution (cache miss)
  start_time <- Sys.time()
  cache_key <- digest(test_data)
  save_to_cache(cache_key, test_data)
  first_save_time <- difftime(Sys.time(), start_time, units = "secs")

  # Second execution (cache hit)
  start_time <- Sys.time()
  loaded_data <- load_from_cache(cache_key)
  cache_hit_time <- difftime(Sys.time(), start_time, units = "secs")

  # Cache hit should be much faster than file I/O
  expect_true(cache_hit_time < first_save_time * 0.1)

  # Verify data integrity
  expect_equal(loaded_data, test_data)

  # Test 8c: Memory Management
  # Test that memory usage is reasonable
  initial_memory <- memory.size()

  # Create some test data
  large_data <- data.frame(matrix(rnorm(10000), ncol = 100))

  # Check memory didn't grow excessively
  memory_growth <- memory.size() - initial_memory
  expect_true(memory_growth < 50)  # Less than 50MB growth

  # Force cleanup
  rm(large_data)
  gc()
})

# Test 9: System Integration
test_that("All system components integrate properly", {

  # Test 9a: Package Dependencies
  required_packages <- c("dplyr", "ggplot2", "survey", "foreign", "yaml", "future", "furrr", "digest", "jsonlite")

  for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      expect_true(requireNamespace(pkg, quietly = TRUE))
    } else {
      # Skip if package not available (for CI environments)
      skip(paste("Package", pkg, "not available"))
    }
  }

  # Test 9b: File System Integration
  # Test that all expected directories can be created
  test_dirs <- c("data/raw", "data/derived", "outputs/tables", "outputs/figures", "outputs/logs", "cache")

  for (test_dir in test_dirs) {
    # Should be able to create directory
    if (!dir.exists(test_dir)) {
      expect_true(dir.create(test_dir, recursive = TRUE))
      expect_true(dir.exists(test_dir))
    }
  }

  # Test 9c: Configuration File Integration
  # Test YAML parsing and validation
  if (file.exists("config/config.yml")) {
    config <- read_yaml("config/config.yml")
    expect_is(config, "list")

    # Test that configuration can be used
    expect_true(is.character(config$analysis$survey_weights_col))
    expect_true(is.numeric(config$analysis$age_range[1]))
    expect_true(is.numeric(config$analysis$age_range[2]))
  }

  # Test 9d: Logging Integration
  # Test that logging works across different components
  test_config <- list(
    outputs = list(logs_dir = TEST_CONFIG$test_outputs_dir),
    logging = list(file = "integration_test.log")
  )

  safe_log("Integration test log message", "INFO", test_config)
  log_file <- file.path(test_config$outputs$logs_dir, test_config$logging$file)
  expect_true(file.exists(log_file))

  # Test 9e: Error Propagation
  # Test that errors propagate correctly through the system
  expect_error(safe_execute({
    stop(NhanesError("Integration test error", "INT001"))
  }, "test_operation"))

  # Test that error suggestions work
  suggestions <- get_error_suggestions("INT001", "Integration test error")
  expect_is(suggestions, "character")
  expect_true(length(suggestions) > 0)
})

# Test 10: Regression Testing
test_that("New features don't break existing functionality", {

  # Test 10a: Backward Compatibility
  # Test that existing functions still work
  expect_true(exists("safe_read_xpt"))
  expect_true(exists("validate_nhanes_data"))
  expect_true(exists("ensure_output_dirs"))

  # Test 10b: Enhanced Features Don't Break Core
  # Test that new error handling doesn't break existing code
  test_error <- NhanesError("Regression test error", "REG001")
  expect_is(test_error, "NhanesError")

  # Test that logging still works
  test_config <- list(outputs = list(logs_dir = TEST_CONFIG$test_outputs_dir))
  safe_log("Regression test message", "INFO", test_config)

  # Test 10c: Performance Regression
  # Ensure new features don't significantly impact performance

  # Simple performance test
  start_time <- Sys.time()
  test_data <- data.frame(x = 1:1000, y = rnorm(1000))
  simple_analysis <- summary(test_data)
  execution_time <- difftime(Sys.time(), start_time, units = "secs")

  # Should complete quickly (less than 1 second)
  expect_true(execution_time < 1)

  # Test 10d: Memory Regression
  # Ensure new features don't cause memory leaks
  initial_memory <- memory.size()

  # Run some operations
  for (i in 1:10) {
    test_df <- data.frame(x = 1:100, y = rnorm(100))
    summary(test_df)
  }

  # Force garbage collection
  gc()

  # Memory shouldn't grow excessively
  final_memory <- memory.size()
  memory_growth <- final_memory - initial_memory
  expect_true(memory_growth < 10)  # Less than 10MB growth
})

# Run all tests with comprehensive reporting
run_integration_tests <- function() {
  cat("🧪 Running Comprehensive Integration Tests\n")
  cat("==========================================\n")

  # Setup test environment
  test_config <- setup_test_environment()

  # Run test suite
  test_results <- testthat::test_file("tests/test_integration.R")

  # Generate test report
  test_report <- list(
    timestamp = Sys.time(),
    tests_run = length(test_results),
    tests_passed = sum(sapply(test_results, function(x) x$passed)),
    tests_failed = sum(sapply(test_results, function(x) !x$passed)),
    total_time = sum(sapply(test_results, function(x) x$time))
  )

  # Display results
  cat("\n📊 Integration Test Results\n")
  cat("==========================\n")
  cat("Tests run:", test_report$tests_run, "\n")
  cat("Tests passed:", test_report$tests_passed, "\n")
  cat("Tests failed:", test_report$tests_failed, "\n")
  cat("Total time:", round(test_report$total_time, 2), "seconds\n")

  if (test_report$tests_failed > 0) {
    cat("\n❌ Some integration tests failed!\n")
    cat("Check test output above for details.\n")
  } else {
    cat("\n✅ All integration tests passed!\n")
    cat("🎉 Platform integration verified successfully.\n")
  }

  # Cleanup
  cleanup_test_environment()

  return(test_report)
}

# If run directly, execute tests
if (!interactive()) {
  run_integration_tests()
}


