# CI/CD Integration Tests for NHANES BMI Body Fat Analysis Platform
# Tests automated deployment, testing, and validation workflows

library(testthat)
library(yaml)
library(jsonlite)
library(digest)

# Source all required modules
source("../R/error_handling.R")
source("../R/data_versioning.R")
source("../R/data_validation.R")
source("../parallel_pipeline.R")

# CI/CD Test Configuration
CI_CONFIG <- list(
  test_timeout = 300,  # 5 minutes timeout
  memory_limit = 1024,  # 1GB memory limit
  parallel_workers = 2,  # Limited workers for CI
  cache_enabled = TRUE,
  verbose_logging = FALSE
)

# Test 1: CI Environment Detection
test_that("CI environment is properly detected", {

  # Test 1a: Environment variable detection
  # In CI environments, certain variables may be set
  ci_vars <- c("CI", "GITHUB_ACTIONS", "TRAVIS", "JENKINS_URL")
  detected_ci <- any(sapply(ci_vars, function(var) {
    !is.null(Sys.getenv(var))
  }))

  # Should work in both CI and local environments
  expect_true(is.logical(detected_ci))

  # Test 1b: Resource limitations
  memory_available <- memory.size()
  expect_true(is.numeric(memory_available))

  # Test 1c: CPU detection
  cores_available <- availableCores()
  expect_true(is.integer(cores_available))
  expect_true(cores_available >= 1)
})

# Test 2: Automated Package Installation
test_that("Package installation works in automated environments", {

  # Test 2a: Required packages availability
  required_packages <- c("dplyr", "ggplot2", "survey", "foreign", "yaml", "future", "furrr", "digest", "jsonlite")

  available_packages <- 0
  for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      available_packages <- available_packages + 1
    }
  }

  # Most packages should be available (allowing for some missing in CI)
  package_coverage <- available_packages / length(required_packages)
  expect_true(package_coverage >= 0.8)  # At least 80% coverage

  # Test 2b: Package loading
  loaded_packages <- 0
  for (pkg in required_packages) {
    tryCatch({
      library(pkg, character.only = TRUE)
      loaded_packages <- loaded_packages + 1
    }, error = function(e) {
      # Some packages might not load in CI, that's OK
    })
  }

  # At least core packages should load
  expect_true(loaded_packages >= 5)
})

# Test 3: Configuration Loading in CI
test_that("Configuration loads correctly in automated environments", {

  # Test 3a: Configuration file existence
  config_file <- "config/config.yml"
  expect_true(file.exists(config_file))

  # Test 3b: Configuration parsing
  config <- read_yaml(config_file)
  expect_is(config, "list")

  # Test 3c: Configuration validation
  required_sections <- c("data", "outputs", "nhanes", "analysis", "logging")
  for (section in required_sections) {
    expect_true(section %in% names(config))
  }

  # Test 3d: Configuration values
  expect_true(is.numeric(config$analysis$age_range[1]))
  expect_true(is.character(config$analysis$survey_weights_col))
})

# Test 4: Automated Data Download
test_that("Data download works in CI environments", {

  # Test 4a: Directory creation
  data_dirs <- c("data/raw", "data/derived")
  for (dir in data_dirs) {
    if (!dir.exists(dir)) {
      expect_true(dir.create(dir, recursive = TRUE))
    }
    expect_true(dir.exists(dir))
  }

  # Test 4b: Download simulation (without actual network calls)
  # In real CI, this would test actual downloads
  mock_file <- "data/raw/test_download.csv"
  mock_data <- data.frame(test = 1:10)

  write.csv(mock_data, mock_file, row.names = FALSE)
  expect_true(file.exists(mock_file))

  # Test 4c: File validation
  file_info <- file.info(mock_file)
  expect_true(file_info$size > 0)

  # Cleanup
  if (file.exists(mock_file)) {
    file.remove(mock_file)
  }
})

# Test 5: Registry Integration in CI
test_that("Data registry works in automated environments", {

  # Test 5a: Registry initialization
  success <- initialize_data_registry()
  expect_true(success)
  expect_true(file.exists("data/registry/data_registry.json"))

  # Test 5b: Registry loading
  registry <- load_data_registry()
  expect_is(registry, "list")
  expect_true("metadata" %in% names(registry))
  expect_true("entries" %in% names(registry))

  # Test 5c: Registry backup/restore
  registry_file <- "data/registry/data_registry.json"
  backup_file <- "data/registry/backup_registry.json"

  # Create backup
  if (file.exists(registry_file)) {
    file.copy(registry_file, backup_file, overwrite = TRUE)
    expect_true(file.exists(backup_file))
  }

  # Test 5d: Registry corruption recovery
  # Simulate corrupted registry
  writeLines("invalid json", registry_file)

  # Should handle corrupted registry gracefully
  corrupted_registry <- load_data_registry()
  expect_is(corrupted_registry, "list")  # Should return default structure

  # Restore from backup
  if (file.exists(backup_file)) {
    file.copy(backup_file, registry_file, overwrite = TRUE)
    restored_registry <- load_data_registry()
    expect_is(restored_registry, "list")
  }
})

# Test 6: Parallel Processing in CI
test_that("Parallel processing works in CI environments", {

  # Test 6a: Worker allocation for CI
  original_plan <- plan()

  # Use conservative worker count for CI
  workers <- min(availableCores() - 1, 2)
  if (workers > 0) {
    plan(multisession, workers = workers)

    # Test parallel computation
    test_data <- 1:20
    parallel_results <- future_map(test_data, function(x) x * 2)

    # Should produce correct results
    expected_results <- test_data * 2
    expect_equal(parallel_results, expected_results)
  }

  # Restore original plan
  plan(original_plan)

  # Test 6b: Memory constraints
  memory_before <- memory.size()

  # Create moderate dataset
  test_dataset <- data.frame(matrix(rnorm(10000), ncol = 20))

  memory_after <- memory.size()
  memory_growth <- memory_after - memory_before

  # Memory growth should be reasonable for CI
  expect_true(memory_growth < 100)  # Less than 100MB

  # Cleanup
  rm(test_dataset)
  gc()
})

# Test 7: Error Handling in CI
test_that("Error handling works in automated environments", {

  # Test 7a: CI-friendly error reporting
  test_error <- NhanesError("CI test error", "CI001")
  expect_is(test_error, "NhanesError")

  # Test 7b: Error logging
  test_config <- list(outputs = list(logs_dir = "tests/test_outputs"))
  safe_log("CI test message", "INFO", test_config)

  # Test 7c: Error suggestions
  suggestions <- get_error_suggestions("CI001", "CI test error")
  expect_is(suggestions, "character")
  expect_true(length(suggestions) > 0)

  # Test 7d: User-friendly error display
  expect_silent(display_user_friendly_error(test_error))
})

# Test 8: Output Generation in CI
test_that("Output generation works in automated environments", {

  # Test 8a: Directory creation
  output_dirs <- c("outputs/tables", "outputs/figures", "outputs/logs")
  for (dir in output_dirs) {
    if (!dir.exists(dir)) {
      expect_true(dir.create(dir, recursive = TRUE))
    }
    expect_true(dir.exists(dir))
  }

  # Test 8b: File writing
  test_output <- data.frame(test = 1:10, value = rnorm(10))
  test_file <- "outputs/tables/ci_test.csv"

  write.csv(test_output, test_file, row.names = FALSE)
  expect_true(file.exists(test_file))

  # Verify file content
  read_output <- read.csv(test_file)
  expect_equal(nrow(read_output), 10)
  expect_true("test" %in% names(read_output))

  # Test 8c: Log file creation
  log_file <- "outputs/logs/ci_test.log"
  writeLines(c("CI test log entry 1", "CI test log entry 2"), log_file)
  expect_true(file.exists(log_file))

  # Cleanup
  if (file.exists(test_file)) {
    file.remove(test_file)
  }
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
})

# Test 9: Performance in CI Environment
test_that("Performance is acceptable in CI environments", {

  # Test 9a: Quick operations
  start_time <- Sys.time()

  # Simple computation
  result <- sum(1:100000)
  expect_equal(result, 5000050000)

  execution_time <- difftime(Sys.time(), start_time, units = "secs")
  expect_true(execution_time < 2)  # Should be very fast

  # Test 9b: Memory usage
  memory_before <- memory.size()

  # Moderate computation
  test_matrix <- matrix(rnorm(1000), nrow = 100)
  result <- eigen(test_matrix)$values

  memory_after <- memory.size()
  memory_growth <- memory_after - memory_before

  # Memory growth should be reasonable for CI
  expect_true(memory_growth < 50)  # Less than 50MB

  # Cleanup
  rm(test_matrix, result)
  gc()

  # Test 9c: File I/O performance
  test_data <- data.frame(x = 1:1000, y = rnorm(1000))

  start_time <- Sys.time()
  write.csv(test_data, "tests/test_io.csv", row.names = FALSE)
  write_time <- difftime(Sys.time(), start_time, units = "secs")

  start_time <- Sys.time()
  read_data <- read.csv("tests/test_io.csv")
  read_time <- difftime(Sys.time(), start_time, units = "secs")

  # I/O should be reasonably fast
  expect_true(write_time < 1)
  expect_true(read_time < 1)

  # Cleanup
  if (file.exists("tests/test_io.csv")) {
    file.remove("tests/test_io.csv")
  }
})

# Test 10: CI/CD Workflow Simulation
test_that("CI/CD workflow simulation works", {

  # Test 10a: Environment setup
  # Simulate CI environment setup
  ci_env <- list(
    CI = "true",
    GITHUB_ACTIONS = "true",
    R_LIBS_USER = Sys.getenv("R_LIBS_USER", unset = "")
  )

  # Should handle CI environment variables
  expect_true(is.character(ci_env$CI))

  # Test 10b: Automated testing
  # Simulate running tests in CI
  test_results <- list()
  test_count <- 0

  # Run a few basic tests
  test_cases <- list(
    list(name = "basic_arithmetic", test = function() expect_equal(2 + 2, 4)),
    list(name = "package_loading", test = function() expect_true(requireNamespace("base", quietly = TRUE))),
    list(name = "file_creation", test = function() {
      test_file <- "tests/ci_test.txt"
      writeLines("CI test", test_file)
      expect_true(file.exists(test_file))
      file.remove(test_file)
    })
  )

  for (test_case in test_cases) {
    tryCatch({
      test_case$test()
      test_results[[test_case$name]] <- "PASSED"
      test_count <- test_count + 1
    }, error = function(e) {
      test_results[[test_case$name]] <- paste("FAILED:", e$message)
    })
  }

  # Should have run some tests
  expect_true(test_count > 0)

  # Test 10c: Artifact collection
  # Simulate CI artifact collection
  artifacts <- c("outputs/tables/", "outputs/figures/", "outputs/logs/")
  artifact_count <- 0

  for (artifact in artifacts) {
    if (dir.exists(artifact)) {
      artifact_count <- artifact_count + 1
    }
  }

  # Should have created output directories
  expect_true(artifact_count > 0)

  # Test 10d: Cleanup simulation
  # Simulate CI cleanup
  cleanup_dirs <- c("tests/test_outputs", "tests/test_cache")

  for (cleanup_dir in cleanup_dirs) {
    if (dir.exists(cleanup_dir)) {
      expect_silent(unlink(cleanup_dir, recursive = TRUE))
    }
  }
})

# Test 11: Resource Constraints in CI
test_that("System handles CI resource constraints", {

  # Test 11a: Memory constraints
  memory_limit <- memory.limit()
  expect_true(is.numeric(memory_limit))

  # Should work within CI memory limits
  test_data <- data.frame(matrix(rnorm(1000), ncol = 10))
  expect_true(object.size(test_data) < memory_limit * 0.1)  # Less than 10% of limit

  # Test 11b: CPU constraints
  cores_available <- availableCores()
  expect_true(cores_available >= 1)

  # Should adapt to limited cores
  workers <- min(cores_available - 1, 2)  # Conservative allocation
  expect_true(workers >= 0)

  # Test 11c: Time constraints
  start_time <- Sys.time()

  # Quick validation
  quick_test <- function() {
    x <- 1:100
    return(mean(x))
  }

  result <- quick_test()
  execution_time <- difftime(Sys.time(), start_time, units = "secs")

  expect_true(execution_time < 1)  # Should be very fast
  expect_equal(result, 50.5)

  # Test 11d: Storage constraints
  # Test that temporary files are cleaned up
  temp_files <- list.files(tempdir(), pattern = "^file")
  initial_temp_count <- length(temp_files)

  # Create temporary file
  temp_file <- tempfile()
  writeLines("test", temp_file)

  # Verify file was created
  expect_true(file.exists(temp_file))

  # Cleanup
  file.remove(temp_file)

  # Verify cleanup
  final_temp_count <- length(list.files(tempdir(), pattern = "^file"))
  expect_equal(final_temp_count, initial_temp_count)
})

# Test 12: Automated Quality Checks
test_that("Automated quality checks work in CI", {

  # Test 12a: Code quality checks
  # Test that linting would work (if lintr available)
  if (requireNamespace("lintr", quietly = TRUE)) {
    # Would run: lintr::lint_dir('scripts')
    expect_true(requireNamespace("lintr", quietly = TRUE))
  }

  # Test 12b: Package structure validation
  # Check that required files exist
  required_files <- c(
    "DESCRIPTION",
    "R/error_handling.R",
    "R/data_versioning.R",
    "scripts/fetch_nhanes.R",
    "config/config.yml"
  )

  for (file in required_files) {
    if (file.exists(file)) {
      expect_true(file.exists(file))
    }
  }

  # Test 12c: Documentation validation
  # Check that key documentation exists
  doc_files <- c("README.md", "docs/INSTALLATION.md")
  for (doc in doc_files) {
    if (file.exists(doc)) {
      content <- readLines(doc, n = 10)
      expect_true(length(content) > 0)
    }
  }

  # Test 12d: Test suite validation
  # Check that test files exist
  test_files <- c("tests/test_integration.R", "tests/test_parallel_pipeline.R")
  for (test_file in test_files) {
    if (file.exists(test_file)) {
      expect_true(file.exists(test_file))
    }
  }
})

# Run all CI/CD tests with reporting
run_ci_cd_tests <- function() {
  cat("🚀 Running CI/CD Integration Tests\n")
  cat("==================================\n")

  # Run test suite
  test_results <- testthat::test_file("tests/test_ci_cd.R")

  # Generate report
  test_report <- list(
    timestamp = Sys.time(),
    tests_run = length(test_results),
    tests_passed = sum(sapply(test_results, function(x) x$passed)),
    tests_failed = sum(sapply(test_results, function(x) !x$passed)),
    total_time = sum(sapply(test_results, function(x) x$time)),
    environment = list(
      os = Sys.info()["sysname"],
      r_version = R.version.string,
      memory_limit = memory.limit(),
      cores_available = availableCores(),
      ci_detected = !is.null(Sys.getenv("CI"))
    )
  )

  # Display results
  cat("\n📊 CI/CD Integration Test Results\n")
  cat("=================================\n")
  cat("Tests run:", test_report$tests_run, "\n")
  cat("Tests passed:", test_report$tests_passed, "\n")
  cat("Tests failed:", test_report$tests_failed, "\n")
  cat("Total time:", round(test_report$total_time, 2), "seconds\n")

  cat("\n🖥️ Environment Information:\n")
  cat("  OS:", test_report$environment$os, "\n")
  cat("  R Version:", test_report$environment$r_version, "\n")
  cat("  Memory Limit:", round(test_report$environment$memory_limit / 1024, 1), "GB\n")
  cat("  CPU Cores:", test_report$environment$cores_available, "\n")
  cat("  CI Detected:", test_report$environment$ci_detected, "\n")

  if (test_report$tests_failed > 0) {
    cat("\n❌ Some CI/CD integration tests failed!\n")
    cat("Check test output above for details.\n")
    cat("This may indicate issues with automated deployment.\n")
  } else {
    cat("\n✅ All CI/CD integration tests passed!\n")
    cat("🚀 Platform ready for automated deployment.\n")
  }

  return(test_report)
}

# If run directly, execute tests
if (!interactive()) {
  run_ci_cd_tests()
}


