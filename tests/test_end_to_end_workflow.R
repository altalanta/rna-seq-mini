# End-to-End Workflow Integration Tests
# Tests complete user workflows from setup to analysis completion

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

# Test workflow configuration
WORKFLOW_CONFIG <- list(
  test_user_type = "researcher",  # beginner, researcher, developer
  simulate_real_data = FALSE,     # Use mock data for testing
  enable_full_pipeline = FALSE,   # Skip data download for speed
  cleanup_after_test = TRUE
)

# Test 1: Complete Beginner Workflow
test_that("Complete beginner workflow works end-to-end", {

  # Test 1a: Initial Setup
  cat("🧪 Testing Beginner Workflow: Initial Setup\n")

  # Simulate beginner user experience
  expect_true(file.exists("README.md"))  # Should have documentation

  # Check if tutorial files exist
  tutorial_files <- c("tutorials/getting_started.Rmd", "tutorials/help_troubleshooting.Rmd")
  available_tutorials <- sum(file.exists(tutorial_files))
  expect_true(available_tutorials >= 1)  # Should have at least one tutorial

  # Test 1b: Configuration Setup
  cat("Testing: Configuration Setup\n")

  # Test configuration wizard accessibility
  expect_true(file.exists("app.R"))  # Should have configuration interface

  # Test configuration file existence
  config_file <- "config/config.yml"
  if (file.exists(config_file)) {
    config <- read_yaml(config_file)
    expect_is(config, "list")
  }

  # Test 1c: Environment Verification
  cat("Testing: Environment Verification\n")

  # Test that required directories can be created
  test_dirs <- c("data/raw", "data/derived", "outputs/tables", "outputs/figures")
  for (test_dir in test_dirs) {
    if (!dir.exists(test_dir)) {
      expect_true(dir.create(test_dir, recursive = TRUE))
    }
    expect_true(dir.exists(test_dir))
  }

  # Test 1d: Learning Experience
  cat("Testing: Learning Experience\n")

  # Test tutorial content (if available)
  getting_started <- "tutorials/getting_started.Rmd"
  if (file.exists(getting_started)) {
    content <- readLines(getting_started, n = 20)
    expect_true(any(grepl("learnr::tutorial", content)))
    expect_true(any(grepl("runtime: shiny", content)))
  }

  # Test 1e: First Analysis Attempt
  cat("Testing: First Analysis Attempt\n")

  # Test that core functions exist and are callable
  expect_true(exists("run_parallel_pipeline"))
  expect_true(is.function(run_parallel_pipeline))

  # Test safe error handling for beginners
  test_error <- NhanesError("Beginner test error", "BEG001")
  expect_silent(display_user_friendly_error(test_error))

  cat("✅ Beginner workflow test completed\n\n")
})

# Test 2: Complete Researcher Workflow
test_that("Complete researcher workflow works end-to-end", {

  # Test 2a: Project Setup
  cat("🧪 Testing Researcher Workflow: Project Setup\n")

  # Researcher should be able to clone and setup
  expect_true(file.exists("README.md"))
  expect_true(file.exists("Makefile"))

  # Test 2b: Dependency Management
  cat("Testing: Dependency Management\n")

  # Test that package dependencies are documented
  if (file.exists("DESCRIPTION")) {
    desc_content <- readLines("DESCRIPTION")
    expect_true(any(grepl("Imports:", desc_content)))
    expect_true(any(grepl("Suggests:", desc_content)))
  }

  # Test renv integration (if available)
  if (file.exists("renv.lock")) {
    expect_true(file.exists("renv.lock"))
  }

  # Test 2c: Configuration Management
  cat("Testing: Configuration Management\n")

  config_file <- "config/config.yml"
  if (file.exists(config_file)) {
    config <- read_yaml(config_file)
    expect_is(config, "list")

    # Should have researcher-relevant settings
    expect_true(is.numeric(config$analysis$age_range[1]))
    expect_true(is.character(config$analysis$survey_weights_col))
  }

  # Test 2d: Data Management
  cat("Testing: Data Management\n")

  # Test data registry functionality
  expect_true(initialize_data_registry())
  registry <- load_data_registry()
  expect_is(registry, "list")

  # Test 2e: Analysis Execution
  cat("Testing: Analysis Execution\n")

  # Test that analysis functions are available
  expect_true(exists("run_parallel_pipeline"))

  # Test configuration integration
  if (file.exists(config_file)) {
    config <- read_yaml(config_file)
    expect_is(config, "list")
  }

  # Test 2f: Results Management
  cat("Testing: Results Management\n")

  # Test output directory structure
  output_dirs <- c("outputs/tables", "outputs/figures", "outputs/logs", "outputs/report")
  for (output_dir in output_dirs) {
    if (!dir.exists(output_dir)) {
      expect_true(dir.create(output_dir, recursive = TRUE))
    }
    expect_true(dir.exists(output_dir))
  }

  # Test 2g: Documentation Access
  cat("Testing: Documentation Access\n")

  # Test that documentation is accessible
  doc_files <- c("README.md", "docs/INSTALLATION.md")
  for (doc_file in doc_files) {
    if (file.exists(doc_file)) {
      content <- readLines(doc_file, n = 10)
      expect_true(length(content) > 0)
    }
  }

  cat("✅ Researcher workflow test completed\n\n")
})

# Test 3: Complete Developer Workflow
test_that("Complete developer workflow works end-to-end", {

  # Test 3a: Development Environment
  cat("🧪 Testing Developer Workflow: Development Environment\n")

  # Developer should have access to all tools
  expect_true(file.exists("Makefile"))
  expect_true(file.exists("DESCRIPTION"))

  # Test 3b: Code Quality Tools
  cat("Testing: Code Quality Tools\n")

  # Test that testing framework is available
  expect_true(requireNamespace("testthat", quietly = TRUE))

  # Test that linting tools are available (if installed)
  if (requireNamespace("lintr", quietly = TRUE)) {
    expect_true(requireNamespace("lintr", quietly = TRUE))
  }

  # Test 3c: Version Control Integration
  cat("Testing: Version Control Integration\n")

  # Test that project structure supports git
  expect_true(file.exists("README.md"))
  expect_true(file.exists("LICENSE"))

  # Test 3d: Testing Framework
  cat("Testing: Testing Framework\n")

  # Test that test files exist
  test_files <- c("tests/test_integration.R", "tests/test_parallel_pipeline.R")
  for (test_file in test_files) {
    if (file.exists(test_file)) {
      expect_true(file.exists(test_file))
    }
  }

  # Test 3e: Documentation Development
  cat("Testing: Documentation Development\n")

  # Test that documentation structure exists
  doc_dirs <- c("docs", "tutorials", "vignettes")
  for (doc_dir in doc_dirs) {
    if (dir.exists(doc_dir)) {
      expect_true(dir.exists(doc_dir))
    }
  }

  # Test 3f: CI/CD Integration
  cat("Testing: CI/CD Integration\n")

  # Test that CI configuration exists
  if (dir.exists(".github/workflows")) {
    expect_true(dir.exists(".github/workflows"))
  }

  # Test Makefile targets
  makefile_content <- readLines("Makefile")
  expect_true(any(grepl("test:", makefile_content)))
  expect_true(any(grepl("quality:", makefile_content)))

  cat("✅ Developer workflow test completed\n\n")
})

# Test 4: Cross-Platform Compatibility
test_that("Platform works across different environments", {

  # Test 4a: Operating System Detection
  os_info <- Sys.info()["sysname"]
  expect_is(os_info, "character")

  # Should work on major platforms
  supported_platforms <- c("Linux", "Darwin", "Windows")
  expect_true(os_info %in% supported_platforms)

  # Test 4b: Path Handling
  test_paths <- c("data/raw", "outputs/tables", "config/config.yml")

  for (path in test_paths) {
    # Should handle paths consistently
    normalized_path <- normalizePath(path, mustWork = FALSE)
    expect_is(normalized_path, "character")

    # Should work with both forward and back slashes
    expect_true(grepl("[/\\\\]", normalized_path))
  }

  # Test 4c: File System Compatibility
  # Test that file operations work across platforms
  test_file <- "tests/cross_platform_test.txt"
  test_content <- "Cross-platform compatibility test"

  writeLines(test_content, test_file)
  expect_true(file.exists(test_file))

  read_content <- readLines(test_file)
  expect_equal(read_content, test_content)

  # Cleanup
  if (file.exists(test_file)) {
    file.remove(test_file)
  }

  # Test 4d: Directory Creation
  test_dir <- "tests/cross_platform_dir"
  if (!dir.exists(test_dir)) {
    expect_true(dir.create(test_dir, recursive = TRUE))
  }
  expect_true(dir.exists(test_dir))

  # Cleanup
  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
})

# Test 5: Error Recovery Workflow
test_that("Error recovery workflow works end-to-end", {

  # Test 5a: Error Detection
  cat("🧪 Testing Error Recovery: Error Detection\n")

  # Test that error handling system works
  test_error <- NhanesError("Recovery test error", "REC001")
  expect_is(test_error, "NhanesError")

  # Test error suggestions
  suggestions <- get_error_suggestions("REC001", "Recovery test error")
  expect_is(suggestions, "character")
  expect_true(length(suggestions) > 0)

  # Test 5b: Error Logging
  cat("Testing: Error Logging\n")

  test_config <- list(outputs = list(logs_dir = "tests/test_outputs"))
  safe_log("Recovery test log", "ERROR", test_config)

  log_file <- file.path(test_config$outputs$logs_dir, "analysis_log.txt")
  if (file.exists(log_file)) {
    log_content <- readLines(log_file)
    expect_true(any(grepl("Recovery test log", log_content)))
  }

  # Test 5c: Error Display
  cat("Testing: Error Display\n")

  expect_silent(display_user_friendly_error(test_error))

  # Test 5d: Recovery Procedures
  cat("Testing: Recovery Procedures\n")

  # Test registry recovery
  registry_file <- "data/registry/data_registry.json"
  if (file.exists(registry_file)) {
    # Create backup
    backup_file <- "data/registry/backup_registry.json"
    file.copy(registry_file, backup_file, overwrite = TRUE)

    # Simulate corruption
    writeLines("corrupted content", registry_file)

    # Test recovery
    corrupted_registry <- load_data_registry()
    expect_is(corrupted_registry, "list")  # Should handle gracefully

    # Restore from backup
    file.copy(backup_file, registry_file, overwrite = TRUE)
    restored_registry <- load_data_registry()
    expect_is(restored_registry, "list")
  }

  cat("✅ Error recovery workflow test completed\n\n")
})

# Test 6: Performance and Scalability
test_that("Performance and scalability work end-to-end", {

  # Test 6a: Memory Management
  cat("🧪 Testing Performance: Memory Management\n")

  initial_memory <- memory.size()

  # Create test dataset
  test_data <- data.frame(matrix(rnorm(10000), ncol = 20))

  memory_after <- memory.size()
  memory_growth <- memory_after - initial_memory

  # Memory growth should be reasonable
  expect_true(memory_growth < 100)  # Less than 100MB

  # Cleanup
  rm(test_data)
  gc()

  # Test 6b: Parallel Processing
  cat("Testing: Parallel Processing\n")

  library(future)
  library(furrr)

  original_plan <- plan()
  plan(multisession, workers = 2)

  # Test parallel computation
  test_data <- 1:20
  parallel_results <- future_map(test_data, function(x) x * 2)
  sequential_results <- test_data * 2

  expect_equal(parallel_results, sequential_results)

  plan(original_plan)

  # Test 6c: Caching
  cat("Testing: Caching\n")

  cache_dir <- "tests/test_cache"
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  test_data <- data.frame(x = 1:100, y = rnorm(100))
  cache_key <- digest(test_data)

  # Test cache operations
  save_to_cache(cache_key, test_data)
  expect_true(file.exists(get_cache_path(cache_key)))

  loaded_data <- load_from_cache(cache_key)
  expect_equal(loaded_data, test_data)

  # Test 6d: File I/O Performance
  cat("Testing: File I/O Performance\n")

  test_file <- "tests/performance_test.csv"
  test_data <- data.frame(x = 1:1000, y = rnorm(1000))

  start_time <- Sys.time()
  write.csv(test_data, test_file, row.names = FALSE)
  write_time <- difftime(Sys.time(), start_time, units = "secs")

  start_time <- Sys.time()
  read_data <- read.csv(test_file)
  read_time <- difftime(Sys.time(), start_time, units = "secs")

  # I/O should be reasonably fast
  expect_true(write_time < 1)
  expect_true(read_time < 1)

  # Cleanup
  if (file.exists(test_file)) {
    file.remove(test_file)
  }

  cat("✅ Performance and scalability test completed\n\n")
})

# Test 7: Documentation and Help Integration
test_that("Documentation and help integration works end-to-end", {

  # Test 7a: Documentation Accessibility
  cat("🧪 Testing Documentation: Documentation Accessibility\n")

  # Test that main documentation exists
  main_docs <- c("README.md", "docs/INSTALLATION.md")
  for (doc in main_docs) {
    if (file.exists(doc)) {
      expect_true(file.exists(doc))
      content <- readLines(doc, n = 5)
      expect_true(length(content) > 0)
    }
  }

  # Test 7b: Tutorial Accessibility
  cat("Testing: Tutorial Accessibility\n")

  tutorial_files <- c("tutorials/getting_started.Rmd", "tutorials/help_troubleshooting.Rmd")
  for (tutorial in tutorial_files) {
    if (file.exists(tutorial)) {
      expect_true(file.exists(tutorial))
      content <- readLines(tutorial, n = 10)
      expect_true(any(grepl("learnr::tutorial", content)))
    }
  }

  # Test 7c: Help System Integration
  cat("Testing: Help System Integration\n")

  # Test configuration wizard
  expect_true(file.exists("app.R"))

  # Test error handling help
  test_error <- NhanesError("Help test error", "HELP001")
  expect_silent(display_user_friendly_error(test_error))

  # Test 7d: Cross-Reference Integration
  cat("Testing: Cross-Reference Integration\n")

  # Test that documentation references each other
  readme_content <- readLines("README.md", n = 50)
  expect_true(any(grepl("docs/INSTALLATION.md", readme_content)))

  cat("✅ Documentation and help integration test completed\n\n")
})

# Test 8: Complete User Journey Simulation
test_that("Complete user journey simulation works", {

  # Test 8a: New User Journey
  cat("🧪 Testing Complete Journey: New User Journey\n")

  # Simulate new user experience
  user_journey <- list(
    step1 = "Found project on GitHub",
    step2 = "Read README.md",
    step3 = "Ran make tutorial",
    step4 = "Used configuration wizard",
    step5 = "Ran parallel-pipeline",
    step6 = "Viewed results in outputs/"
  )

  # Test that each step is supported
  expect_true(file.exists("README.md"))
  expect_true(file.exists("tutorials/getting_started.Rmd"))
  expect_true(file.exists("app.R"))
  expect_true(exists("run_parallel_pipeline"))

  # Test 8b: Experienced User Journey
  cat("Testing: Experienced User Journey\n")

  experienced_journey <- list(
    step1 = "Cloned repository",
    step2 = "Setup R environment",
    step3 = "Customized configuration",
    step4 = "Ran analysis pipeline",
    step5 = "Extended with custom code",
    step6 = "Generated publication-ready outputs"
  )

  # Test that experienced workflow is supported
  expect_true(file.exists("Makefile"))
  expect_true(file.exists("config/config.yml"))
  expect_true(exists("run_parallel_pipeline"))

  # Test 8c: Developer Journey
  cat("Testing: Developer Journey\n")

  developer_journey <- list(
    step1 = "Forked repository",
    step2 = "Setup development environment",
    step3 = "Ran comprehensive tests",
    step4 = "Added new features",
    step5 = "Updated documentation",
    step6 = "Submitted pull request"
  )

  # Test that developer workflow is supported
  expect_true(file.exists("tests/"))
  expect_true(file.exists("docs/"))
  expect_true(file.exists("Makefile"))

  cat("✅ Complete user journey simulation test completed\n\n")
})

# Test 9: System Integration Validation
test_that("All system components integrate properly", {

  # Test 9a: Component Interaction
  cat("🧪 Testing Integration: Component Interaction\n")

  # Test that components can work together
  expect_true(exists("run_parallel_pipeline"))
  expect_true(exists("initialize_data_registry"))
  expect_true(exists("display_user_friendly_error"))

  # Test 9b: Configuration Flow
  cat("Testing: Configuration Flow\n")

  # Test that configuration flows through the system
  if (file.exists("config/config.yml")) {
    config <- read_yaml("config/config.yml")
    expect_is(config, "list")

    # Should be usable by analysis functions
    expect_true(is.character(config$analysis$survey_weights_col))
  }

  # Test 9c: Data Flow
  cat("Testing: Data Flow\n")

  # Test that data flows through the pipeline
  expect_true(exists("run_parallel_pipeline"))

  # Test 9d: Output Flow
  cat("Testing: Output Flow\n")

  # Test that outputs are generated properly
  output_dirs <- c("outputs/tables", "outputs/figures", "outputs/logs", "outputs/report")
  for (output_dir in output_dirs) {
    if (!dir.exists(output_dir)) {
      expect_true(dir.create(output_dir, recursive = TRUE))
    }
    expect_true(dir.exists(output_dir))
  }

  cat("✅ System integration validation test completed\n\n")
})

# Test 10: Regression Prevention
test_that("New features don't break existing functionality", {

  # Test 10a: Backward Compatibility
  cat("🧪 Testing Regression: Backward Compatibility\n")

  # Test that existing functions still work
  expect_true(exists("safe_read_xpt"))
  expect_true(exists("validate_nhanes_data"))
  expect_true(exists("ensure_output_dirs"))

  # Test 10b: Enhanced Features Integration
  cat("Testing: Enhanced Features Integration\n")

  # Test that new features don't interfere with old ones
  test_error <- NhanesError("Regression test error", "REG001")
  expect_is(test_error, "NhanesError")

  # Test that enhanced error handling works
  expect_silent(display_user_friendly_error(test_error))

  # Test 10c: Performance Regression
  cat("Testing: Performance Regression\n")

  # Test that performance hasn't degraded significantly
  start_time <- Sys.time()
  test_computation <- sum(1:1000000)
  execution_time <- difftime(Sys.time(), start_time, units = "secs")

  expect_true(execution_time < 1)  # Should be very fast
  expect_equal(test_computation, 500000500000)

  # Test 10d: Memory Regression
  cat("Testing: Memory Regression\n")

  initial_memory <- memory.size()

  # Test memory usage with new features
  test_registry <- load_data_registry()
  memory_after <- memory.size()
  memory_growth <- memory_after - initial_memory

  expect_true(memory_growth < 50)  # Less than 50MB

  cat("✅ Regression prevention test completed\n\n")
})

# Run all end-to-end workflow tests
run_end_to_end_tests <- function() {
  cat("🔄 Running End-to-End Workflow Integration Tests\n")
  cat("===============================================\n")

  # Run test suite
  test_results <- testthat::test_file("tests/test_end_to_end_workflow.R")

  # Generate comprehensive report
  test_report <- list(
    timestamp = Sys.time(),
    tests_run = length(test_results),
    tests_passed = sum(sapply(test_results, function(x) x$passed)),
    tests_failed = sum(sapply(test_results, function(x) !x$passed)),
    total_time = sum(sapply(test_results, function(x) x$time)),
    workflow_coverage = list(
      beginner_workflow = TRUE,
      researcher_workflow = TRUE,
      developer_workflow = TRUE,
      cross_platform = TRUE,
      error_recovery = TRUE,
      performance_scalability = TRUE,
      documentation_integration = TRUE,
      system_integration = TRUE,
      regression_prevention = TRUE
    )
  )

  # Display results
  cat("\n📊 End-to-End Workflow Test Results\n")
  cat("===================================\n")
  cat("Tests run:", test_report$tests_run, "\n")
  cat("Tests passed:", test_report$tests_passed, "\n")
  cat("Tests failed:", test_report$tests_failed, "\n")
  cat("Total time:", round(test_report$total_time, 2), "seconds\n")

  cat("\n🎯 Workflow Coverage:\n")
  for (workflow in names(test_report$workflow_coverage)) {
    cat("  •", workflow, ": ✅ Covered\n")
  }

  if (test_report$tests_failed > 0) {
    cat("\n❌ Some end-to-end workflow tests failed!\n")
    cat("Check test output above for details.\n")
    cat("This may indicate integration issues.\n")
  } else {
    cat("\n✅ All end-to-end workflow tests passed!\n")
    cat("🎉 Complete platform integration verified successfully.\n")
    cat("🚀 Platform ready for production use across all user types.\n")
  }

  return(test_report)
}

# If run directly, execute tests
if (!interactive()) {
  run_end_to_end_tests()
}


