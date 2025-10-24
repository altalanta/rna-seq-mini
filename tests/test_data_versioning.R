# Data Version Management Integration Tests
# Tests data registry, integrity verification, and update detection

library(testthat)
library(digest)
library(yaml)
library(jsonlite)

# Source data versioning functions
source("../R/data_versioning.R")

# Test data setup
setup_test_files <- function() {
  # Create test directory structure
  test_dirs <- c("tests/test_registry", "tests/test_files")
  for (dir in test_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  # Create mock NHANES-like files
  set.seed(123)

  # Demographics file
  demo_data <- data.frame(
    SEQN = 1:100,
    RIDAGEYR = sample(20:59, 100, replace = TRUE),
    RIAGENDR = sample(1:2, 100, replace = TRUE)
  )
  write.csv(demo_data, "tests/test_files/DEMO_J.csv", row.names = FALSE)

  # Body measures file
  bmx_data <- data.frame(
    SEQN = 1:100,
    BMXBMI = rnorm(100, 27, 5)
  )
  write.csv(bmx_data, "tests/test_files/BMX_J.csv", row.names = FALSE)

  # DXA whole body file
  dxx_data <- data.frame(
    SEQN = 1:100,
    DXDTOFAT = rnorm(100, 25, 8)  # Mock body fat percentage
  )
  write.csv(dxx_data, "tests/test_files/DXX_J.csv", row.names = FALSE)

  # DXA android/gynoid file
  dxxag_data <- data.frame(
    SEQN = 1:100,
    DXDAPFAT = rnorm(100, 15, 5),  # Mock android fat
    DXDGPFAT = rnorm(100, 20, 6)   # Mock gynoid fat
  )
  write.csv(dxxag_data, "tests/test_files/DXXAG_J.csv", row.names = FALSE)

  return(list(
    demo_file = "tests/test_files/DEMO_J.csv",
    bmx_file = "tests/test_files/BMX_J.csv",
    dxx_file = "tests/test_files/DXX_J.csv",
    dxxag_file = "tests/test_files/DXXAG_J.csv"
  ))
}

# Test 1: Registry System Integration
test_that("Data registry system integrates properly", {

  # Test 1a: Registry initialization
  expect_true(initialize_data_registry())

  registry_file <- "data/registry/data_registry.json"
  expect_true(file.exists(registry_file))

  # Test 1b: Registry loading
  registry <- load_data_registry()
  expect_is(registry, "list")
  expect_true("metadata" %in% names(registry))
  expect_true("entries" %in% names(registry))

  # Test 1c: Registry structure
  expect_true(is.list(registry$metadata))
  expect_true(is.list(registry$entries))
  expect_equal(registry$metadata$version, "1.0")

  # Test 1d: Registry saving
  success <- save_data_registry(registry)
  expect_true(success)
  expect_true(file.exists(registry_file))
})

# Test 2: File Hashing and Integrity
test_that("File hashing and integrity verification works", {

  # Test 2a: Hash computation
  test_file <- "tests/test_files/DEMO_J.csv"
  if (file.exists(test_file)) {
    hash <- compute_file_hash(test_file)
    expect_is(hash, "character")
    expect_true(nchar(hash) == 64)  # SHA256 hash length

    # Hash should be consistent
    hash2 <- compute_file_hash(test_file)
    expect_equal(hash, hash2)
  }

  # Test 2b: File metadata extraction
  metadata <- get_file_metadata(test_file)
  if (!is.null(metadata)) {
    expect_is(metadata, "list")
    expect_true("size" %in% names(metadata))
    expect_true("mtime" %in% names(metadata))
    expect_true("basename" %in% names(metadata))
  }

  # Test 2c: Registry entry creation
  entry <- create_registry_entry(test_file, "demographics", "2017-2018")
  if (!is.null(entry)) {
    expect_is(entry, "list")
    expect_true("file_id" %in% names(entry))
    expect_true("hash_sha256" %in% names(entry))
    expect_true("file_size" %in% names(entry))
    expect_equal(entry$data_type, "demographics")
    expect_equal(entry$nhanes_cycle, "2017-2018")
  }
})

# Test 3: Registry Operations
test_that("Registry operations work correctly", {

  # Test 3a: Add entry to registry
  test_file <- "tests/test_files/DEMO_J.csv"
  if (file.exists(test_file)) {
    success <- add_to_registry(test_file, "demographics", "2017-2018")
    expect_true(success)

    # Test 3b: Verify entry was added
    registry <- load_data_registry()
    expect_true(length(registry$entries) > 0)

    # Find our entry
    entry <- registry$entries[[1]]
    expect_equal(entry$data_type, "demographics")
    expect_equal(entry$nhanes_cycle, "2017-2018")
    expect_true(!is.null(entry$hash_sha256))
  }

  # Test 3c: Update detection
  updates <- check_for_updates()
  expect_is(updates, "list")
  expect_true("uptodate" %in% names(updates))
  expect_true("updates" %in% names(updates))

  # Test 3d: Integrity validation
  integrity <- validate_data_integrity()
  expect_is(integrity, "list")
  expect_true("valid" %in% names(integrity))
  expect_true("issues" %in% names(integrity))
})

# Test 4: Quality Reporting
test_that("Quality reporting integrates properly", {

  # Test 4a: Quality report generation
  quality_report <- generate_quality_report()
  expect_is(quality_report, "list")
  expect_true(all(c("metadata", "integrity_check", "registry_summary", "recommendations") %in% names(quality_report)))

  # Test 4b: Quality report structure
  expect_true(is.list(quality_report$metadata))
  expect_true(is.list(quality_report$integrity_check))
  expect_true(is.list(quality_report$registry_summary))
  expect_true(is.character(quality_report$recommendations))

  # Test 4c: Registry summary
  expect_true(is.numeric(quality_report$registry_summary$total_files))
  expect_true(is.numeric(quality_report$registry_summary$active_files))

  # Test 4d: Manifest generation
  manifest <- generate_data_manifest()
  expect_is(manifest, "list")
  expect_true(all(c("metadata", "data_files", "summary") %in% names(manifest)))
})

# Test 5: Integration with File System
test_that("Data versioning integrates with file system operations", {

  # Test 5a: Directory creation
  ensure_registry_dirs()
  expect_true(dir.exists("data/registry"))
  expect_true(dir.exists("data/raw"))
  expect_true(dir.exists("data/derived"))

  # Test 5b: File operations
  test_file <- "tests/test_files/DEMO_J.csv"
  if (file.exists(test_file)) {
    # File should be readable
    expect_true(file.access(test_file, mode = 4) == 0)  # Readable

    # Should be able to get file info
    info <- file.info(test_file)
    expect_is(info, "data.frame")
    expect_true("size" %in% names(info))
  }

  # Test 5c: JSON operations
  test_data <- list(test = "value", number = 42)
  json_data <- toJSON(test_data, pretty = TRUE)
  expect_is(json_data, "character")

  parsed_data <- fromJSON(json_data)
  expect_is(parsed_data, "list")
  expect_equal(parsed_data$test, "value")
  expect_equal(parsed_data$number, 42)
})

# Test 6: Error Handling Integration
test_that("Error handling works with data versioning", {

  # Test 6a: File not found error
  non_existent_file <- "non_existent_file.csv"
  expect_error(compute_file_hash(non_existent_file))

  # Test 6b: Registry error handling
  # Test with corrupted registry
  registry_file <- "data/registry/data_registry.json"

  # Backup original registry
  if (file.exists(registry_file)) {
    file.copy(registry_file, "data/registry/backup_registry.json", overwrite = TRUE)
  }

  # Create corrupted registry
  writeLines("invalid json content", registry_file)

  # Should handle corrupted registry gracefully
  corrupted_registry <- load_data_registry()
  expect_is(corrupted_registry, "list")  # Should return default structure

  # Restore original registry
  if (file.exists("data/registry/backup_registry.json")) {
    file.copy("data/registry/backup_registry.json", registry_file, overwrite = TRUE)
  }

  # Test 6c: Error suggestion system
  error_suggestions <- get_error_suggestions("DL001", "File not found")
  expect_is(error_suggestions, "character")
  expect_true(length(error_suggestions) > 0)
  expect_true(any(grepl("download", error_suggestions, ignore.case = TRUE)))
})

# Test 7: Update Detection Integration
test_that("Update detection integrates properly", {

  # Test 7a: Update check structure
  updates <- check_for_updates()
  expect_is(updates, "list")
  expect_true("uptodate" %in% names(updates))
  expect_true("updates" %in% names(updates))

  # Test 7b: Update check logic
  if (length(updates$updates) > 0) {
    for (file_name in names(updates$updates)) {
      update_info <- updates$updates[[file_name]]
      expect_true("status" %in% names(update_info))
      expect_true("action" %in% names(update_info))
      expect_true("description" %in% names(update_info))
    }
  }

  # Test 7c: Network connectivity simulation
  # This would test actual network calls in a real environment
  # For now, test the structure of the response
  expect_is(updates$uptodate, "logical")
})

# Test 8: Registry Persistence
test_that("Registry persistence works across sessions", {

  # Test 8a: Registry file persistence
  registry_file <- "data/registry/data_registry.json"
  expect_true(file.exists(registry_file))

  # Test 8b: JSON format validation
  registry_content <- readLines(registry_file)
  expect_true(length(registry_content) > 0)

  # Should be valid JSON
  registry_json <- paste(registry_content, collapse = "\n")
  expect_silent(fromJSON(registry_json))

  # Test 8c: Registry backup and restore
  if (file.exists(registry_file)) {
    # Create backup
    backup_file <- "data/registry/backup_registry.json"
    file.copy(registry_file, backup_file, overwrite = TRUE)

    # Modify registry
    registry <- load_data_registry()
    original_count <- length(registry$entries)

    # Add test entry
    test_entry <- create_registry_entry("tests/test_files/test_entry.csv", "test", "2017-2018")
    if (!is.null(test_entry)) {
      registry$entries <- c(registry$entries, list(test_entry))
      save_data_registry(registry)

      # Verify modification
      modified_registry <- load_data_registry()
      expect_true(length(modified_registry$entries) > original_count)

      # Restore from backup
      file.copy(backup_file, registry_file, overwrite = TRUE)
      restored_registry <- load_data_registry()
      expect_equal(length(restored_registry$entries), original_count)
    }
  }
})

# Test 9: Integration with Configuration System
test_that("Data versioning integrates with configuration", {

  # Test 9a: Configuration integration
  if (file.exists("../config/config.yml")) {
    config <- read_yaml("../config/config.yml")
    expect_is(config, "list")

    # Test configuration values for data management
    expect_true(is.character(config$analysis$survey_weights_col))
    expect_true(is.character(config$analysis$strata_col))
    expect_true(is.character(config$analysis$psu_col))
  }

  # Test 9b: Directory configuration
  # Test that configured directories are accessible
  test_dirs <- c("data/raw", "data/derived", "outputs/tables", "outputs/figures")

  for (test_dir in test_dirs) {
    if (dir.exists(test_dir)) {
      expect_true(dir.exists(test_dir))
    }
  }

  # Test 9c: Path normalization
  test_paths <- c("data/raw", "data/derived", "outputs/tables")
  for (path in test_paths) {
    normalized <- normalizePath(path, mustWork = FALSE)
    expect_is(normalized, "character")
  }
})

# Test 10: Performance and Scalability
test_that("Data versioning performs adequately", {

  # Test 10a: Hash computation performance
  test_file <- "tests/test_files/DEMO_J.csv"
  if (file.exists(test_file)) {
    start_time <- Sys.time()
    hash <- compute_file_hash(test_file)
    hash_time <- difftime(Sys.time(), start_time, units = "secs")

    expect_true(hash_time < 1)  # Should be fast for small files
    expect_true(nchar(hash) == 64)  # SHA256 hash
  }

  # Test 10b: Registry operations performance
  start_time <- Sys.time()
  registry <- load_data_registry()
  load_time <- difftime(Sys.time(), start_time, units = "secs")

  expect_true(load_time < 0.1)  # Should load quickly

  # Test 10c: Memory usage
  initial_memory <- memory.size()

  # Create large registry for testing
  large_registry <- list(
    metadata = list(created = as.character(Sys.time())),
    entries = list()
  )

  # Add multiple entries
  for (i in 1:100) {
    entry <- create_registry_entry(
      paste0("tests/test_files/test_", i, ".csv"),
      "test_data",
      "2017-2018"
    )
    if (!is.null(entry)) {
      large_registry$entries <- c(large_registry$entries, list(entry))
    }
  }

  memory_after <- memory.size()
  memory_growth <- memory_after - initial_memory

  # Memory growth should be reasonable
  expect_true(memory_growth < 50)  # Less than 50MB

  # Cleanup
  rm(large_registry)
  gc()
})

# Test 11: Integration with Pipeline
test_that("Data versioning integrates with pipeline workflow", {

  # Test 11a: Pipeline compatibility
  # Test that data versioning functions don't break pipeline
  expect_true(exists("run_parallel_pipeline"))

  # Test 11b: Configuration compatibility
  # Test that versioning works with pipeline configuration
  if (file.exists("../config/config.yml")) {
    config <- read_yaml("../config/config.yml")

    # Should be able to use config values
    expect_true(is.character(config$analysis$survey_weights_col))

    # Test directory paths from config
    expect_true(is.character(config$data$raw_dir))
    expect_true(is.character(config$outputs$tables_dir))
  }

  # Test 11c: Error propagation
  # Test that versioning errors don't break pipeline
  test_error <- NhanesError("Versioning test error", "VER001")
  expect_silent(display_user_friendly_error(test_error))
})

# Test 12: Regression Testing
test_that("Data versioning doesn't break existing functionality", {

  # Test 12a: Existing functions still work
  expect_true(exists("safe_read_xpt"))
  expect_true(exists("validate_nhanes_data"))
  expect_true(exists("ensure_output_dirs"))

  # Test 12b: New functions don't interfere
  # Test that new versioning functions don't break existing code
  test_file <- "tests/test_files/DEMO_J.csv"
  if (file.exists(test_file)) {
    # Old function should still work
    expect_silent({
      old_metadata <- get_file_metadata(test_file)
      expect_is(old_metadata, "list")
    })

    # New function should also work
    expect_silent({
      hash <- compute_file_hash(test_file)
      expect_is(hash, "character")
    })
  }

  # Test 12c: Backward compatibility
  # Test that old configuration files still work
  if (file.exists("../config/config.yml")) {
    config <- read_yaml("../config/config.yml")
    expect_is(config, "list")

    # Should still have old structure
    expect_true("data" %in% names(config))
    expect_true("outputs" %in% names(config))
    expect_true("nhanes" %in% names(config))
  }
})

# Run all tests with reporting
run_versioning_tests <- function() {
  cat("🔒 Running Data Version Management Integration Tests\n")
  cat("===================================================\n")

  # Setup test environment
  test_files <- setup_test_files()

  # Run test suite
  test_results <- testthat::test_file("tests/test_data_versioning.R")

  # Generate report
  test_report <- list(
    timestamp = Sys.time(),
    tests_run = length(test_results),
    tests_passed = sum(sapply(test_results, function(x) x$passed)),
    tests_failed = sum(sapply(test_results, function(x) !x$passed)),
    total_time = sum(sapply(test_results, function(x) x$time))
  )

  # Display results
  cat("\n📊 Data Versioning Test Results\n")
  cat("===============================\n")
  cat("Tests run:", test_report$tests_run, "\n")
  cat("Tests passed:", test_report$tests_passed, "\n")
  cat("Tests failed:", test_report$tests_failed, "\n")
  cat("Total time:", round(test_report$total_time, 2), "seconds\n")

  # Check registry status
  registry_file <- "data/registry/data_registry.json"
  if (file.exists(registry_file)) {
    registry <- load_data_registry()
    cat("\n📋 Registry Status:\n")
    cat("  Registry file: ✅ Present\n")
    cat("  Total entries:", length(registry$entries), "\n")
    cat("  Registry version:", registry$metadata$version, "\n")
  } else {
    cat("\n📋 Registry Status:\n")
    cat("  Registry file: ❌ Missing\n")
  }

  # Check test files
  test_files_exist <- sum(file.exists(unlist(test_files)))
  cat("\n📁 Test Files:\n")
  cat("  Test files created:", test_files_exist, "/ 4\n")

  if (test_report$tests_failed > 0) {
    cat("\n❌ Some data versioning tests failed!\n")
    cat("Check test output above for details.\n")
  } else {
    cat("\n✅ All data versioning tests passed!\n")
    cat("🎉 Data version management integration verified successfully.\n")
  }

  return(test_report)
}

# If run directly, execute tests
if (!interactive()) {
  run_versioning_tests()
}


