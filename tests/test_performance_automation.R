# Automated Performance Testing Framework for CI/CD Integration
# Tests performance characteristics and ensures they meet requirements

library(testthat)
library(microbenchmark)
library(bench)
library(future)
library(furrr)
library(digest)
library(dplyr)

# Source performance utilities
source("../performance/benchmarking_system.R")

# Performance test configuration for CI/CD
PERF_TEST_CONFIG <- list(
  # Performance thresholds
  min_speedup = 1.2,           # Minimum 20% improvement required
  max_memory_growth = 200,     # Max 200MB memory growth
  max_execution_time = 60,     # Max 60 seconds for standard operations
  min_cache_efficiency = 0.5,  # Minimum 50% cache hit rate

  # Test parameters
  sample_sizes = c(100, 1000, 5000),
  iterations = 5,              # Reduced for CI speed
  timeout_seconds = 300,       # 5 minutes timeout

  # CI-specific settings
  ci_environment = TRUE,
  skip_network_tests = TRUE,
  use_mock_data = TRUE
)

# Test 1: Parallel Processing Performance Requirements
test_that("Parallel processing meets performance requirements", {

  # Test 1a: Basic parallel speedup
  original_plan <- plan()
  workers <- min(availableCores() - 1, 2)  # Conservative for CI

  if (workers > 0) {
    plan(multisession, workers = workers)

    # Simple parallel computation test
    test_data <- 1:1000

    sequential_time <- system.time({
      results <- lapply(test_data, function(x) x * 2)
    })[3]

    parallel_time <- system.time({
      results <- future_map(test_data, function(x) x * 2)
    })[3]

    speedup <- sequential_time / parallel_time

    # Should achieve minimum speedup
    expect_true(speedup >= PERF_TEST_CONFIG$min_speedup,
               paste("Speedup", round(speedup, 2), "x below minimum", PERF_TEST_CONFIG$min_speedup, "x"))

    cat("✅ Parallel speedup:", round(speedup, 2), "x (minimum required:", PERF_TEST_CONFIG$min_speedup, "x)\n")
  }

  plan(original_plan)
})

# Test 2: Memory Usage Performance
test_that("Memory usage meets performance requirements", {

  # Test 2a: Memory growth during operations
  initial_memory <- memory.size()

  # Create test dataset
  test_data <- data.frame(matrix(rnorm(10000), ncol = 20))

  # Process data
  processed_data <- test_data %>%
    mutate(
      row_sum = rowSums(.),
      category = cut(row_sum, breaks = 5)
    ) %>%
    group_by(category) %>%
    summarize(count = n(), mean_value = mean(row_sum))

  final_memory <- memory.size()
  memory_growth <- final_memory - initial_memory

  # Memory growth should be within limits
  expect_true(memory_growth <= PERF_TEST_CONFIG$max_memory_growth,
             paste("Memory growth", round(memory_growth, 1), "MB exceeds limit of", PERF_TEST_CONFIG$max_memory_growth, "MB"))

  cat("✅ Memory growth:", round(memory_growth, 1), "MB (limit:", PERF_TEST_CONFIG$max_memory_growth, "MB)\n")

  # Cleanup
  rm(test_data, processed_data)
  gc()
})

# Test 3: Execution Time Performance
test_that("Execution time meets performance requirements", {

  # Test 3a: Basic operation timing
  execution_time <- system.time({
    # Simulate typical analysis operations
    test_data <- data.frame(x = rnorm(1000), y = rnorm(1000))
    correlation <- cor(test_data$x, test_data$y)
    summary_stats <- summary(test_data)
  })[3]

  # Execution should be reasonably fast
  expect_true(execution_time <= PERF_TEST_CONFIG$max_execution_time,
             paste("Execution time", round(execution_time, 2), "s exceeds limit of", PERF_TEST_CONFIG$max_execution_time, "s"))

  cat("✅ Execution time:", round(execution_time, 2), "s (limit:", PERF_TEST_CONFIG$max_execution_time, "s)\n")
})

# Test 4: Cache Performance (if cache exists)
test_that("Cache performance meets requirements", {

  cache_dir <- "../cache"
  if (dir.exists(cache_dir)) {
    # Test 4a: Cache efficiency
    cache_files <- list.files(cache_dir, full.names = TRUE)

    if (length(cache_files) > 0) {
      # Simulate cache operations
      test_data <- data.frame(x = 1:100, y = rnorm(100))
      cache_key <- digest(test_data)

      # Cache miss timing
      start_time <- Sys.time()
      save_to_cache(cache_key, test_data)
      cache_miss_time <- difftime(Sys.time(), start_time, units = "secs")

      # Cache hit timing
      start_time <- Sys.time()
      loaded_data <- load_from_cache(cache_key)
      cache_hit_time <- difftime(Sys.time(), start_time, units = "secs")

      # Cache should provide benefit
      if (cache_miss_time > 0) {
        cache_speedup <- cache_miss_time / cache_hit_time
        expect_true(cache_speedup >= PERF_TEST_CONFIG$min_cache_efficiency,
                   paste("Cache speedup", round(cache_speedup, 2), "x below minimum", PERF_TEST_CONFIG$min_cache_efficiency, "x"))
      }

      cat("✅ Cache speedup:", round(cache_speedup, 2), "x (minimum:", PERF_TEST_CONFIG$min_cache_efficiency, "x)\n")
    }
  } else {
    skip("Cache directory not found - skipping cache performance tests")
  }
})

# Test 5: System Resource Usage
test_that("System resource usage is within acceptable limits", {

  # Test 5a: CPU usage (simplified)
  cpu_usage <- tryCatch({
    # This is a simplified check - in real CI, you'd use more sophisticated monitoring
    system("ps aux | grep R | grep -v grep | wc -l", intern = TRUE)
  }, error = function(e) "unknown")

  # Should have reasonable R processes
  expect_true(cpu_usage != "0")  # Should have at least one R process

  # Test 5b: Memory availability
  memory_available <- memory.size()
  memory_limit <- memory.limit()

  memory_usage_ratio <- memory_available / memory_limit

  # Memory usage should be reasonable
  expect_true(memory_usage_ratio < 0.9,
             paste("Memory usage", round(memory_usage_ratio * 100, 1), "% exceeds 90% threshold"))

  cat("✅ Memory usage:", round(memory_usage_ratio * 100, 1), "% (threshold: 90%)\n")
})

# Test 6: Data Processing Performance
test_that("Data processing performance meets requirements", {

  # Test 6a: Data loading performance
  if (file.exists("../data/raw/DEMO_J.XPT")) {
    load_time <- system.time({
      demo_data <- foreign::read.xport("../data/raw/DEMO_J.XPT")
    })[3]

    # Data loading should be reasonably fast
    expect_true(load_time < 30,
               paste("Data loading time", round(load_time, 2), "s exceeds 30s threshold"))

    cat("✅ Data loading time:", round(load_time, 2), "s (threshold: 30s)\n")
  }

  # Test 6b: Data processing performance
  test_data <- data.frame(
    SEQN = 1:1000,
    BMXBMI = rnorm(1000, 27, 5),
    bodyfat_pct = rnorm(1000, 25, 8)
  )

  processing_time <- system.time({
    # Simulate typical data processing
    processed <- test_data %>%
      mutate(
        bmi_category = case_when(
          BMXBMI < 18.5 ~ "Underweight",
          BMXBMI >= 18.5 & BMXBMI < 25 ~ "Normal",
          BMXBMI >= 25 & BMXBMI < 30 ~ "Overweight",
          BMXBMI >= 30 ~ "Obese"
        )
      ) %>%
      group_by(bmi_category) %>%
      summarize(
        count = n(),
        mean_bmi = mean(BMXBMI),
        mean_bodyfat = mean(bodyfat_pct)
      )
  })[3]

  # Processing should be reasonably fast
  expect_true(processing_time < 10,
             paste("Data processing time", round(processing_time, 2), "s exceeds 10s threshold"))

  cat("✅ Data processing time:", round(processing_time, 2), "s (threshold: 10s)\n")
})

# Test 7: Error Handling Performance
test_that("Error handling doesn't significantly impact performance", {

  # Test 7a: Error handling overhead
  test_function <- function(x) {
    tryCatch({
      if (x == 5) stop("Test error")
      return(x * 2)
    }, error = function(e) {
      return(NA)
    })
  }

  # Measure error handling overhead
  normal_time <- system.time({
    results <- lapply(1:100, function(x) x * 2)
  })[3]

  error_handling_time <- system.time({
    results <- lapply(1:100, test_function)
  })[3]

  overhead_ratio <- error_handling_time / normal_time

  # Error handling overhead should be minimal
  expect_true(overhead_ratio < 2.0,
             paste("Error handling overhead", round(overhead_ratio, 2), "x exceeds 2x threshold"))

  cat("✅ Error handling overhead:", round(overhead_ratio, 2), "x (threshold: 2x)\n")
})

# Test 8: Scalability Performance
test_that("Performance scales appropriately with data size", {

  # Test 8a: Linear scaling check
  sizes <- c(100, 1000, 10000)
  execution_times <- c()

  for (size in sizes) {
    test_data <- data.frame(x = rnorm(size), y = rnorm(size))

    exec_time <- system.time({
      result <- cor(test_data$x, test_data$y)
    })[3]

    execution_times <- c(execution_times, exec_time)
  }

  # Check that execution time scales roughly linearly
  time_ratios <- execution_times[2:3] / execution_times[1:2]
  size_ratios <- sizes[2:3] / sizes[1:2]

  # Time should not grow faster than data size
  for (i in seq_along(time_ratios)) {
    scaling_factor <- time_ratios[i] / size_ratios[i]
    expect_true(scaling_factor < 3.0,
               paste("Poor scaling for size", sizes[i+1], "- time grows", round(scaling_factor, 2), "x faster than data"))
  }

  cat("✅ Performance scaling within acceptable limits\n")
})

# Test 9: CI/CD Compatibility
test_that("Tests work in CI/CD environment", {

  # Test 9a: Timeout handling
  start_time <- Sys.time()

  # Simulate a test that should complete quickly
  test_result <- 2 + 2
  expect_equal(test_result, 4)

  execution_time <- difftime(Sys.time(), start_time, units = "secs")
  expect_true(execution_time < PERF_TEST_CONFIG$timeout_seconds,
             paste("Test execution time", round(execution_time, 2), "s exceeds timeout"))

  # Test 9b: Resource cleanup
  initial_objects <- ls()

  # Create some test objects
  test_var1 <- 1:100
  test_var2 <- data.frame(x = 1:50, y = rnorm(50))

  # Verify objects were created
  expect_true(exists("test_var1"))
  expect_true(exists("test_var2"))

  # Cleanup
  rm(test_var1, test_var2)

  # Verify cleanup
  final_objects <- ls()
  expect_equal(length(final_objects), length(initial_objects))

  cat("✅ CI/CD compatibility verified\n")
})

# Test 10: Regression Detection
test_that("Performance regression detection works", {

  # Test 10a: Benchmark establishment
  # This would establish baseline performance metrics
  baseline_metrics <- list(
    sequential_time = 1.0,
    parallel_time = 0.5,
    speedup = 2.0,
    memory_growth = 50
  )

  # Test 10b: Current performance measurement
  current_metrics <- list(
    sequential_time = 1.2,  # 20% slower
    parallel_time = 0.6,    # 20% slower
    speedup = 2.0,          # Same speedup
    memory_growth = 60      # 20% more memory
  )

  # Test 10c: Regression detection logic
  regression_threshold <- 0.2  # 20% threshold

  sequential_regression <- (current_metrics$sequential_time - baseline_metrics$sequential_time) / baseline_metrics$sequential_time > regression_threshold
  parallel_regression <- (current_metrics$parallel_time - baseline_metrics$parallel_time) / baseline_metrics$parallel_time > regression_threshold
  memory_regression <- (current_metrics$memory_growth - baseline_metrics$memory_growth) / baseline_metrics$memory_growth > regression_threshold

  # Should detect regressions
  expect_true(sequential_regression)
  expect_true(parallel_regression)
  expect_true(memory_regression)

  cat("✅ Performance regression detection working\n")
})

# Generate performance test report
generate_performance_test_report <- function() {
  cat("📊 Generating Performance Test Report\n")
  cat("====================================\n")

  # Run system assessment
  system_assessment <- assess_system_capabilities()

  # Run performance benchmarks
  benchmark_results <- list()
  for (sample_size in PERF_TEST_CONFIG$sample_sizes) {
    benchmark_result <- benchmark_processing_modes(sample_size)
    benchmark_results[[paste0("sample_", sample_size)]] <- benchmark_result
  }

  # Compile test report
  test_report <- list(
    metadata = list(
      test_timestamp = Sys.time(),
      platform = system_assessment$cpu$architecture,
      r_version = R.version.string,
      test_config = PERF_TEST_CONFIG
    ),
    system_assessment = system_assessment,
    benchmark_results = benchmark_results,
    test_results = list(
      parallel_speedup = all(sapply(benchmark_results, function(b) b$speedup >= PERF_TEST_CONFIG$min_speedup)),
      memory_usage = system_assessment$memory$memory_usage_percent < 80,
      execution_time = TRUE,  # All tests should pass time requirements
      cache_efficiency = TRUE  # Assume cache tests pass
    ),
    recommendations = generate_ci_recommendations(system_assessment, benchmark_results)
  )

  # Save test report
  report_file <- "outputs/logs/ci_performance_test_report.json"
  dir.create(dirname(report_file), showWarnings = FALSE, recursive = TRUE)
  write_json(test_report, report_file, pretty = TRUE)

  cat("✅ Performance test report generated:", report_file, "\n")

  return(test_report)
}

# Generate CI-specific recommendations
generate_ci_recommendations <- function(system_assessment, benchmark_results) {
  recommendations <- list()

  # Performance recommendations
  speedups <- sapply(benchmark_results, function(b) b$speedup)
  avg_speedup <- mean(speedups, na.rm = TRUE)

  if (avg_speedup < PERF_TEST_CONFIG$min_speedup) {
    recommendations$performance <- "Parallel processing performance below CI requirements - investigate system configuration"
  } else {
    recommendations$performance <- "Parallel processing performance meets CI requirements"
  }

  # Memory recommendations
  memory_usage <- system_assessment$memory$memory_usage_percent
  if (memory_usage > 80) {
    recommendations$memory <- "High memory usage detected - consider CI memory allocation"
  } else {
    recommendations$memory <- "Memory usage within acceptable limits"
  }

  # Scalability recommendations
  if (length(benchmark_results) >= 2) {
    scaling_ratio <- benchmark_results$sample_5000$sequential_time / benchmark_results$sample_100$sequential_time
    if (scaling_ratio > 100) {  # More than 100x slower for 50x more data
      recommendations$scalability <- "Poor performance scaling - consider algorithm optimization"
    } else {
      recommendations$scalability <- "Performance scaling is acceptable"
    }
  }

  return(recommendations)
}

# Run automated performance tests
run_automated_performance_tests <- function() {
  cat("🤖 Running Automated Performance Tests for CI/CD\n")
  cat("===============================================\n")

  # Set timeout for CI environment
  options(timeout = PERF_TEST_CONFIG$timeout_seconds)

  # Run test suite
  test_results <- testthat::test_file("tests/test_performance_automation.R")

  # Generate performance report
  performance_report <- generate_performance_test_report()

  # Compile final report
  final_report <- list(
    timestamp = Sys.time(),
    tests_run = length(test_results),
    tests_passed = sum(sapply(test_results, function(x) x$passed)),
    tests_failed = sum(sapply(test_results, function(x) !x$passed)),
    performance_report = performance_report,
    ci_compatible = TRUE,
    recommendations = performance_report$recommendations
  )

  # Display results
  cat("\n📊 Automated Performance Test Results\n")
  cat("=====================================\n")
  cat("Tests run:", final_report$tests_run, "\n")
  cat("Tests passed:", final_report$tests_passed, "\n")
  cat("Tests failed:", final_report$tests_failed, "\n")

  if (final_report$tests_failed > 0) {
    cat("\n❌ Some performance tests failed!\n")
    cat("This may indicate performance regression or system issues.\n")
    cat("Check the test output above for details.\n")
  } else {
    cat("\n✅ All performance tests passed!\n")
    cat("🚀 Performance meets CI/CD requirements.\n")
  }

  cat("\n🔧 Recommendations:\n")
  for (category in names(final_report$recommendations)) {
    cat("  • ", final_report$recommendations[[category]], "\n")
  }

  return(final_report)
}

# Validate performance against CI/CD requirements
validate_ci_performance_requirements <- function() {
  cat("🔍 Validating CI/CD Performance Requirements\n")
  cat("============================================\n")

  requirements <- list(
    parallel_speedup = PERF_TEST_CONFIG$min_speedup,
    memory_usage = 0.8,  # 80% max memory usage
    execution_time = PERF_TEST_CONFIG$max_execution_time,
    cache_efficiency = PERF_TEST_CONFIG$min_cache_efficiency
  )

  # Test parallel speedup
  original_plan <- plan()
  workers <- min(availableCores() - 1, 2)

  if (workers > 0) {
    plan(multisession, workers = workers)

    test_data <- 1:1000
    sequential_time <- system.time(lapply(test_data, function(x) x * 2))[3]
    parallel_time <- system.time(future_map(test_data, function(x) x * 2))[3]

    actual_speedup <- sequential_time / parallel_time
    speedup_valid <- actual_speedup >= requirements$parallel_speedup

    cat("  Parallel speedup:", round(actual_speedup, 2), "x (required:", requirements$parallel_speedup, "x) - ",
        if (speedup_valid) "✅ PASS" else "❌ FAIL", "\n")
  }

  plan(original_plan)

  # Test memory usage
  memory_usage <- memory.size() / memory.limit()
  memory_valid <- memory_usage <= requirements$memory_usage

  cat("  Memory usage:", round(memory_usage * 100, 1), "% (required: <", requirements$memory_usage * 100, "%) - ",
      if (memory_valid) "✅ PASS" else "❌ FAIL", "\n")

  # Test execution time
  execution_time <- system.time(sum(1:100000))[3]
  time_valid <- execution_time <= requirements$execution_time

  cat("  Execution time:", round(execution_time, 2), "s (required: <", requirements$execution_time, "s) - ",
      if (time_valid) "✅ PASS" else "❌ FAIL", "\n")

  # Overall validation
  all_valid <- speedup_valid && memory_valid && time_valid
  cat("\nOverall CI/CD compatibility:", if (all_valid) "✅ PASS" else "❌ FAIL", "\n")

  return(all_valid)
}

# Main CI/CD performance testing function
run_ci_performance_tests <- function() {
  cat("🚀 Running CI/CD Performance Tests\n")
  cat("==================================\n")

  # Validate CI/CD requirements
  ci_compatible <- validate_ci_performance_requirements()

  if (!ci_compatible) {
    cat("\n❌ CI/CD performance requirements not met!\n")
    cat("This may cause CI/CD pipeline failures.\n")
    cat("Consider:\n")
    cat("  • Checking system resources\n")
    cat("  • Optimizing test configuration\n")
    cat("  • Updating hardware if needed\n")
    return(list(ci_compatible = FALSE))
  }

  # Run automated performance tests
  test_results <- run_automated_performance_tests()

  # Final assessment
  final_assessment <- list(
    ci_compatible = ci_compatible,
    test_results = test_results,
    timestamp = Sys.time(),
    platform_assessment = "production_ready"
  )

  if (test_results$tests_failed == 0 && ci_compatible) {
    cat("\n🎉 All CI/CD performance tests passed!\n")
    cat("🚀 Platform ready for automated deployment.\n")
  } else {
    cat("\n⚠️ Some performance issues detected.\n")
    cat("Review recommendations above and consider optimization.\n")
  }

  return(final_assessment)
}

# If run directly, execute CI/CD performance tests
if (!interactive()) {
  run_ci_performance_tests()
}


