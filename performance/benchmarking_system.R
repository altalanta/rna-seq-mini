#!/usr/bin/env Rscript
# Comprehensive Performance Benchmarking System for NHANES BMI Body Fat Analysis Platform

library(microbenchmark)
library(profvis)
library(bench)
library(dplyr)
library(ggplot2)
library(future)
library(furrr)
library(digest)
library(yaml)
library(jsonlite)

# Source required modules
source("../R/data_versioning.R")
source("../parallel_pipeline.R")

# Benchmarking configuration
BENCHMARK_CONFIG <- list(
  iterations = 10,
  warmup_iterations = 3,
  sample_sizes = c(100, 500, 1000, 2500, 5000),
  parallel_configs = c(1, 2, 4, 8),
  memory_thresholds = list(
    warning = 512 * 1024 * 1024,  # 512MB
    critical = 1024 * 1024 * 1024  # 1GB
  ),
  time_thresholds = list(
    excellent = 10,  # 10 seconds
    good = 30,       # 30 seconds
    acceptable = 60  # 60 seconds
  )
)

# System performance assessment
assess_system_capabilities <- function() {
  cat("🔍 Assessing system performance capabilities...\n")

  # CPU information
  cpu_cores <- availableCores()
  cpu_info <- list(
    logical_cores = cpu_cores,
    physical_cores = cpu_cores %/% 2,  # Rough estimate
    architecture = Sys.info()["machine"],
    cpu_usage = tryCatch({
      # Get CPU usage (platform-specific)
      if (.Platform$OS.type == "unix") {
        system("top -bn1 | grep 'Cpu(s)' | sed 's/.*, *\\([0-9.]*\\)%* id.*/\\1/' | awk '{print 100 - $1}'", intern = TRUE)
      } else {
        NA
      }
    }, error = function(e) NA)
  )

  # Memory information
  memory_info <- list(
    total_memory = memory.limit(),
    available_memory = memory.size(),
    memory_usage_percent = (memory.size() / memory.limit()) * 100,
    memory_pressure = tryCatch({
      # Check memory pressure
      gc_info <- gc()
      (gc_info[2, 6] / memory.limit()) * 100  # GC overhead
    }, error = function(e) NA)
  )

  # Storage information
  storage_info <- list(
    working_directory = getwd(),
    available_space = tryCatch({
      info <- file.info(getwd())
      info$size / 1024 / 1024 / 1024  # GB
    }, error = function(e) NA),
    disk_io_performance = tryCatch({
      # Simple I/O benchmark
      test_file <- tempfile()
      test_data <- rnorm(100000)

      start_time <- Sys.time()
      saveRDS(test_data, test_file)
      write_time <- difftime(Sys.time(), start_time, units = "secs")

      start_time <- Sys.time()
      loaded_data <- readRDS(test_file)
      read_time <- difftime(Sys.time(), start_time, units = "secs")

      file.remove(test_file)
      list(write_time = write_time, read_time = read_time)
    }, error = function(e) NA)
  )

  # Network information
  network_info <- list(
    has_internet = tryCatch({
      url("https://www.google.com", open = "rb")
      TRUE
    }, error = function(e) FALSE),
    latency_test = tryCatch({
      start_time <- Sys.time()
      url("https://www.google.com", open = "rb")
      difftime(Sys.time(), start_time, units = "secs")
    }, error = function(e) NA)
  )

  system_assessment <- list(
    timestamp = Sys.time(),
    cpu = cpu_info,
    memory = memory_info,
    storage = storage_info,
    network = network_info,
    performance_score = calculate_system_score(cpu_info, memory_info, storage_info, network_info)
  )

  cat("✅ System assessment complete\n")
  return(system_assessment)
}

# Calculate overall system performance score
calculate_system_score <- function(cpu_info, memory_info, storage_info, network_info) {
  score <- 0

  # CPU score (0-30 points)
  cpu_score <- min(cpu_info$logical_cores * 3, 30)
  if (!is.na(cpu_info$cpu_usage) && cpu_info$cpu_usage < 50) {
    cpu_score <- cpu_score + 5  # Bonus for low CPU usage
  }
  score <- score + cpu_score

  # Memory score (0-30 points)
  memory_score <- 30
  if (memory_info$memory_usage_percent > 80) {
    memory_score <- memory_score - 20
  } else if (memory_info$memory_usage_percent > 60) {
    memory_score <- memory_score - 10
  }
  if (!is.na(memory_info$memory_pressure) && memory_info$memory_pressure > 10) {
    memory_score <- memory_score - 5  # Penalty for high GC overhead
  }
  score <- score + max(memory_score, 0)

  # Storage score (0-20 points)
  storage_score <- 20
  if (!is.na(storage_info$available_space) && storage_info$available_space < 10) {
    storage_score <- storage_score - 10
  }
  if (!is.null(storage_info$disk_io_performance)) {
    avg_io_time <- mean(c(storage_info$disk_io_performance$write_time,
                         storage_info$disk_io_performance$read_time), na.rm = TRUE)
    if (avg_io_time > 2) {
      storage_score <- storage_score - 5  # Penalty for slow I/O
    }
  }
  score <- score + max(storage_score, 0)

  # Network score (0-20 points)
  network_score <- 20
  if (!network_info$has_internet) {
    network_score <- network_score - 10
  }
  if (!is.na(network_info$latency_test) && network_info$latency_test > 1) {
    network_score <- network_score - 5  # Penalty for high latency
  }
  score <- score + max(network_score, 0)

  return(min(score, 100))  # Cap at 100
}

# Benchmark sequential vs parallel processing
benchmark_processing_modes <- function(sample_size = 1000) {
  cat("⚡ Benchmarking processing modes with sample size:", sample_size, "\n")

  # Create test dataset
  set.seed(123)
  test_data <- data.frame(
    SEQN = 1:sample_size,
    BMXBMI = rnorm(sample_size, 27, 5),
    bodyfat_pct = rnorm(sample_size, 25, 8),
    sex = sample(c("Male", "Female"), sample_size, replace = TRUE)
  )

  # Sequential processing benchmark
  sequential_benchmark <- microbenchmark({
    # Simulate sequential correlation analysis
    results <- list()
    for (group in c("Overall", "Male", "Female")) {
      if (group == "Overall") {
        subset_data <- test_data
      } else {
        subset_data <- test_data[test_data$sex == group, ]
      }

      if (nrow(subset_data) > 0) {
        corr <- cor(subset_data$BMXBMI, subset_data$bodyfat_pct)
        results[[group]] <- corr
      }
    }
  }, times = BENCHMARK_CONFIG$iterations)

  # Parallel processing benchmark
  original_plan <- plan()
  workers <- min(availableCores() - 1, 4)  # Limit for benchmarking

  if (workers > 0) {
    plan(multisession, workers = workers)

    parallel_benchmark <- microbenchmark({
      results <- future_map(c("Overall", "Male", "Female"), function(group) {
        if (group == "Overall") {
          subset_data <- test_data
        } else {
          subset_data <- test_data[test_data$sex == group, ]
        }

        if (nrow(subset_data) > 0) {
          return(cor(subset_data$BMXBMI, subset_data$bodyfat_pct))
        } else {
          return(NA)
        }
      })
    }, times = BENCHMARK_CONFIG$iterations)

    plan(original_plan)

    # Calculate performance metrics
    sequential_mean <- mean(sequential_benchmark$time) / 1e9
    parallel_mean <- mean(parallel_benchmark$time) / 1e9

    speedup <- sequential_mean / parallel_mean
    efficiency <- speedup / workers * 100

    benchmark_results <- list(
      sample_size = sample_size,
      sequential_time = sequential_mean,
      parallel_time = parallel_mean,
      speedup = speedup,
      efficiency = efficiency,
      workers_used = workers,
      iterations = BENCHMARK_CONFIG$iterations,
      timestamp = Sys.time()
    )

    cat("  Sequential time:", round(sequential_mean, 3), "seconds\n")
    cat("  Parallel time:", round(parallel_mean, 3), "seconds\n")
    cat("  Speedup:", round(speedup, 2), "x\n")
    cat("  Efficiency:", round(efficiency, 1), "%\n")

  } else {
    # Fallback to sequential if parallel not available
    benchmark_results <- list(
      sample_size = sample_size,
      sequential_time = mean(sequential_benchmark$time) / 1e9,
      parallel_time = NA,
      speedup = NA,
      efficiency = NA,
      workers_used = 0,
      iterations = BENCHMARK_CONFIG$iterations,
      timestamp = Sys.time(),
      note = "Parallel processing not available"
    )
  }

  return(benchmark_results)
}

# Benchmark cache performance
benchmark_cache_performance <- function() {
  cat("💾 Benchmarking cache performance...\n")

  cache_dir <- "cache"
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Create test datasets of different sizes
  test_datasets <- list()
  for (size in c(100, 1000, 10000)) {
    test_data <- data.frame(
      x = rnorm(size),
      y = rnorm(size),
      z = rnorm(size)
    )
    test_datasets[[paste0("size_", size)]] <- test_data
  }

  cache_results <- list()

  for (dataset_name in names(test_datasets)) {
    test_data <- test_datasets[[dataset_name]]
    cache_key <- digest(test_data)

    # Cache miss (first save)
    start_time <- Sys.time()
    save_to_cache(cache_key, test_data)
    first_save_time <- difftime(Sys.time(), start_time, units = "secs")

    # Cache hit (load)
    start_time <- Sys.time()
    loaded_data <- load_from_cache(cache_key)
    cache_hit_time <- difftime(Sys.time(), start_time, units = "secs")

    # Verify data integrity
    data_integrity <- all.equal(test_data, loaded_data)

    cache_results[[dataset_name]] <- list(
      dataset_size = nrow(test_data),
      first_save_time = first_save_time,
      cache_hit_time = cache_hit_time,
      speedup_ratio = first_save_time / cache_hit_time,
      data_integrity = data_integrity
    )

    cat("  Dataset size:", nrow(test_data), "rows\n")
    cat("    First save:", round(first_save_time, 4), "seconds\n")
    cat("    Cache hit:", round(cache_hit_time, 4), "seconds\n")
    cat("    Speedup:", round(first_save_time / cache_hit_time, 1), "x\n")
  }

  return(cache_results)
}

# Benchmark memory usage patterns
benchmark_memory_usage <- function() {
  cat("🧠 Benchmarking memory usage patterns...\n")

  memory_results <- list()

  # Test different dataset sizes
  for (size in c(1000, 10000, 100000)) {
    cat("  Testing dataset size:", size, "rows\n")

    # Initial memory state
    initial_memory <- memory.size()
    initial_gc <- gc()

    # Create dataset
    start_time <- Sys.time()
    test_data <- data.frame(
      matrix(rnorm(size * 10), ncol = 10)
    )
    creation_time <- difftime(Sys.time(), start_time, units = "secs")

    # Memory after creation
    after_creation <- memory.size()
    creation_growth <- after_creation - initial_memory

    # Process data
    start_time <- Sys.time()
    processed_data <- test_data %>%
      mutate(
        row_sum = rowSums(.),
        row_mean = rowMeans(.),
        category = cut(row_sum, breaks = 5, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))
      ) %>%
      group_by(category) %>%
      summarize(
        count = n(),
        avg_sum = mean(row_sum),
        sd_sum = sd(row_sum)
      )
    processing_time <- difftime(Sys.time(), start_time, units = "secs")

    # Memory after processing
    after_processing <- memory.size()
    processing_growth <- after_processing - after_creation

    # Final cleanup
    rm(test_data, processed_data)
    final_gc <- gc()
    final_memory <- memory.size()

    memory_results[[paste0("size_", size)]] <- list(
      dataset_size = size,
      creation_time = creation_time,
      processing_time = processing_time,
      initial_memory = initial_memory,
      creation_growth = creation_growth,
      processing_growth = processing_growth,
      final_memory = final_memory,
      total_growth = final_memory - initial_memory,
      gc_efficiency = (final_gc[2, 6] - initial_gc[2, 6]) / (final_gc[1, 6] - initial_gc[1, 6]) * 100
    )
  }

  # Display results
  for (result_name in names(memory_results)) {
    result <- memory_results[[result_name]]
    cat("    Dataset size:", result$dataset_size, "rows\n")
    cat("      Creation time:", round(result$creation_time, 3), "seconds\n")
    cat("      Processing time:", round(result$processing_time, 3), "seconds\n")
    cat("      Memory growth:", round(result$total_growth, 1), "MB\n")
    cat("      GC efficiency:", round(result$gc_efficiency, 1), "%\n")
  }

  return(memory_results)
}

# Comprehensive performance benchmarking
run_comprehensive_benchmark <- function() {
  cat("📊 Running comprehensive performance benchmark...\n")
  cat("==============================================\n")

  start_time <- Sys.time()

  # 1. System assessment
  cat("Step 1: System capability assessment\n")
  system_assessment <- assess_system_capabilities()

  # 2. Processing mode benchmarks
  cat("\nStep 2: Processing mode benchmarks\n")
  processing_benchmarks <- list()
  for (sample_size in BENCHMARK_CONFIG$sample_sizes) {
    benchmark_result <- benchmark_processing_modes(sample_size)
    processing_benchmarks[[paste0("sample_", sample_size)]] <- benchmark_result
  }

  # 3. Cache performance
  cat("\nStep 3: Cache performance analysis\n")
  cache_benchmarks <- benchmark_cache_performance()

  # 4. Memory usage patterns
  cat("\nStep 4: Memory usage analysis\n")
  memory_benchmarks <- benchmark_memory_usage()

  # 5. Compile comprehensive report
  end_time <- Sys.time()
  total_benchmark_time <- difftime(end_time, start_time, units = "secs")

  benchmark_report <- list(
    metadata = list(
      benchmark_timestamp = Sys.time(),
      total_benchmark_time = total_benchmark_time,
      platform = system_assessment$cpu$architecture,
      r_version = R.version.string
    ),
    system_assessment = system_assessment,
    processing_benchmarks = processing_benchmarks,
    cache_benchmarks = cache_benchmarks,
    memory_benchmarks = memory_benchmarks,
    recommendations = generate_benchmark_recommendations(
      system_assessment, processing_benchmarks, cache_benchmarks, memory_benchmarks
    )
  )

  # Save benchmark report
  report_file <- "outputs/logs/comprehensive_benchmark_report.json"
  dir.create(dirname(report_file), showWarnings = FALSE, recursive = TRUE)
  write_json(benchmark_report, report_file, pretty = TRUE)

  cat("✅ Comprehensive benchmark completed in", round(total_benchmark_time, 2), "seconds\n")
  cat("📊 Report saved to:", report_file, "\n")

  return(benchmark_report)
}

# Generate benchmark recommendations
generate_benchmark_recommendations <- function(system_assessment, processing_benchmarks, cache_benchmarks, memory_benchmarks) {
  recommendations <- list()

  # System recommendations
  system_score <- system_assessment$performance_score
  if (system_score < 60) {
    recommendations$system <- "Consider system upgrades for optimal performance"
  } else {
    recommendations$system <- "System performance is adequate"
  }

  # Parallel processing recommendations
  if (length(processing_benchmarks) > 0) {
    # Find best performing configuration
    speedups <- sapply(processing_benchmarks, function(b) b$speedup)
    best_config <- names(processing_benchmarks)[which.max(speedups)]

    if (max(speedups, na.rm = TRUE) < 1.5) {
      recommendations$parallel <- "Parallel processing may not provide significant benefits on this system"
    } else {
      recommendations$parallel <- paste("Use", best_config, "configuration for optimal parallel performance")
    }
  }

  # Cache recommendations
  if (length(cache_benchmarks) > 0) {
    cache_speedups <- sapply(cache_benchmarks, function(b) b$speedup_ratio)
    avg_cache_speedup <- mean(cache_speedups, na.rm = TRUE)

    if (avg_cache_speedup < 10) {
      recommendations$cache <- "Cache performance could be improved - consider SSD storage"
    } else {
      recommendations$cache <- "Cache performance is excellent"
    }
  }

  # Memory recommendations
  if (length(memory_benchmarks) > 0) {
    total_growth <- sum(sapply(memory_benchmarks, function(b) b$total_growth))

    if (total_growth > 500) {
      recommendations$memory <- "High memory usage detected - consider memory optimization"
    } else {
      recommendations$memory <- "Memory usage is within acceptable limits"
    }
  }

  return(recommendations)
}

# Performance regression detection
detect_performance_regression <- function(current_benchmark, historical_benchmarks = NULL) {
  cat("🔍 Detecting performance regression...\n")

  regression_results <- list()

  # Load historical benchmarks if not provided
  if (is.null(historical_benchmarks)) {
    historical_file <- "outputs/logs/historical_benchmarks.json"
    if (file.exists(historical_file)) {
      historical_benchmarks <- fromJSON(historical_file)
    } else {
      cat("📭 No historical benchmarks found - establishing baseline\n")
      regression_results$baseline <- TRUE
      regression_results$recommendations <- "This is the first benchmark - use as baseline for future comparisons"
      return(regression_results)
    }
  }

  # Compare current vs historical performance
  current_processing <- current_benchmark$processing_benchmarks
  historical_processing <- historical_benchmarks$processing_benchmarks

  regression_detected <- FALSE
  regression_details <- list()

  # Compare processing performance
  for (config_name in names(current_processing)) {
    if (config_name %in% names(historical_processing)) {
      current_time <- current_processing[[config_name]]$parallel_time
      historical_time <- historical_processing[[config_name]]$parallel_time

      if (!is.na(current_time) && !is.na(historical_time)) {
        performance_change <- (current_time - historical_time) / historical_time * 100

        if (performance_change > 20) {  # 20% slowdown
          regression_detected <- TRUE
          regression_details[[config_name]] <- list(
            change_percent = performance_change,
            severity = "significant",
            recommendation = "Investigate performance degradation"
          )
        } else if (performance_change > 10) {  # 10% slowdown
          regression_details[[config_name]] <- list(
            change_percent = performance_change,
            severity = "moderate",
            recommendation = "Monitor performance trend"
          )
        }
      }
    }
  }

  regression_results$regression_detected <- regression_detected
  regression_results$details <- regression_details
  regression_results$timestamp <- Sys.time()

  if (regression_detected) {
    cat("⚠️ Performance regression detected\n")
    for (config in names(regression_details)) {
      detail <- regression_details[[config]]
      cat("  •", config, ": ", round(detail$change_percent, 1), "% slower\n")
    }
  } else {
    cat("✅ No performance regression detected\n")
  }

  return(regression_results)
}

# Generate performance optimization report
generate_optimization_report <- function() {
  cat("📈 Generating performance optimization report...\n")

  # Run comprehensive benchmark
  current_benchmark <- run_comprehensive_benchmark()

  # Detect regression
  regression_analysis <- detect_performance_regression(current_benchmark)

  # Compile optimization report
  optimization_report <- list(
    metadata = list(
      report_timestamp = Sys.time(),
      platform_assessment = current_benchmark$system_assessment$performance_score,
      benchmark_version = "1.0"
    ),
    current_performance = current_benchmark,
    regression_analysis = regression_analysis,
    optimization_strategies = generate_optimization_strategies(current_benchmark),
    hardware_recommendations = generate_hardware_recommendations(current_benchmark),
    software_optimizations = generate_software_optimizations(current_benchmark)
  )

  # Save optimization report
  report_file <- "outputs/logs/optimization_report.json"
  dir.create(dirname(report_file), showWarnings = FALSE, recursive = TRUE)
  write_json(optimization_report, report_file, pretty = TRUE)

  cat("✅ Optimization report generated:", report_file, "\n")

  return(optimization_report)
}

# Generate optimization strategies
generate_optimization_strategies <- function(benchmark) {
  strategies <- list()

  # Parallel processing optimization
  processing_benchmarks <- benchmark$processing_benchmarks
  if (length(processing_benchmarks) > 0) {
    # Find optimal configuration
    speedups <- sapply(processing_benchmarks, function(b) b$speedup)
    best_config <- names(processing_benchmarks)[which.max(speedups)]

    strategies$parallel_processing <- list(
      optimal_configuration = best_config,
      expected_speedup = max(speedups, na.rm = TRUE),
      recommendation = paste("Use", best_config, "for optimal parallel performance")
    )
  }

  # Cache optimization
  cache_benchmarks <- benchmark$cache_benchmarks
  if (length(cache_benchmarks) > 0) {
    cache_speedups <- sapply(cache_benchmarks, function(b) b$speedup_ratio)
    avg_cache_speedup <- mean(cache_speedups, na.rm = TRUE)

    strategies$cache_optimization <- list(
      average_speedup = avg_cache_speedup,
      recommendation = if (avg_cache_speedup > 50) {
        "Cache performance is excellent - maintain current configuration"
      } else {
        "Consider SSD storage for improved cache performance"
      }
    )
  }

  # Memory optimization
  memory_benchmarks <- benchmark$memory_benchmarks
  if (length(memory_benchmarks) > 0) {
    total_growth <- sum(sapply(memory_benchmarks, function(b) b$total_growth))

    strategies$memory_optimization <- list(
      total_memory_growth = total_growth,
      recommendation = if (total_growth > 500) {
        "Implement memory optimization strategies - process data in chunks"
      } else {
        "Memory usage is within acceptable limits"
      }
    )
  }

  return(strategies)
}

# Generate hardware recommendations
generate_hardware_recommendations <- function(benchmark) {
  system_assessment <- benchmark$system_assessment
  recommendations <- list()

  # CPU recommendations
  if (system_assessment$cpu$logical_cores < 4) {
    recommendations$cpu <- list(
      current = system_assessment$cpu$logical_cores,
      recommended = 8,
      rationale = "More cores enable better parallel processing performance",
      impact = "High - significant speedup for intensive analyses"
    )
  }

  # Memory recommendations
  memory_usage <- system_assessment$memory$memory_usage_percent
  if (memory_usage > 70) {
    recommendations$memory <- list(
      current = round(system_assessment$memory$available_memory, 1),
      recommended = 16,
      rationale = "More memory reduces swapping and improves performance",
      impact = "Medium - prevents memory-related slowdowns"
    )
  }

  # Storage recommendations
  if (!is.na(system_assessment$storage$available_space) && system_assessment$storage$available_space < 50) {
    recommendations$storage <- list(
      current = round(system_assessment$storage$available_space, 1),
      recommended = 100,
      rationale = "More storage space prevents disk space issues",
      impact = "Low - prevents analysis interruption"
    )
  }

  return(recommendations)
}

# Generate software optimizations
generate_software_optimizations <- function(benchmark) {
  optimizations <- list()

  # R package optimization
  optimizations$r_optimization <- list(
    recommendations = c(
      "Keep R packages updated for performance improvements",
      "Use renv for reproducible package environments",
      "Consider compiled packages for intensive computations"
    )
  )

  # Code optimization
  optimizations$code_optimization <- list(
    recommendations = c(
      "Use vectorized operations instead of loops",
      "Pre-allocate memory for known data structures",
      "Use data.table for large dataset operations",
      "Consider parallel processing for independent computations"
    )
  )

  # Cache optimization
  cache_benchmarks <- benchmark$cache_benchmarks
  if (length(cache_benchmarks) > 0) {
    cache_speedups <- sapply(cache_benchmarks, function(b) b$speedup_ratio)

    optimizations$cache_optimization <- list(
      current_performance = mean(cache_speedups, na.rm = TRUE),
      recommendations = if (mean(cache_speedups, na.rm = TRUE) < 20) {
        c("Consider SSD storage for cache directory", "Optimize cache key generation")
      } else {
        c("Cache performance is optimal")
      }
    )
  }

  return(optimizations)
}

# Display benchmark summary
display_benchmark_summary <- function() {
  cat("📊 Performance Benchmark Summary\n")
  cat("===============================\n")

  # System assessment
  system_assessment <- assess_system_capabilities()
  cat("💻 System Performance Score:", system_assessment$performance_score, "/100\n")

  # Processing benchmarks
  cat("\n⚡ Processing Performance:\n")
  for (sample_size in BENCHMARK_CONFIG$sample_sizes) {
    tryCatch({
      benchmark_result <- benchmark_processing_modes(sample_size)
      if (!is.na(benchmark_result$speedup)) {
        cat("  Sample size", sample_size, ": ", round(benchmark_result$speedup, 2), "x speedup\n")
      }
    }, error = function(e) {
      cat("  Sample size", sample_size, ": Error in benchmarking\n")
    })
  }

  # Cache performance
  cat("\n💾 Cache Performance:\n")
  tryCatch({
    cache_benchmarks <- benchmark_cache_performance()
    cache_speedups <- sapply(cache_benchmarks, function(b) b$speedup_ratio)
    cat("  Average cache speedup:", round(mean(cache_speedups, na.rm = TRUE), 1), "x\n")
  }, error = function(e) {
    cat("  Cache performance: Error in analysis\n")
  })

  # Memory performance
  cat("\n🧠 Memory Performance:\n")
  tryCatch({
    memory_benchmarks <- benchmark_memory_usage()
    total_growth <- sum(sapply(memory_benchmarks, function(b) b$total_growth))
    cat("  Total memory growth:", round(total_growth, 1), "MB\n")
  }, error = function(e) {
    cat("  Memory performance: Error in analysis\n")
  })

  cat("\n🔧 Recommendations:\n")
  tryCatch({
    recommendations <- generate_benchmark_recommendations(
      system_assessment,
      list(),  # Simplified for display
      list(),
      list()
    )

    for (category in names(recommendations)) {
      cat("  • ", recommendations[[category]], "\n")
    }
  }, error = function(e) {
    cat("  • Run full benchmark for detailed recommendations\n")
  })

  cat("\n📋 Available commands:\n")
  cat("  Rscript performance_tools.R benchmark    # Run benchmarks\n")
  cat("  Rscript performance_tools.R optimize     # Auto-optimize settings\n")
  cat("  Rscript performance_tools.R report       # Generate full report\n")
  cat("  Rscript performance_tools.R regression   # Check for regressions\n")
}

# Main benchmarking function
run_performance_benchmark <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  action <- if (length(args) > 0) args[1] else "summary"

  switch(action,
    "full" = {
      cat("🔬 Running comprehensive performance benchmark...\n")
      benchmark_report <- run_comprehensive_benchmark()
      return(benchmark_report)
    },
    "regression" = {
      cat("🔍 Checking for performance regression...\n")
      current_benchmark <- run_comprehensive_benchmark()
      regression_results <- detect_performance_regression(current_benchmark)
      return(regression_results)
    },
    "optimize" = {
      cat("🔧 Generating optimization recommendations...\n")
      optimization_report <- generate_optimization_report()
      return(optimization_report)
    },
    # Default: summary
    {
      display_benchmark_summary()
      return(NULL)
    }
  )
}

# Run benchmarking if script is executed directly
if (!interactive()) {
  run_performance_benchmark()
}


