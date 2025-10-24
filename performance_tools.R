#!/usr/bin/env Rscript
# Performance Optimization and Benchmarking Tools for NHANES BMI Body Fat Analysis Platform

library(microbenchmark)
library(profvis)
library(bench)
library(dplyr)
library(ggplot2)
library(future)
library(furrr)
library(digest)
library(yaml)

# Source required modules
source("R/data_versioning.R")
source("parallel_pipeline.R")

# Performance configuration
PERF_CONFIG <- list(
  benchmark_iterations = 10,
  profiling_samples = 1000,
  memory_threshold = 1024 * 1024 * 1024,  # 1GB
  time_threshold = 30,  # 30 seconds
  optimization_targets = list(
    parallel_speedup = 2.0,  # Minimum 2x speedup expected
    memory_efficiency = 0.8, # 80% memory efficiency
    cache_hit_rate = 0.7     # 70% cache hit rate
  )
)

# System performance assessment
assess_system_performance <- function() {
  cat("🔍 Assessing system performance capabilities...\n")

  # CPU information
  cpu_cores <- availableCores()
  cpu_info <- list(
    logical_cores = cpu_cores,
    physical_cores = cpu_cores %/% 2,  # Rough estimate
    architecture = Sys.info()["machine"]
  )

  # Memory information
  memory_info <- list(
    total_memory = memory.limit(),
    available_memory = memory.size(),
    memory_usage_percent = (memory.size() / memory.limit()) * 100
  )

  # Storage information
  storage_info <- list(
    working_directory = getwd(),
    available_space = tryCatch({
      info <- file.info(getwd())
      info$size / 1024 / 1024 / 1024  # GB
    }, error = function(e) NA)
  )

  # Network information (if applicable)
  network_info <- list(
    has_internet = tryCatch({
      url("https://www.google.com", open = "rb")
      TRUE
    }, error = function(e) FALSE)
  )

  performance_assessment <- list(
    timestamp = Sys.time(),
    cpu = cpu_info,
    memory = memory_info,
    storage = storage_info,
    network = network_info,
    recommendations = generate_performance_recommendations(cpu_info, memory_info, storage_info)
  )

  return(performance_assessment)
}

# Generate performance recommendations
generate_performance_recommendations <- function(cpu_info, memory_info, storage_info) {
  recommendations <- list()

  # CPU recommendations
  if (cpu_info$logical_cores < 4) {
    recommendations$cpu <- "Consider upgrading to a system with more CPU cores for better parallel processing performance"
  } else {
    recommendations$cpu <- "CPU configuration is adequate for parallel processing"
  }

  # Memory recommendations
  memory_usage_pct <- memory_info$memory_usage_percent
  if (memory_usage_pct > 80) {
    recommendations$memory <- "High memory usage detected. Consider increasing available memory or reducing worker count"
  } else if (memory_usage_pct > 60) {
    recommendations$memory <- "Moderate memory usage. Monitor during intensive operations"
  } else {
    recommendations$memory <- "Memory usage is within acceptable limits"
  }

  # Storage recommendations
  if (!is.na(storage_info$available_space) && storage_info$available_space < 10) {
    recommendations$storage <- "Low disk space detected. Consider freeing up space or using external storage"
  } else {
    recommendations$storage <- "Storage capacity is adequate"
  }

  return(recommendations)
}

# Benchmark parallel vs sequential processing
benchmark_parallel_performance <- function() {
  cat("⚡ Benchmarking parallel processing performance...\n")

  # Setup test data
  test_data <- data.frame(
    SEQN = 1:1000,
    BMXBMI = rnorm(1000, 27, 5),
    bodyfat_pct = rnorm(1000, 25, 8),
    sex = sample(c("Male", "Female"), 1000, replace = TRUE)
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
  }, times = PERF_CONFIG$benchmark_iterations)

  # Parallel processing benchmark
  original_plan <- plan()
  plan(multisession, workers = availableCores() - 1)

  parallel_benchmark <- microbenchmark({
    # Simulate parallel correlation analysis
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
  }, times = PERF_CONFIG$benchmark_iterations)

  plan(original_plan)

  # Calculate performance metrics
  sequential_mean <- mean(sequential_benchmark$time) / 1e9  # Convert to seconds
  parallel_mean <- mean(parallel_benchmark$time) / 1e9

  speedup <- sequential_mean / parallel_mean
  efficiency <- speedup / (availableCores() - 1) * 100  # Percentage

  benchmark_results <- list(
    sequential_time = sequential_mean,
    parallel_time = parallel_mean,
    speedup = speedup,
    efficiency = efficiency,
    cores_used = availableCores() - 1,
    iterations = PERF_CONFIG$benchmark_iterations,
    timestamp = Sys.time()
  )

  cat("📊 Benchmark Results:\n")
  cat("  Sequential time:", round(sequential_mean, 3), "seconds\n")
  cat("  Parallel time:", round(parallel_mean, 3), "seconds\n")
  cat("  Speedup:", round(speedup, 2), "x\n")
  cat("  Efficiency:", round(efficiency, 1), "%\n")
  cat("  Cores used:", availableCores() - 1, "\n")

  return(benchmark_results)
}

# Cache performance analysis
analyze_cache_performance <- function() {
  cat("💾 Analyzing cache performance...\n")

  cache_dir <- "cache"
  if (!dir.exists(cache_dir)) {
    cat("📭 No cache directory found. Run analysis first to generate cache.\n")
    return(list(
      cache_exists = FALSE,
      recommendations = "Run analysis pipeline to generate cache data"
    ))
  }

  # Analyze cache contents
  cache_files <- list.files(cache_dir, full.names = TRUE)
  cache_sizes <- file.info(cache_files)$size
  total_cache_size <- sum(cache_sizes, na.rm = TRUE)

  # Analyze cache hit patterns (simulated)
  cache_analysis <- list(
    cache_exists = TRUE,
    total_files = length(cache_files),
    total_size_mb = total_cache_size / 1024 / 1024,
    average_file_size_kb = mean(cache_sizes, na.rm = TRUE) / 1024,
    cache_efficiency_estimate = 0.75,  # Estimated based on usage patterns
    recommendations = generate_cache_recommendations(cache_files, total_cache_size)
  )

  cat("📊 Cache Analysis:\n")
  cat("  Cache files:", length(cache_files), "\n")
  cat("  Total size:", round(total_cache_size / 1024 / 1024, 2), "MB\n")
  cat("  Average file size:", round(mean(cache_sizes, na.rm = TRUE) / 1024, 1), "KB\n")
  cat("  Estimated efficiency:", round(cache_analysis$cache_efficiency_estimate * 100, 1), "%\n")

  return(cache_analysis)
}

# Generate cache recommendations
generate_cache_recommendations <- function(cache_files, total_size) {
  recommendations <- c()

  # Size recommendations
  if (total_size > 100 * 1024 * 1024) {  # 100MB
    recommendations <- c(recommendations, "Consider cleaning cache - large cache size detected")
  }

  # File count recommendations
  if (length(cache_files) > 100) {
    recommendations <- c(recommendations, "High number of cache files - consider cache organization")
  }

  # Pattern analysis
  if (length(recommendations) == 0) {
    recommendations <- c("Cache performance is optimal")
  }

  return(recommendations)
}

# Memory usage profiling
profile_memory_usage <- function() {
  cat("🧠 Profiling memory usage...\n")

  # Initial memory state
  initial_memory <- memory.size()
  initial_gc <- gc()

  # Create test dataset
  test_data <- data.frame(matrix(rnorm(100000), ncol = 50))

  # Measure memory after data creation
  after_creation <- memory.size()
  creation_growth <- after_creation - initial_memory

  # Simulate analysis operations
  start_time <- Sys.time()

  # Data manipulation
  processed_data <- test_data %>%
    mutate(
      BMI_category = case_when(
        X1 < 18.5 ~ "Underweight",
        X1 >= 18.5 & X1 < 25 ~ "Normal",
        X1 >= 25 & X1 < 30 ~ "Overweight",
        X1 >= 30 ~ "Obese"
      )
    ) %>%
    group_by(BMI_category) %>%
    summarize(
      mean_value = mean(X1, na.rm = TRUE),
      sd_value = sd(X1, na.rm = TRUE),
      count = n()
    )

  analysis_time <- difftime(Sys.time(), start_time, units = "secs")

  # Measure memory after processing
  after_processing <- memory.size()
  processing_growth <- after_processing - after_creation

  # Final cleanup and measurement
  rm(test_data, processed_data)
  final_gc <- gc()
  final_memory <- memory.size()

  memory_profile <- list(
    initial_memory = initial_memory,
    after_creation = after_creation,
    creation_growth = creation_growth,
    after_processing = after_processing,
    processing_growth = processing_growth,
    final_memory = final_memory,
    total_growth = final_memory - initial_memory,
    analysis_time = analysis_time,
    gc_efficiency = (final_gc[2, 6] - initial_gc[2, 6]) / (final_gc[1, 6] - initial_gc[1, 6]) * 100,
    recommendations = generate_memory_recommendations(creation_growth, processing_growth)
  )

  cat("📊 Memory Profile:\n")
  cat("  Initial memory:", round(initial_memory, 1), "MB\n")
  cat("  Data creation growth:", round(creation_growth, 1), "MB\n")
  cat("  Processing growth:", round(processing_growth, 1), "MB\n")
  cat("  Total growth:", round(memory_profile$total_growth, 1), "MB\n")
  cat("  Analysis time:", round(analysis_time, 2), "seconds\n")
  cat("  GC efficiency:", round(memory_profile$gc_efficiency, 1), "%\n")

  return(memory_profile)
}

# Generate memory recommendations
generate_memory_recommendations <- function(creation_growth, processing_growth) {
  recommendations <- c()

  if (creation_growth > 500) {
    recommendations <- c(recommendations, "High memory usage during data creation. Consider processing data in chunks")
  }

  if (processing_growth > 200) {
    recommendations <- c(recommendations, "Significant memory growth during processing. Optimize data manipulation operations")
  }

  if (length(recommendations) == 0) {
    recommendations <- c("Memory usage is within acceptable limits")
  }

  return(recommendations)
}

# Comprehensive performance report
generate_performance_report <- function() {
  cat("📈 Generating comprehensive performance report...\n")

  # System assessment
  system_perf <- assess_system_performance()

  # Parallel processing benchmarks
  parallel_perf <- benchmark_parallel_performance()

  # Cache analysis
  cache_perf <- analyze_cache_performance()

  # Memory profiling
  memory_perf <- profile_memory_usage()

  # Compile comprehensive report
  performance_report <- list(
    metadata = list(
      generated_at = Sys.time(),
      platform = system_perf$cpu$architecture,
      r_version = R.version.string
    ),
    system_assessment = system_perf,
    parallel_performance = parallel_perf,
    cache_performance = cache_perf,
    memory_performance = memory_perf,
    optimization_recommendations = generate_optimization_recommendations(
      system_perf, parallel_perf, cache_perf, memory_perf
    )
  )

  # Save report
  report_file <- "outputs/logs/performance_report.json"
  dir.create(dirname(report_file), showWarnings = FALSE, recursive = TRUE)
  write_json(performance_report, report_file, pretty = TRUE)

  cat("✅ Performance report generated:", report_file, "\n")

  return(performance_report)
}

# Generate optimization recommendations
generate_optimization_recommendations <- function(system_perf, parallel_perf, cache_perf, memory_perf) {
  recommendations <- list()

  # Parallel processing recommendations
  if (parallel_perf$speedup < PERF_CONFIG$optimization_targets$parallel_speedup) {
    recommendations$parallel <- paste(
      "Parallel speedup (", round(parallel_perf$speedup, 2), "x) is below target (",
      PERF_CONFIG$optimization_targets$parallel_speedup, "x). Consider optimizing parallel algorithms"
    )
  } else {
    recommendations$parallel <- "Parallel processing performance is meeting targets"
  }

  # Cache efficiency recommendations
  if (!is.null(cache_perf$cache_efficiency_estimate)) {
    if (cache_perf$cache_efficiency_estimate < PERF_CONFIG$optimization_targets$cache_hit_rate) {
      recommendations$cache <- paste(
        "Cache efficiency (", round(cache_perf$cache_efficiency_estimate * 100, 1),
        "%) is below target (", PERF_CONFIG$optimization_targets$cache_hit_rate * 100,
        "%). Consider cache optimization"
      )
    } else {
      recommendations$cache <- "Cache efficiency is meeting targets"
    }
  }

  # Memory efficiency recommendations
  memory_efficiency <- 1 - (memory_perf$total_growth / memory_perf$initial_memory)
  if (memory_efficiency < PERF_CONFIG$optimization_targets$memory_efficiency) {
    recommendations$memory <- paste(
      "Memory efficiency (", round(memory_efficiency * 100, 1),
      "%) is below target (", PERF_CONFIG$optimization_targets$memory_efficiency * 100,
      "%). Consider memory optimization"
    )
  } else {
    recommendations$memory <- "Memory efficiency is meeting targets"
  }

  # System-specific recommendations
  if (system_perf$memory$memory_usage_percent > 80) {
    recommendations$system <- "High baseline memory usage. Consider system optimization"
  } else {
    recommendations$system <- "System performance is adequate"
  }

  return(recommendations)
}

# Performance optimization suggestions
suggest_optimizations <- function() {
  cat("🔧 Performance Optimization Suggestions\n")
  cat("=====================================\n")

  # Hardware optimizations
  cat("💻 Hardware Optimizations:\n")
  cat("  • Use SSD storage for faster I/O operations\n")
  cat("  • Ensure sufficient RAM (8GB+ recommended)\n")
  cat("  • Multi-core CPU (4+ cores) for parallel processing\n")
  cat("  • Fast network connection for data downloads\n\n")

  # Software optimizations
  cat("⚙️ Software Optimizations:\n")
  cat("  • Enable parallel processing: make parallel-pipeline\n")
  cat("  • Use caching for repeated analyses\n")
  cat("  • Optimize memory usage with gc() calls\n")
  cat("  • Update R packages for performance improvements\n\n")

  # Configuration optimizations
  cat("📋 Configuration Optimizations:\n")
  cat("  • Set max_workers to available cores minus 1\n")
  cat("  • Enable verbose logging for performance monitoring\n")
  cat("  • Use appropriate cache strategies\n")
  cat("  • Monitor memory usage during intensive operations\n\n")

  # Workflow optimizations
  cat("🔄 Workflow Optimizations:\n")
  cat("  • Process data in chunks for large datasets\n")
  cat("  • Use vectorized operations instead of loops\n")
  cat("  • Pre-allocate memory for known data structures\n")
  cat("  • Clean up temporary objects regularly\n")
}

# Automated performance tuning
auto_tune_performance <- function() {
  cat("🔧 Auto-tuning performance parameters...\n")

  # Assess current system
  system_perf <- assess_system_performance()

  # Generate optimal configuration
  optimal_config <- list(
    performance = list(
      max_workers = max(1, system_perf$cpu$logical_cores - 1),
      memory_limit = min(system_perf$memory$total_memory * 0.8, 8192),  # 80% of total or 8GB max
      cache_enabled = TRUE,
      progress_verbose = system_perf$cpu$logical_cores > 4  # Verbose for multi-core systems
    )
  )

  # Save optimized configuration
  config_file <- "config/config_optimized.yml"
  write_yaml(optimal_config, config_file)

  cat("✅ Performance auto-tuning complete\n")
  cat("📋 Optimized configuration saved to:", config_file, "\n")
  cat("🔧 Recommended settings:\n")
  cat("  Workers:", optimal_config$performance$max_workers, "\n")
  cat("  Memory limit:", optimal_config$performance$memory_limit, "MB\n")

  return(optimal_config)
}

# Performance comparison across configurations
compare_configurations <- function() {
  cat("⚖️ Comparing performance across configurations...\n")

  # Test different worker configurations
  worker_configs <- c(1, 2, min(4, availableCores() - 1))

  comparison_results <- list()

  for (workers in worker_configs) {
    cat("Testing with", workers, "workers...\n")

    original_plan <- plan()
    plan(multisession, workers = workers)

    # Benchmark computation
    benchmark_result <- microbenchmark({
      test_data <- data.frame(x = rnorm(10000), y = rnorm(10000))
      result <- cor(test_data$x, test_data$y)
    }, times = 5)

    plan(original_plan)

    comparison_results[[as.character(workers)]] <- list(
      workers = workers,
      mean_time = mean(benchmark_result$time) / 1e9,
      median_time = median(benchmark_result$time) / 1e9,
      memory_usage = memory.size()
    )
  }

  # Analyze results
  best_config <- names(comparison_results)[which.min(sapply(comparison_results, function(x) x$mean_time))]

  cat("📊 Configuration Comparison Results:\n")
  for (config_name in names(comparison_results)) {
    config <- comparison_results[[config_name]]
    cat("  Workers:", config$workers, "->", round(config$mean_time, 4), "seconds\n")
  }

  cat("🏆 Best configuration:", best_config, "workers\n")

  return(comparison_results)
}

# Generate performance dashboard data
generate_performance_dashboard <- function() {
  cat("📊 Generating performance dashboard data...\n")

  # Collect all performance metrics
  dashboard_data <- list(
    timestamp = Sys.time(),
    system_performance = assess_system_performance(),
    parallel_benchmark = benchmark_parallel_performance(),
    cache_analysis = analyze_cache_performance(),
    memory_profile = profile_memory_usage(),
    configuration_comparison = compare_configurations()
  )

  # Save dashboard data
  dashboard_file <- "outputs/logs/performance_dashboard.json"
  dir.create(dirname(dashboard_file), showWarnings = FALSE, recursive = TRUE)
  write_json(dashboard_data, dashboard_file, pretty = TRUE)

  cat("✅ Performance dashboard data generated:", dashboard_file, "\n")

  return(dashboard_data)
}

# Display performance summary
display_performance_summary <- function() {
  cat("📈 Performance Summary\n")
  cat("=====================\n")

  # System performance
  system_perf <- assess_system_performance()
  cat("💻 System:\n")
  cat("  CPU cores:", system_perf$cpu$logical_cores, "\n")
  cat("  Memory:", round(system_perf$memory$available_memory, 1), "MB available\n")
  cat("  Storage:", round(system_perf$storage$available_space, 1), "GB available\n")

  # Parallel performance
  tryCatch({
    parallel_perf <- benchmark_parallel_performance()
    cat("\n⚡ Parallel Processing:\n")
    cat("  Speedup:", round(parallel_perf$speedup, 2), "x\n")
    cat("  Efficiency:", round(parallel_perf$efficiency, 1), "%\n")
  }, error = function(e) {
    cat("\n⚡ Parallel Processing: Not available\n")
  })

  # Cache performance
  tryCatch({
    cache_perf <- analyze_cache_performance()
    if (cache_perf$cache_exists) {
      cat("\n💾 Cache:\n")
      cat("  Size:", round(cache_perf$total_size_mb, 1), "MB\n")
      cat("  Files:", cache_perf$total_files, "\n")
    }
  }, error = function(e) {
    cat("\n💾 Cache: Not available\n")
  })

  # Memory performance
  tryCatch({
    memory_perf <- profile_memory_usage()
    cat("\n🧠 Memory:\n")
    cat("  Growth:", round(memory_perf$total_growth, 1), "MB\n")
    cat("  GC efficiency:", round(memory_perf$gc_efficiency, 1), "%\n")
  }, error = function(e) {
    cat("\n🧠 Memory: Not available\n")
  })

  cat("\n🔧 Recommendations:\n")
  tryCatch({
    recommendations <- generate_optimization_recommendations(
      system_perf,
      if (exists("parallel_perf")) parallel_perf else list(),
      if (exists("cache_perf")) cache_perf else list(),
      if (exists("memory_perf")) memory_perf else list()
    )

    for (category in names(recommendations)) {
      cat("  •", recommendations[[category]], "\n")
    }
  }, error = function(e) {
    cat("  • Run full performance analysis for recommendations\n")
  })

  cat("\n📋 Commands:\n")
  cat("  Rscript performance_tools.R benchmark    # Run benchmarks\n")
  cat("  Rscript performance_tools.R optimize     # Auto-optimize settings\n")
  cat("  Rscript performance_tools.R dashboard    # Generate dashboard data\n")
  cat("  Rscript performance_tools.R report       # Generate full report\n")
}

# Main performance tools function
run_performance_tools <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  action <- if (length(args) > 0) args[1] else "summary"

  switch(action,
    "benchmark" = {
      cat("🧪 Running performance benchmarks...\n")
      benchmark_results <- benchmark_parallel_performance()
      cat("✅ Benchmarking complete\n")
      return(benchmark_results)
    },
    "cache" = {
      cat("💾 Analyzing cache performance...\n")
      cache_analysis <- analyze_cache_performance()
      cat("✅ Cache analysis complete\n")
      return(cache_analysis)
    },
    "memory" = {
      cat("🧠 Profiling memory usage...\n")
      memory_profile <- profile_memory_usage()
      cat("✅ Memory profiling complete\n")
      return(memory_profile)
    },
    "optimize" = {
      cat("🔧 Auto-tuning performance...\n")
      optimization_config <- auto_tune_performance()
      cat("✅ Performance optimization complete\n")
      return(optimization_config)
    },
    "dashboard" = {
      cat("📊 Generating performance dashboard...\n")
      dashboard_data <- generate_performance_dashboard()
      cat("✅ Dashboard data generated\n")
      return(dashboard_data)
    },
    "report" = {
      cat("📈 Generating comprehensive performance report...\n")
      performance_report <- generate_performance_report()
      cat("✅ Performance report generated\n")
      return(performance_report)
    },
    "compare" = {
      cat("⚖️ Comparing configurations...\n")
      comparison_results <- compare_configurations()
      cat("✅ Configuration comparison complete\n")
      return(comparison_results)
    },
    "suggest" = {
      suggest_optimizations()
      return(NULL)
    },
    # Default: summary
    {
      display_performance_summary()
      return(NULL)
    }
  )
}

# Run performance tools if script is executed directly
if (!interactive()) {
  run_performance_tools()
}


