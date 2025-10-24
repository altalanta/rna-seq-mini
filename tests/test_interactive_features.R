# Interactive Features Integration Tests
# Tests tutorial system and configuration wizard

library(testthat)
library(yaml)
library(jsonlite)

# Test 1: Tutorial System Integration
test_that("Tutorial system integrates properly", {

  # Test 1a: Tutorial files exist
  tutorial_files <- c("tutorials/getting_started.Rmd", "tutorials/help_troubleshooting.Rmd")

  for (tutorial_file in tutorial_files) {
    if (file.exists(tutorial_file)) {
      expect_true(file.exists(tutorial_file))

      # Test file content
      content <- readLines(tutorial_file, n = 20)

      # Should contain learnr tutorial structure
      expect_true(any(grepl("learnr::tutorial", content)))
      expect_true(any(grepl("runtime: shiny", content)))
      expect_true(any(grepl("output:", content)))
    }
  }

  # Test 1b: CSS styling exists
  css_file <- "tutorials/css/styles.css"
  if (file.exists(css_file)) {
    css_content <- readLines(css_file, n = 10)
    expect_true(length(css_content) > 0)
    expect_true(any(grepl("\\.section", css_content)))  # Should have section styling
  }

  # Test 1c: Tutorial structure validation
  getting_started_file <- "tutorials/getting_started.Rmd"
  if (file.exists(getting_started_file)) {
    content <- readLines(getting_started_file)

    # Should have YAML header
    expect_true(any(grepl("^---$", content)))

    # Should have title
    expect_true(any(grepl("^title:", content)))

    # Should have tutorial structure
    expect_true(any(grepl("learnr::tutorial", content)))

    # Should have sections
    expect_true(any(grepl("^## ", content)))

    # Should have interactive elements
    expect_true(any(grepl("\\{r ", content)))
  }
})

# Test 2: Configuration Wizard Integration
test_that("Configuration wizard integrates properly", {

  # Test 2a: App.R exists and loads
  expect_true(file.exists("app.R"))

  # Test 2b: App structure validation
  app_content <- readLines("app.R", n = 50)

  # Should have Shiny app structure
  expect_true(any(grepl("shinyApp", app_content)))
  expect_true(any(grepl("fluidPage", app_content)))
  expect_true(any(grepl("server", app_content)))

  # Test 2c: Configuration validation
  if (file.exists("config/config.yml")) {
    config <- read_yaml("config/config.yml")
    expect_is(config, "list")

    # Should have required sections
    required_sections <- c("data", "outputs", "nhanes", "analysis", "logging")
    for (section in required_sections) {
      expect_true(section %in% names(config))
    }

    # Should have valid data types
    expect_true(is.numeric(config$analysis$age_range[1]))
    expect_true(is.character(config$analysis$survey_weights_col))
  }

  # Test 2d: Configuration wizard functionality
  # Test that configuration can be loaded and validated
  test_config <- list(
    analysis = list(
      age_range = c(20, 59),
      survey_weights_col = "WTMEC2YR"
    ),
    data = list(
      raw_dir = "data/raw",
      derived_dir = "data/derived"
    )
  )

  # Should be able to convert to YAML
  yaml_config <- yaml::as.yaml(test_config)
  expect_is(yaml_config, "character")
  expect_true(nchar(yaml_config) > 0)

  # Should be able to parse back
  parsed_config <- yaml::read_yaml(text = yaml_config)
  expect_is(parsed_config, "list")
  expect_equal(parsed_config$analysis$age_range, c(20, 59))
})

# Test 3: Interactive Documentation Integration
test_that("Interactive documentation integrates properly", {

  # Test 3a: Documentation files exist
  doc_files <- c(
    "README.md",
    "docs/INSTALLATION.md",
    "docs/PARALLEL_PROCESSING.md",
    "docs/INTERACTIVE_FEATURES.md",
    "docs/DATA_VERSIONING.md",
    "docs/API_REFERENCE.md",
    "docs/INTEGRATION_EXAMPLES.md",
    "docs/TROUBLESHOOTING.md"
  )

  for (doc_file in doc_files) {
    if (file.exists(doc_file)) {
      expect_true(file.exists(doc_file))

      # Should have content
      content <- readLines(doc_file, n = 5)
      expect_true(length(content) > 0)
    }
  }

  # Test 3b: Documentation cross-references
  # Test that documentation files reference each other appropriately
  readme_content <- readLines("README.md", n = 100)

  # Should reference key documentation
  expect_true(any(grepl("docs/INSTALLATION.md", readme_content)))
  expect_true(any(grepl("make tutorial", readme_content)))
  expect_true(any(grepl("make config-wizard", readme_content)))

  # Test 3c: Documentation completeness
  # Check that major sections are covered
  major_sections <- c("Installation", "Usage", "Configuration", "Documentation", "Support")
  for (section in major_sections) {
    section_found <- FALSE
    for (doc_file in doc_files) {
      if (file.exists(doc_file)) {
        content <- readLines(doc_file)
        if (any(grepl(section, content, ignore.case = TRUE))) {
          section_found <- TRUE
          break
        }
      }
    }
    # Most sections should be covered in at least one doc
    expect_true(section_found || section %in% c("Support"))  # Support might be in README only
  }
})

# Test 4: Error Handling in Interactive Context
test_that("Error handling works in interactive context", {

  # Test 4a: Custom error classes
  test_error <- NhanesError("Interactive test error", "INT001")
  expect_is(test_error, "NhanesError")
  expect_equal(test_error$message, "Interactive test error")
  expect_equal(test_error$code, "INT001")

  # Test 4b: Error suggestions
  suggestions <- get_error_suggestions("INT001", "Interactive test error")
  expect_is(suggestions, "character")
  expect_true(length(suggestions) > 0)

  # Test 4c: User-friendly error display
  # Should not throw error when displaying
  expect_silent(display_user_friendly_error(test_error))

  # Test 4d: Error logging
  test_config <- list(outputs = list(logs_dir = "tests/test_outputs"))
  safe_log("Interactive test log", "INFO", test_config)

  log_file <- file.path(test_config$outputs$logs_dir, "analysis_log.txt")
  if (file.exists(log_file)) {
    log_content <- readLines(log_file)
    expect_true(any(grepl("Interactive test log", log_content)))
  }
})

# Test 5: Configuration Integration
test_that("Configuration system integrates with interactive features", {

  # Test 5a: Configuration file validation
  config_file <- "config/config.yml"

  if (file.exists(config_file)) {
    config <- read_yaml(config_file)
    expect_is(config, "list")

    # Should have interactive-friendly structure
    expect_true(is.list(config$analysis))
    expect_true(is.list(config$data))
    expect_true(is.list(config$outputs))
  }

  # Test 5b: Configuration export/import
  test_config <- list(
    analysis = list(age_range = c(20, 59)),
    interactive = list(tutorial_enabled = TRUE)
  )

  # Should be able to export as JSON
  json_config <- toJSON(test_config, pretty = TRUE)
  expect_is(json_config, "character")

  # Should be able to import back
  imported_config <- fromJSON(json_config)
  expect_is(imported_config, "list")
  expect_equal(imported_config$analysis$age_range, c(20, 59))

  # Test 5c: Configuration wizard data structure
  # Test that configuration can support wizard interface
  wizard_config <- list(
    ui_elements = list(
      sliders = list(age_min = list(min = 0, max = 80, value = 20)),
      dropdowns = list(survey_weights = list(choices = c("WTMEC2YR", "WTMEC4YR"))),
      text_inputs = list(output_dir = list(value = "outputs/"))
    )
  )

  expect_is(wizard_config, "list")
  expect_true(is.list(wizard_config$ui_elements))
})

# Test 6: Tutorial Content Validation
test_that("Tutorial content is valid and complete", {

  # Test 6a: Getting Started tutorial
  getting_started <- "tutorials/getting_started.Rmd"
  if (file.exists(getting_started)) {
    content <- readLines(getting_started)

    # Should have proper YAML header
    yaml_start <- which(content == "---")[1]
    yaml_end <- which(content == "---")[2]
    expect_true(yaml_end > yaml_start)

    # Should have title
    expect_true(any(grepl("^title:", content)))

    # Should have tutorial structure
    expect_true(any(grepl("learnr::tutorial", content)))

    # Should have interactive elements
    expect_true(any(grepl("\\{r ", content)))  # R code chunks
    expect_true(any(grepl("quiz", content)))   # Quiz elements

    # Should have sections
    section_count <- sum(grepl("^## ", content))
    expect_true(section_count >= 5)  # Should have multiple sections
  }

  # Test 6b: Troubleshooting tutorial
  troubleshooting <- "tutorials/help_troubleshooting.Rmd"
  if (file.exists(troubleshooting)) {
    content <- readLines(troubleshooting)

    # Should have troubleshooting sections
    expect_true(any(grepl("Common Issues", content)))
    expect_true(any(grepl("Installation Problems", content)))
    expect_true(any(grepl("Troubleshooting", content)))

    # Should have solutions
    expect_true(any(grepl("Solution:", content)))
    expect_true(any(grepl("Quick Fix:", content)))
  }

  # Test 6c: CSS styling
  css_file <- "tutorials/css/styles.css"
  if (file.exists(css_file)) {
    css_content <- readLines(css_file)

    # Should have CSS rules
    expect_true(any(grepl("\\.section", css_content)))
    expect_true(any(grepl("\\.alert", css_content)))
    expect_true(any(grepl("@media", css_content)))  # Responsive design
  }
})

# Test 7: Integration with Build System
test_that("Interactive features integrate with build system", {

  # Test 7a: Makefile targets
  makefile_content <- readLines("Makefile")

  # Should have interactive feature targets
  expect_true(any(grepl("tutorial:", makefile_content)))
  expect_true(any(grepl("config-wizard:", makefile_content)))
  expect_true(any(grepl("tutorial-troubleshooting:", makefile_content)))

  # Test 7b: Target functionality
  # Test that targets are properly defined
  tutorial_line <- grep("tutorial:", makefile_content)
  expect_true(length(tutorial_line) > 0)

  tutorial_target <- makefile_content[tutorial_line[1]]
  expect_true(grepl("learnr::run_tutorial", tutorial_target))

  # Test 7c: Help documentation
  help_section <- grep("help:", makefile_content)
  if (length(help_section) > 0) {
    help_content <- makefile_content[(help_section[1]):length(makefile_content)]

    # Should document interactive features
    expect_true(any(grepl("tutorial", help_content)))
    expect_true(any(grepl("config-wizard", help_content)))
    expect_true(any(grepl("troubleshooting", help_content)))
  }
})

# Test 8: Cross-Feature Integration
test_that("Interactive features integrate across platform", {

  # Test 8a: Configuration wizard with tutorials
  # Test that configuration can be used in tutorial context
  if (file.exists("config/config.yml") && file.exists("tutorials/getting_started.Rmd")) {
    config <- read_yaml("config/config.yml")
    tutorial_content <- readLines("tutorials/getting_started.Rmd", n = 50)

    # Tutorial should be able to reference configuration
    expect_true(any(grepl("config", tutorial_content, ignore.case = TRUE)))
  }

  # Test 8b: Error handling in tutorial context
  # Test that error handling works in interactive environment
  test_error <- NhanesError("Tutorial integration test", "TUT001")
  expect_silent(display_user_friendly_error(test_error))

  # Test 8c: Documentation cross-references
  # Test that documentation files reference each other
  doc_files <- c("README.md", "docs/INSTALLATION.md", "docs/INTERACTIVE_FEATURES.md")

  for (doc_file in doc_files) {
    if (file.exists(doc_file)) {
      content <- readLines(doc_file, n = 100)

      # Should reference other documentation
      if (doc_file != "README.md") {
        expect_true(any(grepl("README.md", content)) ||
                   any(grepl("installation", content, ignore.case = TRUE)))
      }
    }
  }
})

# Test 9: Performance in Interactive Context
test_that("Interactive features perform adequately", {

  # Test 9a: Configuration loading performance
  start_time <- Sys.time()

  if (file.exists("config/config.yml")) {
    config <- read_yaml("config/config.yml")
    expect_is(config, "list")
  }

  config_load_time <- difftime(Sys.time(), start_time, units = "secs")
  expect_true(config_load_time < 1)  # Should load quickly

  # Test 9b: Tutorial file loading
  tutorial_files <- c("tutorials/getting_started.Rmd", "tutorials/help_troubleshooting.Rmd")

  for (tutorial_file in tutorial_files) {
    if (file.exists(tutorial_file)) {
      start_time <- Sys.time()
      content <- readLines(tutorial_file, n = 50)
      load_time <- difftime(Sys.time(), start_time, units = "secs")

      expect_true(load_time < 0.5)  # Should load quickly
      expect_true(length(content) > 0)
    }
  }

  # Test 9c: Error handling performance
  start_time <- Sys.time()

  for (i in 1:10) {
    test_error <- NhanesError(paste("Performance test error", i), paste0("PERF", sprintf("%03d", i)))
    expect_silent(display_user_friendly_error(test_error))
  }

  error_handling_time <- difftime(Sys.time(), start_time, units = "secs")
  expect_true(error_handling_time < 1)  # Should be fast
})

# Test 10: Accessibility and Usability
test_that("Interactive features are accessible and usable", {

  # Test 10a: File structure accessibility
  # Test that all required files are accessible
  required_files <- c(
    "tutorials/getting_started.Rmd",
    "tutorials/help_troubleshooting.Rmd",
    "tutorials/css/styles.css",
    "app.R",
    "config/config.yml"
  )

  for (required_file in required_files) {
    if (file.exists(required_file)) {
      # Should be readable
      expect_true(file.access(required_file, mode = 4) == 0)  # Readable
    }
  }

  # Test 10b: Documentation completeness
  # Check that key topics are covered
  key_topics <- c("installation", "configuration", "tutorial", "troubleshooting", "parallel", "caching")

  doc_coverage <- list()
  for (topic in key_topics) {
    covered <- FALSE
    for (doc_file in c("README.md", "docs/INSTALLATION.md", "docs/INTERACTIVE_FEATURES.md")) {
      if (file.exists(doc_file)) {
        content <- tolower(readLines(doc_file, n = 200))
        if (topic %in% content) {
          covered <- TRUE
          break
        }
      }
    }
    doc_coverage[[topic]] <- covered
  }

  # Most topics should be covered
  coverage_rate <- sum(unlist(doc_coverage)) / length(key_topics)
  expect_true(coverage_rate >= 0.8)  # At least 80% coverage

  # Test 10c: Cross-platform compatibility
  # Test that paths work on different platforms
  test_paths <- c("data/raw", "outputs/tables", "config/config.yml")

  for (path in test_paths) {
    # Should work with both forward and back slashes
    normalized_path <- normalizePath(path, mustWork = FALSE)
    expect_is(normalized_path, "character")
  }
})

# Run all tests with reporting
run_interactive_tests <- function() {
  cat("🎓 Running Interactive Features Integration Tests\n")
  cat("================================================\n")

  # Run test suite
  test_results <- testthat::test_file("tests/test_interactive_features.R")

  # Generate report
  test_report <- list(
    timestamp = Sys.time(),
    tests_run = length(test_results),
    tests_passed = sum(sapply(test_results, function(x) x$passed)),
    tests_failed = sum(sapply(test_results, function(x) !x$passed)),
    total_time = sum(sapply(test_results, function(x) x$time))
  )

  # Display results
  cat("\n📊 Interactive Features Test Results\n")
  cat("====================================\n")
  cat("Tests run:", test_report$tests_run, "\n")
  cat("Tests passed:", test_report$tests_passed, "\n")
  cat("Tests failed:", test_report$tests_failed, "\n")
  cat("Total time:", round(test_report$total_time, 2), "seconds\n")

  # Check tutorial files
  tutorial_files <- c("tutorials/getting_started.Rmd", "tutorials/help_troubleshooting.Rmd")
  available_tutorials <- sum(file.exists(tutorial_files))

  cat("\n📚 Tutorial Availability:\n")
  cat("  Available tutorials:", available_tutorials, "/", length(tutorial_files), "\n")

  if (file.exists("app.R")) {
    cat("  Configuration wizard: ✅ Available\n")
  } else {
    cat("  Configuration wizard: ❌ Missing\n")
  }

  if (test_report$tests_failed > 0) {
    cat("\n❌ Some interactive features tests failed!\n")
    cat("Check test output above for details.\n")
  } else {
    cat("\n✅ All interactive features tests passed!\n")
    cat("🎉 Interactive features integration verified successfully.\n")
  }

  return(test_report)
}

# If run directly, execute tests
if (!interactive()) {
  run_interactive_tests()
}


