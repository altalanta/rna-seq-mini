# NHANES BMI Body Fat Analysis API Tests
# Automated tests for the REST API functionality

library(testthat)
library(httr)
library(jsonlite)

# Skip tests if httr is not available
skip_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    skip(paste(pkg, "not available for testing"))
  }
}

# Test API functionality
test_api <- function(port = 8000) {
  base_url <- paste0("http://localhost:", port)

  context("NHANES BMI Body Fat API")

  test_that("API health endpoint works", {
    skip_if_not_installed("httr")
    response <- GET(paste0(base_url, "/health"))
    expect_equal(status_code(response), 200)
    content <- content(response, "parsed")
    expect_true(content$status == "healthy")
  })

  test_that("API info endpoint works", {
    skip_if_not_installed("httr")
    response <- GET(paste0(base_url, "/api/info"))
    expect_equal(status_code(response), 200)
    content <- content(response, "parsed")
    expect_true(!is.null(content$name))
    expect_true(!is.null(content$version))
  })

  test_that("Correlations endpoint works", {
    skip_if_not_installed("httr")
    response <- GET(paste0(base_url, "/api/correlations"))
    expect_equal(status_code(response), 200)
    content <- content(response, "parsed")
    expect_true(!is.null(content$data))
    expect_true(content$count > 0)
  })

  test_that("Correlations filtering works", {
    skip_if_not_installed("httr")
    response <- GET(paste0(base_url, "/api/correlations?group=Male"))
    expect_equal(status_code(response), 200)
    content <- content(response, "parsed")
    expect_true(all(content$data$group == "Male"))
  })

  test_that("Body fat BMI endpoint works", {
    skip_if_not_installed("httr")
    response <- GET(paste0(base_url, "/api/bodyfat/bmi/Normal"))
    expect_equal(status_code(response), 200)
    content <- content(response, "parsed")
    expect_true(!is.null(content$data))
    expect_true(all(content$data$bmi_cat == "Normal"))
  })

  test_that("Population endpoint works", {
    skip_if_not_installed("httr")
    response <- GET(paste0(base_url, "/api/population"))
    expect_equal(status_code(response), 200)
    content <- content(response, "parsed")
    expect_true(!is.null(content$data))
    expect_true(content$count > 0)
  })

  test_that("Statistics endpoint works", {
    skip_if_not_installed("httr")
    response <- GET(paste0(base_url, "/api/statistics"))
    expect_equal(status_code(response), 200)
    content <- content(response, "parsed")
    expect_true(!is.null(content$data_summary))
  })
}

# Run tests if called directly
if (!interactive()) {
  # Assume API is running on port 8000 for testing
  test_api(8000)
}

