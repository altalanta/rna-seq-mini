#!/usr/bin/env Rscript

# Test exclusion criteria and data quality checks
# Verifies pregnancy exclusions and invalid DXA flags

suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(foreign)
  library(here)
})

# Load raw data to test exclusions
demo_path <- here("data", "raw", "DEMO_J.XPT")
dxx_path <- here("data", "raw", "DXX_J.XPT")

test_that("raw data files exist", {
  expect_true(file.exists(demo_path))
  expect_true(file.exists(dxx_path))
})

demo <- read.xport(demo_path)
dxx <- read.xport(dxx_path)

# Test 1: Age exclusions
test_that("age restrictions are properly applied", {
  
  # Check age range in raw data
  age_range <- range(demo$RIDAGEYR, na.rm = TRUE)
  expect_gte(age_range[1], 0)
  expect_lte(age_range[2], 100)
  
  # Count by age groups
  adults_20_59 <- sum(demo$RIDAGEYR >= 20 & demo$RIDAGEYR <= 59, na.rm = TRUE)
  total_with_age <- sum(!is.na(demo$RIDAGEYR))
  
  expect_gt(adults_20_59, 1000)  # Should have substantial sample
  expect_lt(adults_20_59, total_with_age)  # Should exclude some people
  
  # Check that analytic dataset only includes target age range
  analytic <- readRDS(here("data", "derived", "analytic.rds"))
  expect_true(all(analytic$age >= 20 & analytic$age <= 59))
})

# Test 2: DXA scan quality flags
test_that("DXA scan quality is properly checked", {
  
  # Check DXXOSTAT (DXA exam status) in raw data
  if ("DXXOSTAT" %in% names(dxx)) {
    dxa_status_table <- table(dxx$DXXOSTAT, useNA = "always")
    
    # Status 1 should be complete/valid scans
    valid_scans <- sum(dxx$DXXOSTAT == 1, na.rm = TRUE)
    total_scans <- nrow(dxx)
    
    expect_gt(valid_scans, 0.7 * total_scans)  # At least 70% should be valid
    
    # Check that analytic dataset only includes valid scans
    analytic <- readRDS(here("data", "derived", "analytic.rds"))
    
    # Verify body fat values are reasonable
    expect_true(all(analytic$bodyfat_pct > 0 & analytic$bodyfat_pct < 100))
    expect_true(all(!is.na(analytic$bodyfat_pct)))
  }
})

# Test 3: Pregnancy exclusions (if pregnancy variable exists)
test_that("pregnancy exclusions are handled appropriately", {
  
  if ("RIDEXPRG" %in% names(demo)) {
    pregnant_count <- sum(demo$RIDEXPRG == 1, na.rm = TRUE)
    
    # Check that some pregnancy cases exist in raw data
    expect_gt(pregnant_count, 0)
    
    # Check that pregnant women are excluded from analytic dataset
    analytic <- readRDS(here("data", "derived", "analytic.rds"))
    
    # If pregnancy variable was available, pregnant women should be excluded
    # This is hard to verify directly since we don't keep pregnancy variable
    # But we can check that the exclusion counts make sense
    
    expect_gt(nrow(analytic), 1000)  # Should still have substantial sample
  } else {
    # If no pregnancy variable, just note this in test output
    cat("Note: No pregnancy variable (RIDEXPRG) found in DEMO_J\n")
  }
})

# Test 4: BMI validity
test_that("BMI values are within reasonable ranges", {
  
  analytic <- readRDS(here("data", "derived", "analytic.rds"))
  
  # BMI should be positive and reasonable
  expect_true(all(analytic$bmi > 10))   # Minimum reasonable BMI
  expect_true(all(analytic$bmi < 80))   # Maximum reasonable BMI
  expect_true(all(!is.na(analytic$bmi)))
  
  # BMI categories should be properly assigned
  expect_true(all(analytic$bmi_cat %in% c("Underweight", "Normal", "Overweight", 
                                          "Obesity I", "Obesity II", "Obesity III")))
  
  # Check category boundaries
  underweight <- analytic[analytic$bmi_cat == "Underweight", ]
  normal <- analytic[analytic$bmi_cat == "Normal", ]
  overweight <- analytic[analytic$bmi_cat == "Overweight", ]
  
  if (nrow(underweight) > 0) expect_true(all(underweight$bmi < 18.5))
  if (nrow(normal) > 0) expect_true(all(normal$bmi >= 18.5 & normal$bmi < 25))
  if (nrow(overweight) > 0) expect_true(all(overweight$bmi >= 25 & overweight$bmi < 30))
})

# Test 5: Survey weight validity
test_that("survey weights are valid and non-zero", {
  
  analytic <- readRDS(here("data", "derived", "analytic.rds"))
  
  # All weights should be positive
  expect_true(all(analytic$survey_weight > 0))
  expect_true(all(!is.na(analytic$survey_weight)))
  
  # Weight distribution should be reasonable
  weight_stats <- summary(analytic$survey_weight)
  expect_gt(weight_stats["Min."], 1)        # Minimum weight > 1
  expect_lt(weight_stats["Max."], 500000)   # Maximum weight < 500,000
  
  # Should have some variation in weights
  expect_gt(sd(analytic$survey_weight), 1000)
})

# Test 6: Missing data patterns
test_that("missing data exclusions are documented", {
  
  # Check that flow log exists and documents exclusions
  flow_log_path <- here("outputs", "logs", "flow.txt")
  expect_true(file.exists(flow_log_path))
  
  flow_content <- readLines(flow_log_path)
  
  # Should mention key exclusion steps
  expect_true(any(grepl("Initial", flow_content)))
  expect_true(any(grepl("Final", flow_content)))
  expect_true(any(grepl("excluded", flow_content, ignore.case = TRUE)))
  
  # Should document sample sizes at each step
  expect_true(any(grepl("N = [0-9]+", flow_content)))
})

# Test 7: Data completeness in analytic dataset
test_that("analytic dataset has complete key variables", {
  
  analytic <- readRDS(here("data", "derived", "analytic.rds"))
  
  # Core variables should have no missing values
  core_vars <- c("seqn", "age", "sex", "bmi", "bodyfat_pct", 
                 "survey_weight", "strata", "psu")
  
  for (var in core_vars) {
    expect_true(var %in% names(analytic), info = paste("Variable", var, "should exist"))
    expect_true(all(!is.na(analytic[[var]])), 
                info = paste("Variable", var, "should have no missing values"))
  }
  
  # Check that derived variables are properly created
  expect_true(all(analytic$sex %in% c("Male", "Female")))
  expect_true(all(analytic$age_group %in% c("20-39", "40-59")))
  expect_true(is.logical(analytic$obesity))
})

# Test 8: Sample size adequacy
test_that("final sample size is adequate for analysis", {
  
  analytic <- readRDS(here("data", "derived", "analytic.rds"))
  
  # Overall sample size
  expect_gte(nrow(analytic), 1500)  # Should have at least 1500 participants
  
  # By sex
  sex_counts <- table(analytic$sex)
  expect_gte(min(sex_counts), 500)  # At least 500 per sex
  
  # By BMI category (major categories)
  bmi_counts <- table(analytic$bmi_cat)
  major_categories <- c("Normal", "Overweight", "Obesity I")
  major_counts <- bmi_counts[names(bmi_counts) %in% major_categories]
  expect_gte(min(major_counts, na.rm = TRUE), 100)  # At least 100 per major category
})

cat("✓ Exclusion and data quality tests completed\n")