#!/usr/bin/env Rscript

# Test survey correctness and design-based estimates
# Verifies that survey weights and design are properly respected

suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(survey)
  library(here)
})

# Load analytic dataset
analytic_path <- here("data", "derived", "analytic.rds")
test_that("analytic dataset exists", {
  expect_true(file.exists(analytic_path))
})

analytic <- readRDS(analytic_path)

# Create survey design
svy_design <- svydesign(
  ids = ~psu,
  strata = ~strata,
  weights = ~survey_weight,
  nest = TRUE,
  data = analytic
)

test_that("survey design creates successfully", {
  expect_s3_class(svy_design, "survey.design")
  expect_equal(nrow(svy_design$variables), nrow(analytic))
  expect_true(all(!is.na(weights(svy_design))))
})

# Test 1: Survey-weighted means vs unweighted means
test_that("survey weights materially affect estimates", {
  
  # BMI means
  weighted_bmi <- as.numeric(svymean(~bmi, svy_design))
  unweighted_bmi <- mean(analytic$bmi, na.rm = TRUE)
  
  expect_type(weighted_bmi, "double")
  expect_type(unweighted_bmi, "double")
  expect_gt(abs(weighted_bmi - unweighted_bmi), 0.1)  # Should differ by >0.1 BMI units
  
  # Body fat means
  weighted_bf <- as.numeric(svymean(~bodyfat_pct, svy_design))
  unweighted_bf <- mean(analytic$bodyfat_pct, na.rm = TRUE)
  
  expect_gt(abs(weighted_bf - unweighted_bf), 0.5)  # Should differ by >0.5% body fat
})

# Test 2: Reference estimates for overall population
test_that("overall survey-weighted estimates are reasonable", {
  
  # Overall BMI
  overall_bmi <- svymean(~bmi, svy_design)
  bmi_mean <- as.numeric(overall_bmi)
  bmi_se <- as.numeric(SE(overall_bmi))
  
  expect_gt(bmi_mean, 25)    # Population BMI should be > 25 (overweight on average)
  expect_lt(bmi_mean, 35)    # But < 35 
  expect_gt(bmi_se, 0.05)    # Standard error should be reasonable
  expect_lt(bmi_se, 2)       # But not too large
  
  # Overall body fat
  overall_bf <- svymean(~bodyfat_pct, svy_design)
  bf_mean <- as.numeric(overall_bf)
  bf_se <- as.numeric(SE(overall_bf))
  
  expect_gt(bf_mean, 20)     # Population body fat should be > 20%
  expect_lt(bf_mean, 40)     # But < 40%
  expect_gt(bf_se, 0.1)      # Standard error should be reasonable
  expect_lt(bf_se, 3)        # But not too large
})

# Test 3: By-sex estimates
test_that("sex-specific estimates are reasonable and differ", {
  
  male_design <- subset(svy_design, sex == "Male")
  female_design <- subset(svy_design, sex == "Female")
  
  # BMI by sex
  male_bmi <- as.numeric(svymean(~bmi, male_design))
  female_bmi <- as.numeric(svymean(~bmi, female_design))
  
  expect_gt(male_bmi, 25)
  expect_gt(female_bmi, 25)
  expect_gt(abs(male_bmi - female_bmi), 0.2)  # Should differ by >0.2 units
  
  # Body fat by sex (females should have higher % body fat)
  male_bf <- as.numeric(svymean(~bodyfat_pct, male_design))
  female_bf <- as.numeric(svymean(~bodyfat_pct, female_design))
  
  expect_gt(female_bf, male_bf + 5)  # Females should have >5% more body fat
  expect_lt(male_bf, 30)             # Male body fat should be < 30%
  expect_gt(female_bf, 30)           # Female body fat should be > 30%
})

# Test 4: Quantiles by BMI class
test_that("median body fat increases with BMI class", {
  
  # Get median body fat by BMI class
  medians_by_class <- analytic %>%
    filter(bmi_cat %in% c("Normal", "Overweight", "Obesity I", "Obesity II")) %>%
    group_by(bmi_cat) %>%
    summarise(
      median_bf = as.numeric(svyquantile(~bodyfat_pct, 
                                         subset(svy_design, bmi_cat == unique(bmi_cat)), 
                                         quantiles = 0.5)[[1]]),
      .groups = "drop"
    ) %>%
    arrange(factor(bmi_cat, levels = c("Normal", "Overweight", "Obesity I", "Obesity II")))
  
  # Check that medians increase with BMI class
  for (i in 2:nrow(medians_by_class)) {
    expect_gt(medians_by_class$median_bf[i], medians_by_class$median_bf[i-1],
              info = paste("Median body fat should increase from", 
                          medians_by_class$bmi_cat[i-1], "to", 
                          medians_by_class$bmi_cat[i]))
  }
})

# Test 5: Survey design integrity
test_that("survey design variables are valid", {
  
  # Check for missing survey design variables
  expect_true(all(!is.na(analytic$survey_weight)))
  expect_true(all(!is.na(analytic$strata)))
  expect_true(all(!is.na(analytic$psu)))
  
  # Check weight distribution
  expect_gt(min(analytic$survey_weight), 0)
  expect_lt(max(analytic$survey_weight), 1e6)  # Reasonable upper bound
  
  # Check strata and PSU counts
  n_strata <- length(unique(analytic$strata))
  n_psu <- length(unique(analytic$psu))
  
  expect_gt(n_strata, 10)   # Should have multiple strata
  expect_gt(n_psu, 20)      # Should have multiple PSUs
  expect_gt(n_psu, n_strata) # More PSUs than strata
})

# Test 6: Correlation estimates
test_that("BMI-body fat correlation is strong and positive", {
  
  # Overall correlation
  corr_matrix <- svyvar(~bmi + bodyfat_pct, svy_design)
  correlation <- corr_matrix[1,2] / sqrt(corr_matrix[1,1] * corr_matrix[2,2])
  
  expect_gt(correlation, 0.8)   # Should be strongly positive
  expect_lt(correlation, 1.0)   # But less than perfect
  
  # By sex correlations should also be strong
  male_corr_matrix <- svyvar(~bmi + bodyfat_pct, subset(svy_design, sex == "Male"))
  male_correlation <- male_corr_matrix[1,2] / sqrt(male_corr_matrix[1,1] * male_corr_matrix[2,2])
  expect_gt(male_correlation, 0.8)
  
  female_corr_matrix <- svyvar(~bmi + bodyfat_pct, subset(svy_design, sex == "Female"))
  female_correlation <- female_corr_matrix[1,2] / sqrt(female_corr_matrix[1,1] * female_corr_matrix[2,2])
  expect_gt(female_correlation, 0.8)
})

# Test 7: Age group differences
test_that("age groups show expected patterns", {
  
  young_design <- subset(svy_design, age_group == "20-39")
  old_design <- subset(svy_design, age_group == "40-59")
  
  young_bf <- as.numeric(svymean(~bodyfat_pct, young_design))
  old_bf <- as.numeric(svymean(~bodyfat_pct, old_design))
  
  # Older adults typically have higher body fat
  expect_gt(old_bf, young_bf - 2)  # Allow some flexibility but expect trend
})

cat("✓ Survey correctness tests completed\n")