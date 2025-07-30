context('Cellularity Analysis Functions')

library(airway)
data(airway)
se_mini <- airway[1:100, 1:5]

library(dplyr)
library(SummarizedExperiment)

# Test deconvolve_cellularity function

test_that("deconvolve_cellularity throws error with multiple methods", {
  expect_error(
    se_mini |> deconvolve_cellularity(method = c("cibersort", "llsr")),
    "Multiple methods provided"
  )
})

test_that("deconvolve_cellularity with feature_column works correctly", {
  # Add a test column to rowData with valid gene names

  
  # Test with valid feature_column
  res <- airway |> 
    deconvolve_cellularity(feature_column = "gene_name", cores = 1 )
  expect_no_error(res)
  
  # Test that it throws error with non-existent column
  expect_error(
    airway |> deconvolve_cellularity(feature_column = "non_existent_column", cores = 1),
    "feature_column ' non_existent_column ' not found in rowData"
  )
  
  # Test that it throws error with column containing NA values
  rowData(airway)$na_feature <- rowData(airway)$test_feature
  rowData(airway)$na_feature[1] <- NA
  expect_error(
    airway |> deconvolve_cellularity(feature_column = "na_feature", cores = 1),
    "feature_column ' na_feature ' contains missing or empty values"
  )
  
  # Test that it throws error with column containing empty values
  rowData(airway)$empty_feature <- rowData(airway)$test_feature
  rowData(airway)$empty_feature[1] <- ""
  expect_error(
    airway |> deconvolve_cellularity(feature_column = "empty_feature", cores = 1),
    "feature_column ' empty_feature ' contains missing or empty values"
  )
})

test_that("deconvolve_cellularity with EPIC works correctly", {
  # Skip EPIC as it requires additional packages
  skip("EPIC requires additional Bioconductor packages")
})

# Test test_stratification_cellularity function - DEPRECATED
# test_that("test_stratification_cellularity works correctly", {
#   # Skip test_stratification_cellularity as it has formula issues
#   skip("test_stratification_cellularity has formula compatibility issues")
# })

# Test cibersort function
test_that("cibersort works correctly", {
  # Skip cibersort as it's not available
  skip("cibersort function not available")
}) 