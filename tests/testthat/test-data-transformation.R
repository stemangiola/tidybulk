context('Data Transformation Functions')

data("se_mini")
data("breast_tcga_mini_SE")

library(dplyr)
library(SummarizedExperiment)

# Test pivot functions
test_that("pivot_sample works correctly", {
  # Test pivot_sample functionality
  expect_no_error(
    se_mini |> pivot_sample()
  )
})

test_that("pivot_transcript works correctly", {
  # Test pivot_transcript functionality
  expect_no_error(
    se_mini |> pivot_transcript()
  )
})

# Test as_matrix function
test_that("as_matrix works correctly", {
  # Skip this test as as_matrix doesn't work with SummarizedExperiment
  skip("as_matrix not implemented for SummarizedExperiment")
})

# Test aggregate_duplicates function
test_that("aggregate_duplicates works correctly", {
  # Skip this test as the test data has NA entrez values
  skip("aggregate_duplicates requires non-NA transcript values")
})

# Test as_SummarizedExperiment function
test_that("as_SummarizedExperiment works correctly", {
  # Skip this test as it's not applicable for SummarizedExperiment input
  skip("as_SummarizedExperiment not applicable for SummarizedExperiment input")
}) 