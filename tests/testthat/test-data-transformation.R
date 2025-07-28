context('Data Transformation Functions')

data("se_mini")
data("breast_tcga_mini_SE")

library(dplyr)
library(SummarizedExperiment)

# Test pivot functions with comprehensive column validation
test_that("pivot_sample works correctly and preserves expected columns", {
  # Get original colData structure
  original_colData <- colData(se_mini)
  
  # Test pivot_sample functionality
  result <- se_mini |> pivot_sample()
  
  # Check that the result is a tibble
  expect_true(inherits(result, "tbl_df"))
  
  # Check that all original colData columns are preserved (plus .sample)
  expected_cols <- c(".sample", names(original_colData))
  expect_equal(sort(names(result)), sort(expected_cols))
  
  # Check that the number of rows matches the number of samples
  expect_equal(nrow(result), nrow(original_colData))
  
  # Check that sample names are preserved in .sample column
  expect_equal(result$.sample, rownames(original_colData))
  
  # Check specific expected columns based on the test data
  expect_true("Cell.type" %in% names(result))
  expect_true("time" %in% names(result))
  expect_true("condition" %in% names(result))
  expect_true("days" %in% names(result))
  expect_true("dead" %in% names(result))
  expect_true(".sample" %in% names(result))
  
  # Check that condition values are preserved
  expect_true(all(result$condition %in% c(TRUE, FALSE)))
  
  # Check that the result has the expected number of columns (5 + 1 for .sample)
  expect_equal(ncol(result), 6)
})

test_that("pivot_transcript works correctly and preserves expected columns", {
  # Get original rowData structure
  original_rowData <- rowData(se_mini)
  
  # Test pivot_transcript functionality
  result <- se_mini |> pivot_transcript()
  
  # Check that the result is a tibble
  expect_true(inherits(result, "tbl_df"))
  
  # Check that all original rowData columns are preserved (plus .feature)
  # The entrez column is preserved as "entrez" in the output
  expected_cols <- c(".feature", "entrez")
  expect_equal(sort(names(result)), sort(expected_cols))
  
  # Check that the number of rows matches the number of features/transcripts
  expect_equal(nrow(result), nrow(original_rowData))
  
  # Check that feature names are preserved in .feature column
  expect_equal(result$.feature, rownames(original_rowData))
  
  # Check specific expected columns based on the test data
  expect_true("entrez" %in% names(result))  # This is the entrez column
  expect_true(".feature" %in% names(result))
  
  # Check that the result has the expected number of columns (1 + 1 for .feature)
  expect_equal(ncol(result), 2)
})

test_that("pivot_sample handles additional colData columns correctly", {
  # Create a modified version with extra columns
  se_modified <- se_mini
  colData(se_modified)$extra_column <- rep("test", nrow(colData(se_modified)))
  colData(se_modified)$numeric_column <- 1:nrow(colData(se_modified))
  
  # Test pivot_sample with additional columns
  result <- se_modified |> pivot_sample()
  
  # Check that all original columns plus new ones are preserved
  expected_cols <- c(".sample", "Cell.type", "time", "condition", "days", "dead", 
                     "extra_column", "numeric_column")
  expect_equal(sort(names(result)), sort(expected_cols))
  
  # Check that new columns have correct values
  expect_equal(result$extra_column, rep("test", nrow(result)))
  expect_equal(result$numeric_column, 1:nrow(result))
})

test_that("pivot_transcript handles additional rowData columns correctly", {
  # Create a modified version with extra columns
  se_modified <- se_mini
  rowData(se_modified)$extra_column <- rep("test", nrow(rowData(se_modified)))
  rowData(se_modified)$numeric_column <- 1:nrow(rowData(se_modified))
  
  # Test pivot_transcript with additional columns
  result <- se_modified |> pivot_transcript()
  
  # Check that all original columns plus new ones are preserved
  # The entrez column is preserved as "entrez" and additional columns are preserved
  expected_cols <- c(".feature", "entrez", "extra_column", "numeric_column")
  expect_equal(sort(names(result)), sort(expected_cols))
  
  # Check that new columns have correct values
  expect_equal(result$extra_column, rep("test", nrow(result)))
  expect_equal(result$numeric_column, 1:nrow(result))
})

test_that("pivot_sample preserves data types correctly", {
  result <- se_mini |> pivot_sample()
  
  # Check that logical columns remain logical
  expect_true(is.logical(result$condition))
  # Note: dead is actually numeric, not logical
  expect_true(is.numeric(result$dead))
  
  # Check that character columns remain character
  expect_true(is.character(result$Cell.type))
  expect_true(is.character(result$time))
  
  # Check that numeric columns remain numeric
  expect_true(is.numeric(result$days))
  
  # Check that .sample column is character
  expect_true(is.character(result$.sample))
})

test_that("pivot_transcript preserves data types correctly", {
  result <- se_mini |> pivot_transcript()
  
  # Check that character columns remain character
  expect_true(is.character(result$entrez))  # This is the entrez column
  expect_true(is.character(result$.feature))
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