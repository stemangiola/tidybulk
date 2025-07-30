context('Test Setup and Common Functions')

# Load test data
library(airway)
data(airway)
se_mini <- airway[1:100, 1:5]
data("breast_tcga_mini_SE")

# Load required libraries
library(dplyr)
library(SummarizedExperiment)
library(testthat)

# Helper function to check if SummarizedExperiment has expected assays
check_assays <- function(se, expected_assays) {
  actual_assays <- names(SummarizedExperiment::assays(se))
  expect_true(all(expected_assays %in% actual_assays))
}

# Helper function to check if colData has expected columns
check_colData_columns <- function(se, expected_columns) {
  actual_columns <- names(SummarizedExperiment::colData(se))
  expect_true(all(expected_columns %in% actual_columns))
}

# Helper function to check if rowData has expected columns
check_rowData_columns <- function(se, expected_columns) {
  actual_columns <- names(SummarizedExperiment::rowData(se))
  expect_true(all(expected_columns %in% actual_columns))
}

# Helper function to check if result has fewer rows than input
check_reduced_rows <- function(result, input) {
  expect_true(nrow(result) <= nrow(input))
}

# Helper function to check if result has fewer columns than input
check_reduced_cols <- function(result, input) {
  expect_true(ncol(result) <= ncol(input))
} 