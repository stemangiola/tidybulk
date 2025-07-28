context('Abundance Normalization Functions')

data("se_mini")
data("breast_tcga_mini_SE")

library(dplyr)
library(SummarizedExperiment)

# Test scale_abundance function
test_that("scale_abundance works correctly", {
  res <- se_mini |> identify_abundant() |> scale_abundance()
  
  expect_equal(
    names(SummarizedExperiment::assays(res)),
    c("counts", "counts_scaled")
  )
})

test_that("scale_abundance with subset works correctly", {
  res <- se_mini |> identify_abundant() |> scale_abundance(
    .subset_for_scaling = .abundant & grepl("^A", .feature)
  )
  
  expect_true("counts_scaled" %in% names(SummarizedExperiment::assays(res)))
})

# Test quantile_normalise_abundance function
test_that("quantile_normalise_abundance works correctly", {
  res <- se_mini |> quantile_normalise_abundance()
  
  # Check if any normalized assay exists
  assay_names <- names(SummarizedExperiment::assays(res))
  expect_true(length(assay_names) > 1 || any(grepl("normalised", assay_names)))
})

test_that("quantile_normalise_abundance with preprocessCore works correctly", {
  res <- se_mini |> quantile_normalise_abundance(method = "preprocesscore_normalize_quantiles_use_target")
  
  # Check if any normalized assay exists
  assay_names <- names(SummarizedExperiment::assays(res))
  expect_true(length(assay_names) > 1 || any(grepl("normalised", assay_names)))
})

# Test adjust_abundance function

test_that("adjust_abundance adds adjusted assay correctly", {
  # Add a batch column to se_mini
  se_mini2 <- se_mini
  # Create a batch variable with two groups
  colData(se_mini2)$batch <- rep(1:2, length.out = ncol(se_mini2))
  # Run identify_abundant and adjust_abundance
  res <- se_mini2 |> identify_abundant() |> adjust_abundance(
    .factor_unwanted = batch,
    .factor_of_interest = condition,
    method = "combat_seq"
  )
  # Check that an adjusted assay is present
  assay_names <- names(SummarizedExperiment::assays(res))
  expect_true(any(grepl("_adjusted$", assay_names)))
})

# Test fill_missing_abundance function
test_that("fill_missing_abundance works correctly", {
  # This function doesn't exist for SummarizedExperiment
  # Skip this test for now
  skip("fill_missing_abundance not implemented for SummarizedExperiment")
})

# Test impute_missing_abundance function
test_that("impute_missing_abundance works correctly", {
  res <- se_mini |> impute_missing_abundance(.formula = ~ condition)
  
  expect_no_error(res)
}) 

# Test scale_abundance with custom suffix

test_that("scale_abundance uses custom suffix correctly", {
  res <- se_mini |> identify_abundant() |> scale_abundance(suffix = "_custom")
  expect_true("counts_custom" %in% names(SummarizedExperiment::assays(res)))
  expect_false("counts_scaled" %in% names(SummarizedExperiment::assays(res)))
})

test_that("scale_abundance default suffix still works", {
  res <- se_mini |> identify_abundant() |> scale_abundance()
  expect_true("counts_scaled" %in% names(SummarizedExperiment::assays(res)))
}) 

# Test adjust_abundance on a custom assay

test_that("adjust_abundance on custom assay creates correct adjusted assay name", {
  se_mini2 <- se_mini
  # Create a batch variable with two groups
  colData(se_mini2)$batch <- rep(1:2, length.out = ncol(se_mini2))
  # Create a custom assay
  SummarizedExperiment::assays(se_mini2)[["my_custom_assay"]] <- SummarizedExperiment::assay(se_mini2, "counts") + 1
  # Run identify_abundant and adjust_abundance on the custom assay
  res <- se_mini2 |> identify_abundant(abundance = "my_custom_assay") |> adjust_abundance(
    abundance = "my_custom_assay",
    .factor_unwanted = batch,
    .factor_of_interest = condition,
    method = "combat_seq"
  )
  # Check that the adjusted assay name is correct
  assay_names <- names(SummarizedExperiment::assays(res))
  expect_true("my_custom_assay_adjusted" %in% assay_names)
}) 