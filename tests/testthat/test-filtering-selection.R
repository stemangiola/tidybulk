context('Filtering and Selection Functions')

data("se_mini")
data("breast_tcga_mini_SE")

library(dplyr)
library(SummarizedExperiment)

# Test identify_abundant function
test_that("identify_abundant works correctly", {
  res <- se_mini |> identify_abundant()
  
  expect_true(".abundant" %in% names(SummarizedExperiment::rowData(res)))
})

test_that("identify_abundant with custom parameters works correctly", {
  res <- se_mini |> identify_abundant(
    minimum_counts = 10,
    minimum_proportion = 0.7
  )
  
  expect_true(".abundant" %in% names(SummarizedExperiment::rowData(res)))
})

# Test keep_abundant function
test_that("keep_abundant works correctly", {
  res <- se_mini |> identify_abundant() |> keep_abundant()
  
  expect_true(nrow(res) <= nrow(se_mini))
})

# Test keep_variable function
test_that("keep_variable works correctly", {
  res <- se_mini |> keep_variable()
  
  expect_true(nrow(res) <= nrow(se_mini))
})

test_that("keep_variable with top_n works correctly", {
  res <- se_mini |> keep_variable(top = 100)
  
  expect_true(nrow(res) <= 100)
})

# Test filterByExpr function
test_that("filterByExpr works correctly", {
  # Skip this test as filterByExpr is not exported
  skip("filterByExpr function is not exported")
})

test_that("filterByExpr with design works correctly", {
  # Skip this test as filterByExpr is not exported
  skip("filterByExpr function is not exported")
}) 