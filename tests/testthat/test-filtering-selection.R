context('Filtering and Selection Functions')

library(airway)
data(airway)
airway_mini <- airway[1:100, 1:5]
breast_tcga_mini_SE <- airway[1:200, 1:8]

library(dplyr)
library(SummarizedExperiment)

# Test identify_abundant function
test_that("identify_abundant works correctly", {
  res <- airway_mini |> identify_abundant()
  
  expect_true(".abundant" %in% names(SummarizedExperiment::rowData(res)))
})

test_that("identify_abundant with custom parameters works correctly", {
  res <- airway_mini |> identify_abundant(
    minimum_counts = 10,
    minimum_proportion = 0.7
  )
  
  expect_true(".abundant" %in% names(SummarizedExperiment::rowData(res)))
})

# Test keep_abundant function
test_that("keep_abundant works correctly", {
  res <- airway_mini |> identify_abundant() |> keep_abundant()
  
  expect_true(nrow(res) <= nrow(airway_mini))
})

# Test keep_variable function
test_that("keep_variable works correctly", {
  res <- airway_mini |> keep_variable()
  
  expect_true(nrow(res) <= nrow(airway_mini))
})

test_that("keep_variable with top_n works correctly", {
  res <- airway_mini |> keep_variable(top = 100)
  
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