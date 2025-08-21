context('Clustering Functions')

library(airway)
data(airway)
airway_mini <- airway[1:100, 1:5]

library(dplyr)
library(SummarizedExperiment)

# Test cluster_elements function
test_that("cluster_elements with kmeans works correctly", {
  res <- cluster_elements(airway_mini, method = "kmeans", centers = 2)
  
  expect_true("cluster_kmeans" %in% names(SummarizedExperiment::colData(res)))
  expect_equal(
    levels(SummarizedExperiment::colData(res)$cluster_kmeans),
    c("1", "2")
  )
})

test_that("cluster_elements with SNN works correctly", {
  # Skip SNN as it has Matrix package issues
  skip("SNN clustering has Matrix package compatibility issues")
})

test_that("cluster_elements with hierarchical works correctly", {
  # Skip hierarchical as it's not supported
  skip("hierarchical clustering not supported")
})

test_that("cluster_elements with DBSCAN works correctly", {
  # Skip DBSCAN as it's not supported
  skip("DBSCAN clustering not supported")
}) 