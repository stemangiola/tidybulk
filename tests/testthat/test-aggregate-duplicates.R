# Test that aggregate_duplicates() preserves numeric columns and does not coerce them to character

test_that("aggregate_duplicates preserves numeric columns", {
  library(SummarizedExperiment)
  library(tidybulk)

  # minimal SE with duplicate transcript ids and mixed column types
  test_rowdata <- data.frame(
    entrez = c("A", "A", "B"),
    gene_name = c("geneA", "geneA", "geneB"),
    logFC = c(0.1, 0.2, -0.3),                 # numeric
    pvalue = c(0.05, 0.04, 0.9),               # numeric
    stringsAsFactors = FALSE
  )

  se <- SummarizedExperiment(
    assays = list(counts = matrix(1:6, nrow = 3, ncol = 2)),
    rowData = test_rowdata
  )

  res <- aggregate_duplicates(se, feature = "entrez")

  # Identify numeric columns in original rowData
  numeric_cols <- names(Filter(is.numeric, test_rowdata))

  # After aggregation, those columns should still be numeric
  expect_true(all(sapply(as.data.frame(rowData(res))[numeric_cols], is.numeric)))

  # Character columns should remain character
  char_cols <- setdiff(colnames(test_rowdata), numeric_cols)
  expect_true(all(sapply(as.data.frame(rowData(res))[char_cols], is.character)))

  # merged_transcripts column added should be integer
  expect_true(is.integer(rowData(res)$merged_transcripts))
}) 

test_that("aggregate_duplicates keeps integer type when keep_integer=TRUE", {
  library(SummarizedExperiment)
  library(tidybulk)

  # Create test data with integer counts and duplicates
  test_rowdata <- data.frame(
    entrez = c("A", "A", "B", "B", "C"),
    gene_name = c("geneA", "geneA", "geneB", "geneB", "geneC"),
    stringsAsFactors = FALSE
  )

  # Create integer matrix with duplicates
  se <- SummarizedExperiment(
    assays = list(counts = matrix(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L), nrow = 5, ncol = 2)),
    rowData = test_rowdata
  )

  # Test that original data is integer
  expect_true(is.integer(assay(se)))

  # Test with keep_integer=TRUE (default)
  res_true <- aggregate_duplicates(se, feature = "entrez", aggregation_function = mean, keep_integer = TRUE)
  expect_true(is.integer(assay(res_true)))
  
  # Test with keep_integer=FALSE
  res_false <- aggregate_duplicates(se, feature = "entrez", aggregation_function = mean, keep_integer = FALSE)
  expect_false(is.integer(assay(res_false)))
  expect_true(is.numeric(assay(res_false)))

  # Test with sum aggregation (should preserve integers naturally)
  res_sum <- aggregate_duplicates(se, feature = "entrez", aggregation_function = sum, keep_integer = TRUE)
  expect_true(is.integer(assay(res_sum)))

  # Test with sum aggregation and keep_integer=FALSE
  res_sum_false <- aggregate_duplicates(se, feature = "entrez", aggregation_function = sum, keep_integer = FALSE)
  expect_true(is.integer(assay(res_sum_false))) # sum of integers is still integer

  # Verify that the aggregation actually happened (fewer rows)
  expect_true(nrow(res_true) < nrow(se))
  expect_true(nrow(res_false) < nrow(se))
})

test_that("aggregate_duplicates keep_integer works with airway dataset", {
  library(SummarizedExperiment)
  library(tidybulk)
  library(airway)
  
  data('airway', package = 'airway')
  SummarizedExperiment::rowData(airway)$gene_name = rownames(airway)
  
  # Create duplicates by duplicating some genes
  rowData(airway)$gene_name[1:10] = 'DUPLICATE_GENE'
  
  # Test that original data is integer
  expect_true(is.integer(assay(airway)))
  
  # Test with keep_integer=TRUE
  result_true <- aggregate_duplicates(airway, feature = "gene_name", aggregation_function = mean, keep_integer = TRUE)
  expect_true(is.integer(assay(result_true)))
  
  # Test with keep_integer=FALSE
  result_false <- aggregate_duplicates(airway, feature = "gene_name", aggregation_function = mean, keep_integer = FALSE)
  expect_false(is.integer(assay(result_false)))
  expect_true(is.numeric(assay(result_false)))
  
  # Verify aggregation occurred
  expect_true(nrow(result_true) < nrow(airway))
  expect_true(nrow(result_false) < nrow(airway))
})

test_that("aggregate_duplicates shows deprecation warning for .transcript parameter", {
  library(SummarizedExperiment)
  library(tidybulk)
  
  # Create test data
  test_rowdata <- data.frame(
    entrez = c("A", "A", "B"),
    gene_name = c("geneA", "geneA", "geneB"),
    stringsAsFactors = FALSE
  )
  
  se <- SummarizedExperiment(
    assays = list(counts = matrix(1:6, nrow = 3, ncol = 2)),
    rowData = test_rowdata
  )
  
  # Test that using .transcript shows deprecation warning
  # We need to use the correct syntax for .transcript parameter
  expect_warning(
    aggregate_duplicates(se, .transcript = as.symbol("entrez")),
    "deprecated"
  )
  
  # Test that using feature parameter does not show warning
  expect_no_warning(
    aggregate_duplicates(se, feature = "entrez")
  )
  
  # Test that both approaches produce the same result
  suppressWarnings({
    result_old <- aggregate_duplicates(se, .transcript = as.symbol("entrez"))
  })
  result_new <- aggregate_duplicates(se, feature = "entrez")
  
  expect_identical(result_old, result_new)
})