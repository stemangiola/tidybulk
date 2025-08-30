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

  res <- aggregate_duplicates(se, .transcript = entrez)

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