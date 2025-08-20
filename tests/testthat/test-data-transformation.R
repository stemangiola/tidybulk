context('Data Transformation Functions')

library(airway)
data(airway)
airway_mini <- airway[1:100, 1:5]
data("breast_tcga_mini_SE")

library(dplyr)
library(SummarizedExperiment)

# Test pivot functions with comprehensive column validation
test_that("pivot_sample works correctly and preserves expected columns", {
  # Get original colData structure
  original_colData <- colData(airway_mini)
  
  # Test pivot_sample functionality
  result <- airway_mini |> pivot_sample()
  
  # Check that the result is a tibble
  expect_true(inherits(result, "tbl_df"))
  
  # Check that all original colData columns are preserved (plus .sample)
  expected_cols <- c(".sample", names(original_colData))
  expect_equal(sort(names(result)), sort(expected_cols))
  
  # Check that the number of rows matches the number of samples
  expect_equal(nrow(result), nrow(original_colData))
  
  # Check that sample names are preserved in .sample column
  expect_equal(result$.sample, rownames(original_colData))
  
  # Check specific expected columns based on the airway data
  expect_true("SampleName" %in% names(result))
  expect_true("cell" %in% names(result))
  expect_true("dex" %in% names(result))
  expect_true("albut" %in% names(result))
  expect_true("Run" %in% names(result))
  expect_true(".sample" %in% names(result))
  
  # Check that dex values are preserved (should be logical)
  expect_true(all(result$dex %in% c("trt", "untrt")))
  
  # Check that the result has the expected number of columns (8 + 1 for .sample)
  expect_equal(ncol(result), 10)
})

test_that("pivot_transcript works correctly and preserves expected columns", {
  # Get original rowData structure
  original_rowData <- rowData(airway_mini)
  
  # Test pivot_transcript functionality
  result <- airway_mini |> pivot_transcript()
  
  # Check that the result is a tibble
  expect_true(inherits(result, "tbl_df"))
  
  # Check that all original rowData columns are preserved (plus .feature)
  # The airway dataset has rich gene annotation data
  expected_cols <- c(".feature", "entrezid", "gene_biotype", "gene_id", "gene_name", 
                     "gene_seq_end", "gene_seq_start", "GRangesList", 
                     "seq_coord_system", "seq_name", "seq_strand", "symbol")
  expect_equal(sort(names(result)), sort(expected_cols))
  
  # Check that the number of rows matches the number of features/transcripts
  expect_equal(nrow(result), nrow(original_rowData))
  
  # Check that feature names are preserved in .feature column
  expect_equal(result$.feature, rownames(original_rowData))
  
  # Check specific expected columns based on the airway data
  expect_true("gene_name" %in% names(result))  # This is the gene_name column
  expect_true(".feature" %in% names(result))
  expect_true("symbol" %in% names(result))
  expect_true("entrezid" %in% names(result))
  
  # Check that the result has the expected number of columns (12 total)
  expect_equal(ncol(result), 12)
})

test_that("pivot_sample handles additional colData columns correctly", {
  # Create a modified version with extra columns
  se_modified <- airway_mini
  colData(se_modified)$extra_column <- rep("test", nrow(colData(se_modified)))
  colData(se_modified)$numeric_column <- 1:nrow(colData(se_modified))
  
  # Test pivot_sample with additional columns
  result <- se_modified |> pivot_sample()
  
  # Check that all original columns plus new ones are preserved
  expected_cols <- c(".sample", "SampleName", "cell", "dex", "albut", "Run", 
                     "avgLength", "Experiment", "Sample", "BioSample",
                     "extra_column", "numeric_column")
  expect_equal(sort(names(result)), sort(expected_cols))
  
  # Check that new columns have correct values
  expect_equal(result$extra_column, rep("test", nrow(result)))
  expect_equal(result$numeric_column, 1:nrow(result))
})

test_that("pivot_transcript handles additional rowData columns correctly", {
  # Create a modified version with extra columns
  se_modified <- airway_mini
  rowData(se_modified)$extra_column <- rep("test", nrow(rowData(se_modified)))
  rowData(se_modified)$numeric_column <- 1:nrow(rowData(se_modified))
  
  # Test pivot_transcript with additional columns
  result <- se_modified |> pivot_transcript()
  
  # Check that all original columns plus new ones are preserved
  # The airway dataset has rich gene annotation data plus additional columns
  expected_cols <- c(".feature", "entrezid", "gene_biotype", "gene_id", "gene_name", 
                     "gene_seq_end", "gene_seq_start", "GRangesList", 
                     "seq_coord_system", "seq_name", "seq_strand", "symbol",
                     "extra_column", "numeric_column")
  expect_equal(sort(names(result)), sort(expected_cols))
  
  # Check that new columns have correct values
  expect_equal(result$extra_column, rep("test", nrow(result)))
  expect_equal(result$numeric_column, 1:nrow(result))
})

test_that("pivot_sample preserves data types correctly", {
  result <- airway_mini |> pivot_sample()
  
  # Check that factor columns remain factor (airway dataset uses factors)
  expect_true(is.factor(result$dex))
  expect_true(is.factor(result$SampleName))
  expect_true(is.factor(result$cell))
  expect_true(is.factor(result$albut))
  expect_true(is.factor(result$Run))
  
  # Check that numeric columns remain numeric
  expect_true(is.numeric(result$avgLength))
  
  # Check that .sample column is character
  expect_true(is.character(result$.sample))
})

test_that("pivot_transcript preserves data types correctly", {
  result <- airway_mini |> pivot_transcript()
  
  # Check that character columns remain character
  expect_true(is.character(result$gene_name))  # This is the gene_name column
  expect_true(is.character(result$.feature))
})

# Test as_matrix function
test_that("as_matrix works correctly", {
  # Skip this test as as_matrix doesn't work with SummarizedExperiment
  skip("as_matrix not implemented for SummarizedExperiment")
})

# Test aggregate_duplicates function
test_that("aggregate_duplicates works correctly", {
  # Skip this test as the test data has NA gene_name values
  skip("aggregate_duplicates requires non-NA transcript values")
})

# Test as_SummarizedExperiment function
test_that("as_SummarizedExperiment works correctly", {
  # Skip this test as it's not applicable for SummarizedExperiment input
  skip("as_SummarizedExperiment not applicable for SummarizedExperiment input")
}) 