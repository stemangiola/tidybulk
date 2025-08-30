# Tidybulk Test Suite Organization

This directory contains unit tests organized by functionality groups for the tidybulk package.

## Test Organization

### 1. Data Transformation Functions (`test-data-transformation.R`)
- **pivot_sample()** - Pivot data by sample
- **pivot_transcript()** - Pivot data by transcript
- **as_matrix()** - Convert to matrix format
- **aggregate_duplicates()** - Aggregate duplicate entries
- **as_SummarizedExperiment()** - Convert to SummarizedExperiment

### 2. Abundance Normalization Functions (`test-abundance-normalization.R`)
- **scale_abundance()** - Scale abundance data
- **quantile_normalise_abundance()** - Quantile normalization
- **adjust_abundance()** - Adjust abundance for unwanted variation
- **fill_missing_abundance()** - Fill missing abundance values
- **impute_missing_abundance()** - Impute missing abundance values

### 3. Filtering and Selection Functions (`test-filtering-selection.R`)
- **identify_abundant()** - Identify abundant transcripts
- **keep_abundant()** - Keep abundant transcripts
- **keep_variable()** - Keep variable transcripts
- **filterByExpr()** - Filter by expression

### 4. Dimensionality Reduction Functions (`test-dimensionality-reduction.R`)
- **reduce_dimensions()** - Reduce dimensions with PCA/MDS/tSNE/UMAP
- **rotate_dimensions()** - Rotate dimensions
- **remove_redundancy()** - Remove redundant features

### 5. Clustering Functions (`test-clustering.R`)
- **cluster_elements()** - Cluster elements with various methods
- **kmeans clustering** - K-means clustering
- **SNN clustering** - Shared nearest neighbor clustering
- **hierarchical clustering** - Hierarchical clustering
- **DBSCAN clustering** - Density-based clustering

### 6. Differential Analysis Functions (`test-differential-analysis.R`)
- **test_differential_abundance()** - Test differential abundance with various methods
- **test_differential_cellularity()** - Test differential cellularity
- **glmmSeq()** - Generalized linear mixed models

### 7. Cellularity Analysis Functions (`test-cellularity-analysis.R`)
- **deconvolve_cellularity()** - Deconvolve cellularity with various methods
- **cibersort()** - CIBERSORT analysis

### 8. Gene Enrichment Functions (`test-gene-enrichment.R`)
- **test_gene_enrichment()** - Test gene enrichment
- **test_gene_overrepresentation()** - Test gene overrepresentation
- **test_gene_rank()** - Test gene rank

### 9. Utility Functions (`test-utility-functions.R`)
- **describe_transcript()** - Describe transcript characteristics
- **get_bibliography()** - Get bibliography
- **resolve_complete_confounders_of_non_interest()** - Resolve confounders

### 10. Validation and Utility Functions (`test-validation-utilities.R`)
- **Validation functions** - Check data integrity
- **Utility functions** - Log transformations and other utilities

## Test Setup

### Common Setup (`test-setup.R`)
- Common test data loading
- Helper functions for testing
- Shared utilities across test files

## Running Tests

To run all tests:
```r
library(testthat)
library(tidybulk)
test_dir('tests/testthat/', reporter = 'summary')
```

To run specific test files:
```r
test_file('tests/testthat/test-data-transformation.R')
test_file('tests/testthat/test-abundance-normalization.R')
```

## Test Coverage

The organized test suite provides comprehensive coverage of:
- Core tidybulk functionality
- Data transformation operations
- Abundance normalization methods
- Filtering and selection operations
- Dimensionality reduction techniques
- Clustering methods
- Differential analysis workflows
- Cellularity analysis
- Gene enrichment analysis
- Utility and validation functions

## Benefits of Organization

1. **Clear Functional Grouping** - Tests are organized by related functionality
2. **Easy Maintenance** - Related tests are in the same file
3. **Better Debugging** - Failures are grouped by functionality
4. **Improved Coverage** - Systematic testing of all major functions
5. **Documentation** - Test files serve as usage examples 