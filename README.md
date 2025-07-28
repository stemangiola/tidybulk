tidybulk: An R tidy framework for modular transcriptomic data analysis
================
Stefano Mangiola
2025-07-28

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidybulk/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidybulk/actions/)
[![Bioconductor
status](https://bioconductor.org/shields/build/release/bioc/tidybulk.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/tidybulk/)
<!-- badges: end -->

# <img src="inst/logo.svg" height="139px" width="120px"/>

**tidybulk** is a powerful R package designed for modular transcriptomic
data analysis that brings transcriptomics to the tidyverse.

## Why tidybulk?

Tidybulk provides a unified interface for comprehensive transcriptomic
data analysis with seamless integration of SummarizedExperiment objects
and tidyverse principles. It streamlines the entire workflow from raw
data to biological insights.

# <img src="inst/new_SE_usage-01.png" width="100%"/>

## Functions/utilities available

### Abundance Normalization Functions

| Function                         | Description                             |
|----------------------------------|-----------------------------------------|
| `scale_abundance()`              | Scale abundance data                    |
| `quantile_normalise_abundance()` | Quantile normalization                  |
| `adjust_abundance()`             | Adjust abundance for unwanted variation |
| `fill_missing_abundance()`       | Fill missing abundance values           |
| `impute_missing_abundance()`     | Impute missing abundance values         |

### Filtering and Selection Functions

| Function              | Description                                         |
|-----------------------|-----------------------------------------------------|
| `identify_abundant()` | Identify abundant transcripts without removing them |
| `keep_abundant()`     | Keep abundant transcripts                           |
| `keep_variable()`     | Keep variable transcripts                           |
| `filterByExpr()`      | Filter by expression                                |

### Dimensionality Reduction Functions

| Function              | Description                              |
|-----------------------|------------------------------------------|
| `reduce_dimensions()` | Reduce dimensions with PCA/MDS/tSNE/UMAP |
| `rotate_dimensions()` | Rotate dimensions                        |
| `remove_redundancy()` | Remove redundant features                |

### Clustering Functions

| Function                  | Description                           |
|---------------------------|---------------------------------------|
| `cluster_elements()`      | Cluster elements with various methods |
| `kmeans clustering`       | K-means clustering                    |
| `SNN clustering`          | Shared nearest neighbor clustering    |
| `hierarchical clustering` | Hierarchical clustering               |
| `DBSCAN clustering`       | Density-based clustering              |

### Differential Analysis Functions

| Function | Description |
|----|----|
| `test_differential_expression()` | Test differential expression with various methods |

### Cellularity Analysis Functions

| Function | Description |
|----|----|
| `deconvolve_cellularity()` | Deconvolve cellularity with various methods |
| `test_stratification_cellularity()` | Test stratification cellularity (DEPRECATED) |
| `cibersort()` | CIBERSORT analysis |

### Gene Enrichment Functions

| Function                         | Description                  |
|----------------------------------|------------------------------|
| `test_gene_enrichment()`         | Test gene enrichment         |
| `test_gene_overrepresentation()` | Test gene overrepresentation |
| `test_gene_rank()`               | Test gene rank               |

### Utility Functions

| Function | Description |
|----|----|
| `describe_transcript()` | Describe transcript characteristics |
| `get_bibliography()` | Get bibliography |
| `resolve_complete_confounders_of_non_interest()` | Resolve confounders |

### Validation and Utility Functions

| Function                      | Description                       |
|-------------------------------|-----------------------------------|
| `check_if_counts_is_na()`     | Check if counts contain NA values |
| `check_if_duplicated_genes()` | Check for duplicated genes        |
| `check_if_wrong_input()`      | Validate input data               |
| `log10_reverse_trans()`       | Log10 reverse transformation      |
| `logit_trans()`               | Logit transformation              |

All functions are directly compatible with `SummarizedExperiment`
objects and follow tidyverse principles for seamless integration with
the tidyverse ecosystem.

### Scientific Citation

Mangiola, Stefano, Ramyar Molania, Ruining Dong, Maria A. Doyle, and
Anthony T. Papenfuss. 2021. “Tidybulk: An R tidy framework for modular
transcriptomic data analysis.” Genome Biology 22 (42).
<https://doi.org/10.1186/s13059-020-02233-7>

[Genome Biology - tidybulk: an R tidy framework for modular
transcriptomic data
analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02233-7)

# Installation Guide

**Bioconductor**

``` r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("tidybulk")
```

**Github**

``` r
devtools::install_github("stemangiola/tidybulk")
```

# Comprehensive Example Pipeline

This vignette demonstrates a complete transcriptomic analysis workflow
using tidybulk, with special emphasis on differential expression
analysis.

## Data Overview

We will use a `SummarizedExperiment` object containing RNA-seq data:

``` r
se_mini
```

    ## class: SummarizedExperiment 
    ## dim: 527 5 
    ## metadata(0):
    ## assays(1): counts
    ## rownames(527): ABCB4 ABCB9 ... ZNF324 ZNF442
    ## rowData names(1): entrez
    ## colnames(5): SRR1740034 SRR1740035 SRR1740043 SRR1740058 SRR1740067
    ## colData names(5): Cell.type time condition days dead

Loading `tidySummarizedExperiment` automatically abstracts this object
as a `tibble`, making it compatible with tidyverse tools while
maintaining its `SummarizedExperiment` nature:

``` r
class(se_mini)
```

    ## [1] "SummarizedExperiment"
    ## attr(,"package")
    ## [1] "SummarizedExperiment"

### Prepare Data for Analysis

Before analysis, we need to ensure our variables are in the correct
format:

``` r
# Convert condition to factor for proper differential expression analysis
colData(se_mini)$condition = as.factor(colData(se_mini)$condition)
```

### Visualize Raw Counts

Visualize the distribution of raw counts before any filtering:

``` r
ggplot(as_tibble(se_mini), aes(counts + 1, group = .sample, color = `Cell.type`)) +
  geom_density() +
  scale_x_log10() +
  my_theme +
  labs(title = "Raw counts by cell type (before any filtering)")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/plot-raw-counts-1.png)<!-- -->

## Step 1: Data Preprocessing

### Aggregate Duplicated Transcripts (optional)

Aggregate duplicated transcripts (e.g., isoforms, ensembl IDs):

> Transcript aggregation is a standard bioinformatics approach for
> gene-level summarization.

``` r
# Add gene names to rowData
rowData(se_mini)$gene_name = rownames(se_mini)

# Aggregate duplicates
se_mini = se_mini |> aggregate_duplicates(.transcript = gene_name, aggregation_function = mean)
```

    ## tidybulk says: your object does not have duplicates along the gene_name column. The input dataset is returned.

### Abundance Filtering: tidybulk approaches only

Abundance filtering can be performed using tidybulk’s built-in methods
([Robinson, McCarthy, and Smyth 2010](#ref-robinson2010edger); [Chen,
Lun, and Smyth 2016](#ref-chen2016edgeR)).

#### 1. tidybulk: Default, formula_design, and CPM threshold

``` r
# Default (simple filtering)
se_abundant_default = se_mini |> keep_abundant()
```

    ## Warning in filterByExpr.DGEList(y, design = design, group = group, lib.size =
    ## lib.size, : All samples appear to belong to the same group.

``` r
# With factor_of_interest (recommended for complex designs)
se_abundant_formula = se_mini |> keep_abundant(minimum_counts = 10, minimum_proportion = 0.5, factor_of_interest = condition)
```

    ## Warning: The `factor_of_interest` argument of `keep_abundant()` is deprecated as of
    ## tidybulk 2.0.0.
    ## ℹ Please use the `formula_design` argument instead.
    ## ℹ The argument 'factor_of_interest' is deprecated and will be removed in a
    ##   future release. Please use the 'design' or 'formula_design' argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The `factor_of_interest` argument of `identify_abundant()` is deprecated as of
    ## tidybulk 2.0.0.
    ## ℹ Please use the `formula_design` argument instead.
    ## ℹ The argument 'factor_of_interest' is deprecated and will be removed in a
    ##   future release. Please use the 'design' or 'formula_design' argument instead.
    ## ℹ The deprecated feature was likely used in the tidybulk package.
    ##   Please report the issue at <https://github.com/stemangiola/tidybulk/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
# With CPM threshold (using design parameter)
se_abundant_cpm = se_mini |> keep_abundant(minimum_counts = 10, minimum_proportion = 0.5)
```

    ## Warning in filterByExpr.DGEList(y, design = design, group = group, lib.size =
    ## lib.size, : All samples appear to belong to the same group.

#### 2. Summary statistics and density plots

``` r
# Example: summary for default tidybulk filtering
# Before filtering
se_mini |> as_tibble() |> summarise(
  n_features = n_distinct(.feature),
  min_count = min(counts),
  median_count = median(counts),
  max_count = max(counts)
)
```

    ## # A tibble: 1 × 4
    ##   n_features min_count median_count max_count
    ##        <int>     <dbl>        <dbl>     <dbl>
    ## 1        527         0           26    134561

``` r
# After filtering
se_abundant_default |> as_tibble() |> summarise(
  n_features = n_distinct(.feature),
  min_count = min(counts),
  median_count = median(counts),
  max_count = max(counts)
)
```

    ## # A tibble: 1 × 4
    ##   n_features min_count median_count max_count
    ##        <int>     <dbl>        <dbl>     <dbl>
    ## 1        182         7         370.    134561

``` r
se_abundant_formula |> as_tibble() |> summarise(
  n_features = n_distinct(.feature),
  min_count = min(counts),
  median_count = median(counts),
  max_count = max(counts)
)
```

    ## # A tibble: 1 × 4
    ##   n_features min_count median_count max_count
    ##        <int>     <dbl>        <dbl>     <dbl>
    ## 1        394         0         120.    134561

``` r
se_abundant_cpm |> as_tibble() |> summarise(
  n_features = n_distinct(.feature),
  min_count = min(counts),
  median_count = median(counts),
  max_count = max(counts)
)
```

    ## # A tibble: 1 × 4
    ##   n_features min_count median_count max_count
    ##        <int>     <dbl>        <dbl>     <dbl>
    ## 1        182         7         370.    134561

``` r
# Merge all methods into a single tibble
se_abundant_all = 
  bind_rows(
    se_mini |> assay() |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "no filter"),
    se_abundant_default |> assay() |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "default"),
    se_abundant_formula |> assay() |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "formula"),
    se_abundant_cpm |> assay() |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "cpm")
  )

# Density plot across methods
se_abundant_all |> 
  as_tibble() |> 
  ggplot(aes(counts + 1, group = .sample, color = method)) +
    geom_density() +
    scale_x_log10() +
    facet_wrap(~method) +
    my_theme +
    labs(title = "Counts after abundance filtering (tidybulk default)")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/filtering-density-plot-comparison-1.png)<!-- -->

Update the `se_mini` object with the filtered data:

``` r
se_mini = se_abundant_formula
```

> **Tip:** Use `formula_design` for complex designs, and use the CPM
> threshold for library-size-aware filtering.

### Remove Redundant Transcripts

Redundancy removal is a standard approach for reducing highly correlated
features.

``` r
se_mini_non_redundant = 
  se_mini |> 
  remove_redundancy(method = "correlation", top = 100) 
```

    ## Getting the 100 most variable genes

``` r
  # Make  

se_mini |> as_tibble() |> summarise(
  n_features = n_distinct(.feature),
  min_count = min(counts),
  median_count = median(counts),
  max_count = max(counts)
)
```

    ## # A tibble: 1 × 4
    ##   n_features min_count median_count max_count
    ##        <int>     <dbl>        <dbl>     <dbl>
    ## 1        394         0         120.    134561

``` r
# Summary statistics
se_mini_non_redundant |> as_tibble() |> summarise(
  n_features = n_distinct(.feature),
  min_count = min(counts),
  median_count = median(counts),
  max_count = max(counts)
)
```

    ## # A tibble: 1 × 4
    ##   n_features min_count median_count max_count
    ##        <int>     <dbl>        <dbl>     <dbl>
    ## 1        394         0         120.    134561

``` r
# Plot before and after
# Merge before and after into a single tibble
se_mini_all = bind_rows(
  se_mini |> assay() |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |>  mutate(method = "before"),
  se_mini_non_redundant |> assay() |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |>  mutate(method = "after")
)

# Density plot
ggplot(as_tibble(se_mini_all), aes(counts + 1, group = .sample, color = method)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~method) +
  my_theme +
  labs(title = "Counts after removing redundant transcripts")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/preprocessing-remove-redundancy-1.png)<!-- -->

### Filter Variable Transcripts

Keep only the most variable transcripts for downstream analysis.

Variance-based feature selection using edgeR methodology ([Robinson,
McCarthy, and Smyth 2010](#ref-robinson2010edger)) is used for selecting
informative features.

``` r
se_mini_variable = se_mini |> keep_variable()
```

    ## Getting the 394 most variable genes

### Visualize After Variable Filtering Variable Transcripts (optional)

``` r
# Before filtering
se_mini |> as_tibble() |> summarise(
  n_features = n_distinct(.feature),
  min_count = min(counts),
  median_count = median(counts),
  max_count = max(counts)
)
```

    ## # A tibble: 1 × 4
    ##   n_features min_count median_count max_count
    ##        <int>     <dbl>        <dbl>     <dbl>
    ## 1        394         0         120.    134561

``` r
# After filtering
se_mini_variable |> as_tibble() |> summarise(
  n_features = n_distinct(.feature),
  min_count = min(counts),
  median_count = median(counts),
  max_count = max(counts)
)
```

    ## # A tibble: 1 × 4
    ##   n_features min_count median_count max_count
    ##        <int>     <dbl>        <dbl>     <dbl>
    ## 1        394         0         120.    134561

``` r
# Density plot
# Merge before and after into a single tibble
se_mini_all = bind_rows(
  se_mini |> assay() |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "before"),
  se_mini_variable |> assay() |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "after")
)

# Density plot
ggplot(as_tibble(se_mini_all), aes(counts + 1, group = .sample, color = method)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~method) +
  my_theme +
  labs(title = "Counts after variable filtering")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/filtering-variable-summary-and-plot-1.png)<!-- -->

### Scale Abundance

Scale for sequencing depth using TMM ([Robinson, McCarthy, and Smyth
2010](#ref-robinson2010edger)), upper quartile ([Bullard et al.
2010](#ref-bullard2010uq)), and RLE ([Anders and Huber
2010](#ref-anders2010rle)) normalization.

``` r
se_mini = 
se_mini |> 
    scale_abundance(method = "TMM", suffix = "_tmm") |>
    scale_abundance(method = "upperquartile", suffix = "_upperquartile") |>
    scale_abundance(method = "RLE", suffix = "_RLE")
```

    ## tidybulk says: the sample with largest library size SRR1740035 was chosen as reference for scaling
    ## tidybulk says: the sample with largest library size SRR1740035 was chosen as reference for scaling
    ## tidybulk says: the sample with largest library size SRR1740035 was chosen as reference for scaling

### Visualize After Scaling

``` r
# Before scaling
se_mini |> assay("counts") |> as.matrix() |> rowMeans() |> summary()
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     6.6   115.8   540.5  1769.1  1608.3 48505.2

``` r
se_mini |> assay("counts_tmm") |> as.matrix() |> rowMeans() |> summary()
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ##     8.296   126.854   618.203  1964.138  1919.999 50581.378

``` r
se_mini |> assay("counts_upperquartile") |> as.matrix() |> rowMeans() |> summary()
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ##     8.501   144.420   660.466  2095.482  1971.146 51299.288

``` r
se_mini |> assay("counts_RLE") |> as.matrix() |> rowMeans() |> summary()
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ##     8.922   145.443   700.123  2190.865  2213.889 52881.398

``` r
# Merge all methods into a single tibble
se_mini_scaled_all = bind_rows(
  se_mini |> assay("counts") |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "no_scaling"),
  se_mini |> assay("counts_tmm") |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "TMM"),
  se_mini |> assay("counts_upperquartile") |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "upperquartile"),
  se_mini |> assay("counts_RLE") |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts") |> mutate(method = "RLE")
)


# Density plot
ggplot(as_tibble(se_mini_scaled_all), aes(counts + 1, group = .sample, color = method)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~method) +
  my_theme +
  labs(title = "Scaled counts by method (after scaling)")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/normalization-visualize-scaling-1.png)<!-- -->

## Step 2: Exploratory Data Analysis

### Remove Zero-Variance Features (required for PCA)

Variance filtering is a standard preprocessing step for dimensionality
reduction.

``` r
library(matrixStats)
# Remove features with zero variance across samples
se_mini = se_mini[rowVars(assay(se_mini)) > 0, ]
```

### Dimensionality Reduction

MDS ([Kruskal 1964](#ref-kruskal1964mds)) using limma::plotMDS
([**ritchie2015limma?**](#ref-ritchie2015limma)) and PCA ([Hotelling
1933](#ref-hotelling1933pca)) are used for dimensionality reduction.

``` r
se_mini = se_mini |>
  reduce_dimensions(method="MDS", .dims = 2)
```

    ## Warning in reduce_dimensions(se_mini, method = "MDS", .dims = 2): tidybulk
    ## says: the "top" argument 500 is higher than the number of features 394

    ## Getting the 394 most variable genes

    ## [1] "MDS result_df colnames: sample, 1, 2"

    ## tidybulk says: to access the raw results do `metadata(.)$tidybulk$MDS`

``` r
se_mini = se_mini |>
  reduce_dimensions(method="PCA", .dims = 2)
```

    ## Warning in reduce_dimensions(se_mini, method = "PCA", .dims = 2): tidybulk
    ## says: the "top" argument 500 is higher than the number of features 394

    ## Getting the 394 most variable genes

    ## Fraction of variance explained by the selected principal components

    ## # A tibble: 2 × 2
    ##   `Fraction of variance`    PC
    ##                    <dbl> <int>
    ## 1                 0.0878     1
    ## 2                 0.0466     2

    ## tidybulk says: to access the raw results do `metadata(.)$tidybulk$PCA`

### Visualize Dimensionality Reduction Results

``` r
# MDS plot
se_mini |>
    pivot_sample() |>
    ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell.type`)) +
  geom_point() +
    my_theme +
    labs(title = "MDS Analysis")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/eda-plot-dimensionality-reduction-1.png)<!-- -->

``` r
# PCA plot
    se_mini |>
    pivot_sample() |>
    ggplot(aes(x=`PC1`, y=`PC2`, color=`Cell.type`)) +
    geom_point() +
    my_theme +
    labs(title = "PCA Analysis")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/eda-plot-dimensionality-reduction-2.png)<!-- -->

### Clustering Analysis

K-means clustering ([MacQueen 1967](#ref-macqueen1967kmeans)) is used
for unsupervised grouping.

``` r
se_mini = se_mini |>
  cluster_elements(method="kmeans", centers = 2)

# Visualize clustering
    se_mini |>
    ggplot(aes(x=`Dim1`, y=`Dim2`, color=`cluster_kmeans`)) +
  geom_point() +
  my_theme +
  labs(title = "K-means Clustering")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/eda-clustering-analysis-1.png)<!-- -->

## Step 3: Differential Expression Analysis

This workflow uses

**edgeR quasi-likelihood** ([Robinson, McCarthy, and Smyth
2010](#ref-robinson2010edger); [Chen, Lun, and Smyth
2016](#ref-chen2016edgeR)) **edgeR robust likelihood ratio** ([Chen,
Lun, and Smyth 2016](#ref-chen2016edgeR)), **DESeq2** ([Love, Huber, and
Anders 2014](#ref-love2014deseq2)), **limma-voom** ([Law et al.
2014](#ref-law2014voom)), and **limma-voom with sample weights** ([Liu
et al. 2015](#ref-liu2015voomweights)) for differential expression
analysis.

### Basic Differential Expression

**Methods:**

- **edgeR quasi-likelihood:** Quasi-likelihood F-tests for differential
  expression

- **edgeR robust likelihood ratio:** Robust likelihood ratio tests

- **DESeq2:** Negative binomial distribution with dispersion estimation

- **limma-voom:** Linear modeling with empirical Bayes moderation

- **limma-voom with sample weights:** Enhanced voom with quality weights
  **References:**

- Robinson et al. (2010) edgeR: a Bioconductor package for differential
  expression analysis

- Chen et al. (2016) From reads to genes to pathways: differential
  expression analysis of RNA-Seq experiments using Rsubread and the
  edgeR quasi-likelihood pipeline

- Love et al. (2014) Moderated estimation of fold change and dispersion
  for RNA-seq data with DESeq2

- Law et al. (2014) voom: precision weights unlock linear model analysis
  tools for RNA-seq read counts

- Liu et al. (2015) Why weight? Modelling sample and observational level
  variability improves power in RNA-seq analyses

``` r
# Standard differential expression analysis
se_mini = se_mini |>

# Use QL method
    test_differential_expression(~ condition, method = "edgeR_quasi_likelihood", prefix = "ql__") |>
    
    # Use edger_robust_likelihood_ratio
    test_differential_expression(~ condition, method = "edger_robust_likelihood_ratio", prefix = "lr_robust__") |>
    
# Use DESeq2 method
    test_differential_expression(~ condition, method = "DESeq2", prefix = "deseq2__") |>
    
    # Use limma_voom
    test_differential_expression(~ condition, method = "limma_voom", prefix = "voom__") |>

# Use limma_voom_sample_weights
    test_differential_expression(~ condition, method = "limma_voom_sample_weights", prefix = "voom_weights__") 
```

    ## Warning: The `.abundance` argument of `test_differential_abundance()` is deprecated as
    ## of tidybulk 2.0.0.
    ## ℹ Please use the `abundance` argument instead.
    ## ℹ The deprecated feature was likely used in the tidybulk package.
    ##   Please report the issue at <https://github.com/stemangiola/tidybulk/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## =====================================
    ## tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance
    ## or adjust_abundance have been calculated. Therefore, it is essential to add covariates
    ## such as batch effects (if applicable) in the formula.
    ## =====================================
    ## tidybulk says: The design column names are "(Intercept), conditionTRUE"
    ## 
    ## tidybulk says: to access the DE object do `metadata(.)$tidybulk$edgeR_quasi_likelihood_object`
    ## tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$edgeR_quasi_likelihood_fit`
    ## tidybulk says: The design column names are "(Intercept), conditionTRUE"
    ## 
    ## tidybulk says: to access the DE object do `metadata(.)$tidybulk$edger_robust_likelihood_ratio_object`
    ## tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$edger_robust_likelihood_ratio_fit`
    ## converting counts to integer mode
    ## 
    ## estimating size factors
    ## 
    ## estimating dispersions
    ## 
    ## gene-wise dispersion estimates
    ## 
    ## mean-dispersion relationship
    ## 
    ## -- note: fitType='parametric', but the dispersion trend was not well captured by the
    ##    function: y = a/x + b, and a local regression fit was automatically substituted.
    ##    specify fitType='local' or 'mean' to avoid this message next time.
    ## 
    ## final dispersion estimates
    ## 
    ## fitting model and testing
    ## 
    ## tidybulk says: to access the DE object do `metadata(.)$tidybulk$DESeq2_object`
    ## tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$DESeq2_fit`
    ## tidybulk says: The design column names are "(Intercept), conditionTRUE"
    ## 
    ## tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
    ## tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
    ## tidybulk says: The design column names are "(Intercept), conditionTRUE"
    ## 
    ## tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_sample_weights_object`
    ## tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_sample_weights_fit`
    ## This message is displayed once per session.

### Quality Control of the Fit

It is important to check the quality of the fit. All methods produce a
fit object that can be used for quality control. The fit object produced
by each underlying method are stored in as attributes of the `se_mini`
object. We can use them for example to perform quality control of the
fit.

#### For edgeR

Plot the biological coefficient of variation (BCV) trend. This plot is
helpful to understant the dispersion of the data.

``` r
library(edgeR)
```

    ## Warning: package 'edgeR' was built under R version 4.5.1

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
metadata(se_mini)$tidybulk$edgeR_quasi_likelihood_object |>
  plotBCV()
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-edgeR-object-1.png)<!-- -->

Plot the log-fold change vs mean plot.

``` r
library(edgeR)

metadata(se_mini)$tidybulk$edgeR_quasi_likelihood_fit |>
  plotMD()
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-edgeR-fit-1.png)<!-- -->

#### For DESeq2

Plot the mean-variance trend.

``` r
library(DESeq2)

metadata(se_mini)$tidybulk$DESeq2_object |>
  plotDispEsts()
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-DESeq2-object-1.png)<!-- -->

Plot the log-fold change vs mean plot.

``` r
library(DESeq2)

metadata(se_mini)$tidybulk$DESeq2_object |>
  plotMA()
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-DESeq2-fit-1.png)<!-- -->

### Histograms of p-values across methods

Inspection of the raw p-value histogram provides a rapid check of
differential-expression results. When no gene is truly differentially
expressed, the p-values follow a uniform U(0,1) distribution across the
interval 0–1, so the histogram appears flat
[Source](https://bioconductor.org/help/course-materials/2014/useR2014/Workflows.html).
In a more realistic scenario where only a subset of genes changes, this
uniform background is still present but an obvious spike emerges close
to zero, created by the genuine signals.

Thanks to the modularity of the `tidybulk` workflow, that can multiplex
different methods, we can easily compare the p-values across methods.

``` r
    se_mini |>
  rowData() |> 
  as_tibble() |> 
  select(
    ql__PValue, 
    lr_robust__PValue, 
    voom__P.Value, 
    voom_weights__P.Value, 
    deseq2__pvalue
  ) |> 
  pivot_longer(everything(), names_to = "method", values_to = "pvalue") |>
  ggplot(aes(x = pvalue, fill = method)) +
  geom_histogram(binwidth = 0.01) +
  facet_wrap(~method) +
  my_theme +
  labs(title = "Histogram of p-values across methods")
```

    ## Warning: Removed 36 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-pvalue-histograms-1.png)<!-- -->

### Compare Results Across Methods

``` r
# Summay statistics
se_mini |> rowData() |> as_tibble() |> select(contains("ql|lr_robust|voom|voom_weights|deseq2")) |> select(contains("logFC")) |> 
summarise(across(everything(), list(min = min, median = median, max = max), na.rm = TRUE))
```

    ## Warning: There was 1 warning in `summarise()`.
    ## ℹ In argument: `across(...)`.
    ## Caused by warning:
    ## ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
    ## Supply arguments directly to `.fns` through an anonymous function instead.
    ## 
    ##   # Previously
    ##   across(a:b, mean, na.rm = TRUE)
    ## 
    ##   # Now
    ##   across(a:b, \(x) mean(x, na.rm = TRUE))

    ## # A tibble: 1 × 0

### Pairplot of pvalues across methods (GGpairs)

``` r
library(GGally)
```

    ## Warning: package 'GGally' was built under R version 4.5.1

``` r
se_mini |> 
  rowData() |> 
  as_tibble() |> 
  select(ql__PValue, lr_robust__PValue, voom__P.Value, voom_weights__P.Value, deseq2__pvalue) |> 
  ggpairs(columns = 1:5) +
  scale_x_continuous(trans = tidybulk::log10_reverse_trans()) +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) +
  my_theme +
  labs(title = "Pairplot of p-values across methods")
```

    ## Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
    ## Removed 36 rows containing missing values

    ## Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
    ## Removed 36 rows containing missing values
    ## Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
    ## Removed 36 rows containing missing values
    ## Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
    ## Removed 36 rows containing missing values

    ## Warning: Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_point()`).
    ## Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 36 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-pvalue-pairplot-1.png)<!-- -->

### Pairplot of effect sizes across methods (GGpairs)

``` r
library(GGally)
se_mini |> 
  rowData() |> 
  as_tibble() |> 
  select(ql__logFC, lr_robust__logFC, voom__logFC, voom_weights__logFC, deseq2__log2FoldChange) |> 
  ggpairs(columns = 1:5) +
  my_theme +
  labs(title = "Pairplot of effect sizes across methods")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-effectsize-pairplot-1.png)<!-- -->

### Volcano Plots for Each Method

Visualising the significance and effect size of the differential
expression results as a volcano plots we appreciate that DESeq2 has much
lower p-values than other methods, for the same model.

``` r
# Create volcano plots
se_mini |>

    # Select the columns we want to plot
    rowData() |> 
    as_tibble(rownames = ".feature") |> 
    select(
            .feature,
      ql__logFC, ql__PValue,
      lr_robust__logFC, lr_robust__PValue,
      voom__logFC, voom__P.Value,
      voom_weights__logFC, voom_weights__P.Value,
      deseq2__log2FoldChange, deseq2__pvalue
    ) |>

    # Pivot longer to get a tidy data frame
    pivot_longer(
      - .feature,
      names_to = c("method", "stat"),
      values_to = "value", names_sep = "__"
    ) |>

    # Harmonize column names
    mutate(stat  = case_when(
        stat %in% c("logFC", "log2FoldChange") ~ "logFC",
        stat %in% c("PValue", "pvalue", "P.Value", "p.value") ~ "PValue"
    )) |>
  pivot_wider(names_from = "stat", values_from = "value") |>
  unnest( logFC, PValue) |> 

    # Plot
  ggplot(aes(x = logFC, y = PValue)) +
  geom_point(aes(color = PValue < 0.05, size = PValue < 0.05)) +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_size_manual(values = c("TRUE" = 0.5, "FALSE" = 0.1)) +
  facet_wrap(~method) +
  my_theme +
  labs(title = "Volcano Plots by Method")
```

    ## Warning: `unnest()` has a new interface. See `?unnest` for details.
    ## ℹ Try `df %>% unnest(c(logFC, PValue))`, with `mutate()` if needed.

    ## Warning: Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-volcano-plots-1-1.png)<!-- -->

Plotting independent y-axis scales for the p-values and effect sizes
allows us to compare the top genes across methods.

``` r
# Create volcano plots
se_mini |>

    # Select the columns we want to plot
    rowData() |> 
    as_tibble(rownames = ".feature") |> 
    select(
            .feature,
      ql__logFC, ql__PValue,
      lr_robust__logFC, lr_robust__PValue,
      voom__logFC, voom__P.Value,
      voom_weights__logFC, voom_weights__P.Value,
      deseq2__log2FoldChange, deseq2__pvalue
    ) |>

    # Pivot longer to get a tidy data frame
    pivot_longer(
      - .feature,
      names_to = c("method", "stat"),
      values_to = "value", names_sep = "__"
    ) |>

    # Harmonize column names
    mutate(stat  = case_when(
        stat %in% c("logFC", "log2FoldChange") ~ "logFC",
        stat %in% c("PValue", "pvalue", "P.Value", "p.value") ~ "PValue"
    )) |>
  pivot_wider(names_from = "stat", values_from = "value") |>
  unnest( logFC, PValue) |> 

    # Plot
  ggplot(aes(x = logFC, y = PValue)) +
  geom_point(aes(color = PValue < 0.05, size = PValue < 0.05)) +
  ggrepel::geom_text_repel(aes(label = .feature), size = 2, max.overlaps = 10) +
  scale_y_continuous(trans = tidybulk::log10_reverse_trans()) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_size_manual(values = c("TRUE" = 0.5, "FALSE" = 0.1)) +
  facet_wrap(~method, scales = "free_y") +
  my_theme +
  labs(title = "Volcano Plots by Method")
```

    ## Warning: `unnest()` has a new interface. See `?unnest` for details.
    ## ℹ Try `df %>% unnest(c(logFC, PValue))`, with `mutate()` if needed.

    ## Warning: Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 36 rows containing missing values or values outside the scale range
    ## (`geom_text_repel()`).

    ## Warning: ggrepel: 329 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 381 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 385 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 379 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 374 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/differential-expression-volcano-plots-2-1.png)<!-- -->

### Differential Expression with Contrasts

Contrast-based differential expression analysis using **edgeR**
([Robinson, McCarthy, and Smyth 2010](#ref-robinson2010edger); [Chen,
Lun, and Smyth 2016](#ref-chen2016edgeR)) is a standard statistical
approach for testing specific comparisons in complex designs.

``` r
# Using contrasts for more complex comparisons
se_mini |>
    test_differential_expression(
        ~ 0 + condition,                  
        .contrasts = c("conditionTRUE - conditionFALSE"),
        method = "edgeR_quasi_likelihood", 
        prefix = "contrasts__"
    ) |> 

    # Print the gene statistics
  pivot_transcript() |>
  select(contains("contrasts"))
```

    ## Warning: The `.contrasts` argument of `test_differential_abundance()` is deprecated as
    ## of tidybulk 1.7.4.
    ## ℹ The argument .contrasts is now deprecated please use contrasts (without the
    ##   dot).
    ## ℹ The deprecated feature was likely used in the tidybulk package.
    ##   Please report the issue at <https://github.com/stemangiola/tidybulk/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## tidybulk says: The design column names are "conditionFALSE, conditionTRUE"

    ## tidybulk says: to access the DE object do `metadata(.)$tidybulk$edgeR_quasi_likelihood_object`
    ## tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$edgeR_quasi_likelihood_fit`

    ## # A tibble: 394 × 5
    ##    contrasts__logFC___conditionT…¹ contrasts__logCPM___…² contrasts__F___condi…³
    ##                              <dbl>                  <dbl>                  <dbl>
    ##  1                           3.60                    9.31                  2.93 
    ##  2                           2.47                    5.94                 12.5  
    ##  3                           2.47                   12.9                   7.65 
    ##  4                           2.07                   10.7                   3.83 
    ##  5                           3.42                   12.7                   2.11 
    ##  6                           1.26                    5.65                  0.390
    ##  7                          -0.983                   8.98                  1.11 
    ##  8                          -8.43                   12.8                  29.0  
    ##  9                           2.21                   10.9                   1.97 
    ## 10                           0.784                   4.75                  0.390
    ## # ℹ 384 more rows
    ## # ℹ abbreviated names: ¹​contrasts__logFC___conditionTRUE...conditionFALSE,
    ## #   ²​contrasts__logCPM___conditionTRUE...conditionFALSE,
    ## #   ³​contrasts__F___conditionTRUE...conditionFALSE
    ## # ℹ 2 more variables: contrasts__PValue___conditionTRUE...conditionFALSE <dbl>,
    ## #   contrasts__FDR___conditionTRUE...conditionFALSE <dbl>

### Differential Expression with minimum fold change (TREAT method)

TREAT method ([McCarthy and Smyth 2009](#ref-mccarthy2009treat)) is used
for testing significance relative to a fold-change threshold.

``` r
# Using contrasts for more complex comparisons
se_mini |>
    test_differential_expression(
        ~ 0 + condition,                  
        .contrasts = c("conditionTRUE - conditionFALSE"),
        method = "edgeR_quasi_likelihood", 
        test_above_log2_fold_change = 2, 
        prefix = "treat__"
    ) |> 

    # Print the gene statistics
  pivot_transcript() |>
  select(contains("treat"))
```

    ## tidybulk says: The design column names are "conditionFALSE, conditionTRUE"

    ## tidybulk says: to access the DE object do `metadata(.)$tidybulk$edgeR_quasi_likelihood_object`
    ## tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$edgeR_quasi_likelihood_fit`

    ## # A tibble: 394 × 5
    ##    treat__logFC___conditionTRUE.…¹ treat__unshrunk.logF…² treat__logCPM___cond…³
    ##                              <dbl>                  <dbl>                  <dbl>
    ##  1                           3.60                   3.60                    9.31
    ##  2                           2.47                   2.49                    5.94
    ##  3                           2.47                   2.47                   12.9 
    ##  4                           2.07                   2.07                   10.7 
    ##  5                           3.42                   3.42                   12.7 
    ##  6                           1.26                   1.27                    5.65
    ##  7                          -0.983                 -0.983                   8.98
    ##  8                          -8.43                  -8.43                   12.8 
    ##  9                           2.21                   2.21                   10.9 
    ## 10                           0.784                  0.790                   4.75
    ## # ℹ 384 more rows
    ## # ℹ abbreviated names: ¹​treat__logFC___conditionTRUE...conditionFALSE,
    ## #   ²​treat__unshrunk.logFC___conditionTRUE...conditionFALSE,
    ## #   ³​treat__logCPM___conditionTRUE...conditionFALSE
    ## # ℹ 2 more variables: treat__PValue___conditionTRUE...conditionFALSE <dbl>,
    ## #   treat__FDR___conditionTRUE...conditionFALSE <dbl>

### Mixed Models for Complex Designs

glmmSeq ([Ma et al. 2020](#ref-ma2020glmmseq)) is used for generalized
linear mixed models for RNA-seq data.

``` r
# Using glmmSeq for mixed models
se_mini = se_mini |>
  keep_abundant(formula_design = ~ condition) |>
  test_differential_expression(
    ~ condition + (1|Cell.type), 
    method = "glmmseq_lme4", 
    prefix = "glmmseq__"
  ) 
```

    ## 
    ## n = 5 samples, 4 individuals

    ## Time difference of 1.962443 mins

    ## Errors in 1 gene(s): PTGER2

    ## tidybulk says: to access the DE object do `attr(..., "internals")$glmmseq_lme4_object`
    ## tidybulk says: to access the raw results (fitted GLM) do `attr(..., "internals")$glmmseq_lme4_fit`

``` r
  se_mini |>
  pivot_transcript() 
```

    ## # A tibble: 391 × 56
    ##    .feature entrez gene_name .abundant ql__logFC ql__logCPM  ql__F ql__PValue
    ##    <chr>    <chr>  <chr>     <lgl>         <dbl>      <dbl>  <dbl>      <dbl>
    ##  1 ABCB4    5244   ABCB4     TRUE          3.60        9.31  2.93     0.139  
    ##  2 ABCB9    23457  ABCB9     TRUE          2.47        5.94 12.5      0.0112 
    ##  3 ACAP1    9744   ACAP1     TRUE          2.47       12.9   7.65     0.0341 
    ##  4 ACP5     54     ACP5      TRUE          2.07       10.7   3.83     0.1000 
    ##  5 ADAM28   10863  ADAM28    TRUE          3.42       12.7   2.11     0.199  
    ##  6 ADAMDEC1 27299  ADAMDEC1  TRUE          1.26        5.65  0.390    0.554  
    ##  7 ADRB2    154    ADRB2     TRUE         -0.983       8.98  1.11     0.333  
    ##  8 AIF1     199    AIF1      TRUE         -8.43       12.8  29.0      0.00168
    ##  9 AIM2     9447   AIM2      TRUE          2.21       10.9   1.97     0.212  
    ## 10 ALOX15   246    ALOX15    TRUE          0.784       4.75  0.390    0.553  
    ## # ℹ 381 more rows
    ## # ℹ 48 more variables: ql__FDR <dbl>, lr_robust__logFC <dbl>,
    ## #   lr_robust__logCPM <dbl>, lr_robust__LR <dbl>, lr_robust__PValue <dbl>,
    ## #   lr_robust__FDR <dbl>, deseq2__baseMean <dbl>, deseq2__log2FoldChange <dbl>,
    ## #   deseq2__lfcSE <dbl>, deseq2__stat <dbl>, deseq2__pvalue <dbl>,
    ## #   deseq2__padj <dbl>, voom__logFC <dbl>, voom__AveExpr <dbl>, voom__t <dbl>,
    ## #   voom__P.Value <dbl>, voom__adj.P.Val <dbl>, voom__B <dbl>, …

### Gene Description

With tidybulk, retrieving gene descriptions is straightforward, making
it easy to enhance the interpretability of your differential expression
results.

``` r
# Add gene descriptions using the original SummarizedExperiment
se_mini |> 

    describe_transcript() |>

    # Filter top significant genes
    filter(ql__FDR < 0.05) |>

    # Print the gene statistics
    pivot_transcript() |> 
    dplyr::select(.feature, description, contains("ql")) |> 
    head()
```

    ## 

    ## 

    ## # A tibble: 6 × 7
    ##   .feature description             ql__logFC ql__logCPM ql__F ql__PValue ql__FDR
    ##   <chr>    <chr>                       <dbl>      <dbl> <dbl>      <dbl>   <dbl>
    ## 1 ABCB9    ATP binding cassette s…      2.47       5.94  12.5 0.0112     0.0332 
    ## 2 AIF1     allograft inflammatory…     -8.43      12.8   29.0 0.00168    0.0107 
    ## 3 ANGPT4   angiopoietin 4              -2.96       5.62  14.4 0.00764    0.0257 
    ## 4 APOBEC3A apolipoprotein B mRNA …     -8.92      11.1   80.0 0.0000692  0.00317
    ## 5 AQP9     aquaporin 9                -11.1       10.3   93.1 0.00000500 0.00197
    ## 6 ASGR1    asialoglycoprotein rec…     -9.29       9.07  24.1 0.00105    0.00787

## Step 4: Batch Effect Correction

ComBat-seq ([Zhang, Parmigiani, and Johnson
2020](#ref-zhang2017combatseq)) is used for batch effect correction in
RNA-seq data.

``` r
# Adjust for batch effects
se_mini = se_mini |>
  adjust_abundance(
      .factor_unwanted = time, 
      .factor_of_interest = condition, 
    method = "combat_seq", 
      abundance = "counts_tmm"
  )
```

    ## Found 2 batches
    ## Using null model in ComBat-seq.
    ## Adjusting for 1 covariate(s) or covariate level(s)
    ## Estimating dispersions
    ## Fitting the GLM model
    ## Shrinkage off - using GLM estimates for parameters
    ## Adjusting the data

``` r
# Scatter plot of adjusted vs unadjusted
left_join(
    se_mini |> assay("counts_tmm") |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts_tmm") ,
    se_mini |> assay("counts_tmm_adjusted") |> as_tibble(rownames = ".feature") |> pivot_longer(cols = -.feature, names_to = ".sample", values_to = "counts_tmm_adjusted") ,
    by = c(".feature", ".sample")
  ) |>
  ggplot(aes(x = counts_tmm + 1, y = counts_tmm_adjusted + 1)) +
  geom_point(aes(color = .sample), size = 0.1) +
  ggrepel::geom_text_repel(aes(label = .feature), size = 2, max.overlaps = 10) +
  scale_x_log10() +
  scale_y_log10() +
  my_theme +
  labs(title = "Scatter plot of adjusted vs unadjusted")
```

    ## Warning: ggrepel: 1927 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/batch-correction-adjust-abundance-1.png)<!-- -->

## Step 5: Cellularity Analysis

CIBERSORT ([Newman et al. 2015](#ref-newman2015cibersort)) is used for
cell type deconvolution.

Cellularity deconvolution is a standard approach for estimating the
cellular composition of a sample.

### Available Deconvolution Methods

The `tidybulk` package provides several methods for deconvolution:

- **CIBERSORT** ([Newman et al. 2015](#ref-newman2015cibersort)): Uses
  support vector regression to deconvolve cell type proportions.
  Requires the `class`, `e1071`, and `preprocessCore` packages.
- **LLSR** ([Newman et al. 2015](#ref-newman2015cibersort)): Linear
  Least Squares Regression for deconvolution.
- **EPIC** ([**racle2017epic?**](#ref-racle2017epic)): Uses a
  reference-based approach to estimate cell fractions.
- **MCP-counter** ([**becht2016mcp?**](#ref-becht2016mcp)): Quantifies
  the abundance of immune and stromal cell populations.
- **quanTIseq**
  ([**finotello2019quantiseq?**](#ref-finotello2019quantiseq)): A
  computational framework for inferring the immune contexture of tumors.
- **xCell** ([**aran2017xcell?**](#ref-aran2017xcell)): Performs cell
  type enrichment analysis.

### Example Usage

``` r
se_mini = 

se_mini |> 
deconvolve_cellularity(method = "cibersort", cores = 1, prefix = "cibersort__") 
```

For the rest of the methods, you need to install the `immunedeconv`
package.

``` r
if (!requireNamespace("immunedeconv")) BiocManager::install("immunedeconv")
```

``` r
se_mini = 

se_mini |> 
# Example using LLSR
deconvolve_cellularity(method = "llsr", prefix = "llsr__") |> 

# Example using EPIC
deconvolve_cellularity(method = "epic", prefix = "epic__") |> 

# Example using MCP-counter
deconvolve_cellularity(method = "mcp_counter", prefix = "mcp__") |> 

# Example using quanTIseq
deconvolve_cellularity(method = "quantiseq", prefix = "quantiseq__") |> 

# Example using xCell
deconvolve_cellularity(method = "xcell", prefix = "xcell__")
```

### Plotting Results

Visualize the cell type proportions as a stacked barplot for each
method:

``` r
# Visualize CIBERSORT results
se_mini  |>
  pivot_sample() |>
  select(.sample, contains("cibersort__")) |>
  pivot_longer(cols = -1, names_to = "Cell_type_inferred", values_to = "proportion") |>
  ggplot(aes(x = .sample, y = proportion, fill = Cell_type_inferred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "CIBERSORT Cell Type Proportions")
```

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/deconvolution-plotting-1.png)<!-- -->

``` r
 # Repeat similar plotting for LLSR, EPIC, MCP-counter, quanTIseq, and xCell
se_mini  |>
  pivot_sample() |>
  select(.sample, contains("llsr__")) |>
  pivot_longer(cols = -1, names_to = "Cell_type_inferred", values_to = "proportion") |>
  ggplot(aes(x = .sample, y = proportion, fill = Cell_type_inferred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "LLSR Cell Type Proportions")

  se_mini    |>
  pivot_sample() |>
  select(.sample, contains("epic__")) |>
  pivot_longer(cols = -1, names_to = "Cell_type_inferred", values_to = "proportion") |>
  ggplot(aes(x = .sample, y = proportion, fill = Cell_type_inferred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "EPIC Cell Type Proportions")

  se_mini    |>
  pivot_sample() |>
  select(.sample, contains("mcp__")) |>
  pivot_longer(cols = -1, names_to = "Cell_type_inferred", values_to = "proportion") |>
  ggplot(aes(x = .sample, y = proportion, fill = Cell_type_inferred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "MCP-counter Cell Type Proportions")

  se_mini    |>
  pivot_sample() |>
  select(.sample, contains("quantiseq__")) |>
  pivot_longer(cols = -1, names_to = "Cell_type_inferred", values_to = "proportion") |>
  ggplot(aes(x = .sample, y = proportion, fill = Cell_type_inferred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "quanTIseq Cell Type Proportions")

  se_mini    |>
  pivot_sample() |>
  select(.sample, contains("xcell__")) |>
  pivot_longer(cols = -1, names_to = "Cell_type_inferred", values_to = "proportion") |>
  ggplot(aes(x = .sample, y = proportion, fill = Cell_type_inferred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "xCell Cell Type Proportions")
```

## Step 6: Gene Enrichment Analysis

Gene Set Enrichment Analysis (GSEA) ([Subramanian et al.
2005](#ref-subramanian2005gsea)) is used for gene set enrichment.

``` r
# Run gene rank enrichment (GSEA style)
gene_rank_res =
  se_mini |>

    # Filter for genes with entrez IDs
  filter(!entrez |> is.na()) |>

  # Test gene rank
  test_gene_rank(
    .entrez = entrez,
    .arrange_desc = lr_robust__logFC,
    species = "Homo sapiens",
    gene_sets = c("H", "C2", "C5")
  )
```

    ## 

    ## using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `fit = map(...)`.
    ## Caused by warning in `fgseaMultilevel()`:
    ## ! For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.

``` r
# Inspect significant gene sets (example for C2 collection)
gene_rank_res |>
  filter(gs_collection == "C2") |>
  dplyr::select(-fit) |>
  unnest(test) |>
  filter(p.adjust < 0.05)
```

    ## # A tibble: 65 × 13
    ##    gs_collection idx_for_plotting ID         Description setSize enrichmentScore
    ##    <chr>                    <int> <chr>      <chr>         <int>           <dbl>
    ##  1 C2                           1 REACTOME_… REACTOME_N…      32          -0.770
    ##  2 C2                           2 SMID_BREA… SMID_BREAS…     105           0.539
    ##  3 C2                           3 FULCHER_I… FULCHER_IN…      42          -0.684
    ##  4 C2                           4 REACTOME_… REACTOME_I…      73          -0.586
    ##  5 C2                           5 RUTELLA_R… RUTELLA_RE…      36          -0.679
    ##  6 C2                           6 CHEN_META… CHEN_METAB…      67          -0.570
    ##  7 C2                           7 SMID_BREA… SMID_BREAS…      82           0.504
    ##  8 C2                           8 LIU_OVARI… LIU_OVARIA…     130          -0.475
    ##  9 C2                           9 RUTELLA_R… RUTELLA_RE…      39          -0.644
    ## 10 C2                          10 SMIRNOV_C… SMIRNOV_CI…      23          -0.737
    ## # ℹ 55 more rows
    ## # ℹ 7 more variables: NES <dbl>, pvalue <dbl>, p.adjust <dbl>, qvalue <dbl>,
    ## #   rank <dbl>, leading_edge <chr>, core_enrichment <chr>

# Visualize enrichment

``` r
  library(enrichplot)
```

    ## Warning: package 'enrichplot' was built under R version 4.5.1

    ## enrichplot v1.28.4 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
    ## clusterProfiler: an R package for comparing biological themes among
    ## gene clusters. OMICS: A Journal of Integrative Biology. 2012,
    ## 16(5):284-287

    ## 
    ## Attaching package: 'enrichplot'

    ## The following object is masked from 'package:GGally':
    ## 
    ##     ggtable

``` r
  library(patchwork)
  gene_rank_res |>
    unnest(test) |>
    head() |>
    mutate(plot = pmap(
      list(fit, ID, idx_for_plotting, p.adjust),
      ~ enrichplot::gseaplot2(
        ..1, geneSetID = ..3,
        title = sprintf("%s \nadj pvalue %s", ..2, round(..4, 2)),
        base_size = 6, rel_heights = c(1.5, 0.5), subplots = c(1, 2)
      )
    )) |>
    pull(plot) 
```

    ## Warning: There were 2 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `plot = pmap(...)`.
    ## Caused by warning:
    ## ! `aes_()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`
    ## ℹ The deprecated feature was likely used in the enrichplot package.
    ##   Please report the issue at
    ##   <https://github.com/GuangchuangYu/enrichplot/issues>.
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

    ## [[1]]

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/enrichment-visualize-gsea-plots-1.png)<!-- -->

    ## 
    ## [[2]]

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/enrichment-visualize-gsea-plots-2.png)<!-- -->

    ## 
    ## [[3]]

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/enrichment-visualize-gsea-plots-3.png)<!-- -->

    ## 
    ## [[4]]

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/enrichment-visualize-gsea-plots-4.png)<!-- -->

    ## 
    ## [[5]]

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/enrichment-visualize-gsea-plots-5.png)<!-- -->

    ## 
    ## [[6]]

![](/Users/a1234450/Documents/GitHub/tidybulk/README_files/figure-gfm/enrichment-visualize-gsea-plots-6.png)<!-- -->

Gene Ontology overrepresentation analysis ([Ashburner et al.
2000](#ref-ashburner2000go)) is used for functional enrichment.

``` r
# Test gene overrepresentation
se_mini_overrep = 
  se_mini |>
  
  # Label genes to test overrepresentation of
  mutate(genes_to_test = ql__FDR < 0.05) |>
  
    # Filter for genes with entrez IDs
  filter(!entrez |> is.na()) |>
  
  test_gene_overrepresentation(
    .entrez = entrez,
    species = "Homo sapiens",
    .do_test = genes_to_test,
    gene_sets = c("H", "C2", "C5")
  )

  se_mini_overrep
```

    ## # A tibble: 1,450 × 13
    ##    gs_collection ID      Description GeneRatio BgRatio RichFactor FoldEnrichment
    ##    <chr>         <chr>   <chr>       <chr>     <chr>        <dbl>          <dbl>
    ##  1 C2            SMID_B… SMID_BREAS… 36/152    484/22…     0.0744          11.0 
    ##  2 C5            GOBP_I… GOBP_IMMUN… 31/148    380/19…     0.0816          10.8 
    ##  3 C2            JAATIN… JAATINEN_H… 24/152    235/22…     0.102           15.0 
    ##  4 C2            MCLACH… MCLACHLAN_… 24/152    264/22…     0.0909          13.4 
    ##  5 C5            GOCC_E… GOCC_EXTER… 29/148    403/19…     0.0720           9.51
    ##  6 C2            POOLA_… POOLA_INVA… 24/152    292/22…     0.0822          12.1 
    ##  7 C5            GOMF_I… GOMF_IMMUN… 20/148    157/19…     0.127           16.8 
    ##  8 C2            BLANCO… BLANCO_MEL… 19/152    174/22…     0.109           16.1 
    ##  9 C5            GOBP_L… GOBP_LEUKO… 22/148    236/19…     0.0932          12.3 
    ## 10 C5            GOBP_L… GOBP_LEUKO… 26/148    396/19…     0.0657           8.67
    ## # ℹ 1,440 more rows
    ## # ℹ 6 more variables: zScore <dbl>, pvalue <dbl>, p.adjust <dbl>, qvalue <dbl>,
    ## #   Count <int>, entrez <list>

EGSEA ([Alhamdoosh et al. 2017](#ref-alhamdoosh2017egsea)) is used for
ensemble gene set enrichment analysis. EGSEA is a method that combines
multiple gene set enrichment analysis methods to provide a more robust
and comprehensive analysis of gene set enrichment. It creates a
web-based interactive tool that allows you to explore the results of the
gene set enrichment analysis.

``` r
library(EGSEA)
```

    ## Loading required package: gage

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: topGO

    ## Warning: package 'topGO' was built under R version 4.5.1

    ## Loading required package: graph

    ## 
    ## Attaching package: 'graph'

    ## The following object is masked from 'package:stringr':
    ## 
    ##     boundary

    ## Loading required package: GO.db

    ## Loading required package: SparseM

    ## 
    ## groupGOTerms:    GOBPTerm, GOMFTerm, GOCCTerm environments built.

    ## 
    ## Attaching package: 'topGO'

    ## The following object is masked from 'package:gage':
    ## 
    ##     geneData

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     members

    ## Loading required package: pathview

    ## ##############################################################################
    ## Pathview is an open source software package distributed under GNU General
    ## Public License version 3 (GPLv3). Details of GPLv3 is available at
    ## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    ## formally cite the original Pathview paper (not just mention it) in publications
    ## or products. For details, do citation("pathview") within R.
    ## 
    ## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    ## license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ## ##############################################################################

    ## 

    ## 

    ## 

``` r
# Test gene enrichment
  se_mini |> 

  # Filter for genes with entrez IDs
  filter(!entrez |> is.na()) |>

  # Test gene enrichment
  test_gene_enrichment(
    .formula = ~condition,
    .entrez = entrez,
    species = "human", 
    gene_sets = "h"
  )
```

    ## tidybulk says: The design column names are "(Intercept), conditionTRUE"

    ## [1] "Loading MSigDB Gene Sets ... "
    ## [1] "Loaded gene sets for the collection h ..."
    ## [1] "Indexed the collection h ..."
    ## [1] "Created annotation for the collection h ..."

    ## EGSEA analysis has started

    ## ##------ Mon Jul 28 13:59:55 2025 ------##

    ## The argument 'contrast' is recommended to be a matrix object.
    ## See Vignette or Help.

    ## Log fold changes are estimated using limma package ...

    ## limma DE analysis is carried out ...

    ## EGSEA is running on the provided data and h collection

    ## 

    ## ##------ Mon Jul 28 13:59:56 2025 ------##

    ## EGSEA analysis took 1.55 seconds.

    ## EGSEA analysis has completed

    ## EGSEA HTML report is being generated ...

    ## ##------ Mon Jul 28 13:59:56 2025 ------##

    ## Report pages and figures are being generated for the h collection ...

    ##    Heat maps are being generated for top-ranked gene sets 
    ## based on logFC ...

    ##    Summary plots are being generated ...

    ##    Comparison summary plots are being generated  ...

    ## ##------ Mon Jul 28 14:01:28 2025 ------##

    ## EGSEA report generation took 91.398 seconds.

    ## EGSEA report has been generated.

    ## # A tibble: 0 × 0

## Bibliography

`tidybulk` allows you to get the bibliography of all methods used in our
workflow.

``` r
# Get bibliography of all methods used in our workflow
se_mini |> get_bibliography()
```

    ##  @Article{tidybulk,
    ##   title = {tidybulk: an R tidy framework for modular transcriptomic data analysis},
    ##   author = {Stefano Mangiola and Ramyar Molania and Ruining Dong and Maria A. Doyle & Anthony T. Papenfuss},
    ##   journal = {Genome Biology},
    ##   year = {2021},
    ##   volume = {22},
    ##   number = {42},
    ##   url = {https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02233-7},
    ##   }
    ## @article{wickham2019welcome,
    ##   title={Welcome to the Tidyverse},
    ##   author={Wickham, Hadley and Averick, Mara and Bryan, Jennifer and Chang, Winston and McGowan, Lucy D'Agostino and Francois, Romain and Grolemund, Garrett and Hayes, Alex and Henry, Lionel and Hester, Jim and others},
    ##   journal={Journal of Open Source Software},
    ##   volume={4},
    ##   number={43},
    ##   pages={1686},
    ##   year={2019}
    ##  }
    ## @article{leek2012sva,
    ##   title={The sva package for removing batch effects and other unwanted variation in high-throughput experiments},
    ##   author={Leek, Jeffrey T and Johnson, W Evan and Parker, Hilary S and Jaffe, Andrew E and Storey, John D},
    ##   journal={Bioinformatics},
    ##   volume={28},
    ##   number={6},
    ##   pages={882--883},
    ##   year={2012},
    ##   publisher={Oxford University Press}
    ##  }
    ## @article{newman2015robust,
    ##   title={Robust enumeration of cell subsets from tissue expression profiles},
    ##   author={Newman, Aaron M and Liu, Chih Long and Green, Michael R and Gentles, Andrew J and Feng, Weiguo and Xu, Yue and Hoang, Chuong D and Diehn, Maximilian and Alizadeh, Ash A},
    ##   journal={Nature methods},
    ##   volume={12},
    ##   number={5},
    ##   pages={453--457},
    ##   year={2015},
    ##   publisher={Nature Publishing Group}
    ##  }

``` r
sessionInfo()
```

    ## R version 4.5.0 (2025-04-11)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Sonoma 14.6.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Australia/Adelaide
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] EGSEA_1.36.0                    pathview_1.48.0                
    ##  [3] topGO_2.60.1                    SparseM_1.84-2                 
    ##  [5] GO.db_3.21.0                    graph_1.86.0                   
    ##  [7] AnnotationDbi_1.70.0            gage_2.58.0                    
    ##  [9] patchwork_1.3.1                 enrichplot_1.28.4              
    ## [11] GGally_2.3.0                    DESeq2_1.48.1                  
    ## [13] edgeR_4.6.3                     limma_3.64.1                   
    ## [15] tidySummarizedExperiment_1.19.4 SummarizedExperiment_1.38.1    
    ## [17] Biobase_2.68.0                  GenomicRanges_1.60.0           
    ## [19] GenomeInfoDb_1.44.1             IRanges_2.42.0                 
    ## [21] S4Vectors_0.46.0                BiocGenerics_0.54.0            
    ## [23] generics_0.1.4                  MatrixGenerics_1.20.0          
    ## [25] matrixStats_1.5.0               tidybulk_1.99.1                
    ## [27] ttservice_0.5.3                 ggrepel_0.9.6                  
    ## [29] magrittr_2.0.3                  lubridate_1.9.4                
    ## [31] forcats_1.0.0                   stringr_1.5.1                  
    ## [33] dplyr_1.1.4                     purrr_1.1.0                    
    ## [35] readr_2.1.5                     tidyr_1.3.1                    
    ## [37] tibble_3.3.0                    ggplot2_3.5.2.9002             
    ## [39] tidyverse_2.0.0                 knitr_1.50                     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fs_1.6.6                    GSVA_2.2.0                 
    ##   [3] bitops_1.0-9                R2HTML_2.3.4               
    ##   [5] httr_1.4.7                  RColorBrewer_1.1-3         
    ##   [7] numDeriv_2016.8-1.1         Rgraphviz_2.52.0           
    ##   [9] doRNG_1.8.6.2               tools_4.5.0                
    ##  [11] backports_1.5.0             utf8_1.2.6                 
    ##  [13] R6_2.6.1                    DT_0.33                    
    ##  [15] HDF5Array_1.36.0            sn_2.1.1                   
    ##  [17] lazyeval_0.2.2              mgcv_1.9-3                 
    ##  [19] rhdf5filters_1.20.0         withr_3.0.2                
    ##  [21] preprocessCore_1.70.0       cli_3.6.5                  
    ##  [23] sandwich_3.1-1              labeling_0.4.3             
    ##  [25] KEGGgraph_1.68.0            mvtnorm_1.3-3              
    ##  [27] S7_0.2.0                    genefilter_1.90.0          
    ##  [29] PADOG_1.50.0                proxy_0.4-27               
    ##  [31] yulab.utils_0.2.0           gson_0.1.0                 
    ##  [33] DOSE_4.2.0                  R.utils_2.13.0             
    ##  [35] HTMLUtils_0.1.9             plotrix_3.8-4              
    ##  [37] rstudioapi_0.17.1           RSQLite_2.4.2              
    ##  [39] gridGraphics_0.5-1          hwriter_1.3.2.1            
    ##  [41] gtools_3.9.5                Matrix_1.7-3               
    ##  [43] abind_1.4-8                 R.methodsS3_1.8.2          
    ##  [45] lifecycle_1.0.4             multcomp_1.4-28            
    ##  [47] yaml_2.3.10                 mathjaxr_1.8-0             
    ##  [49] KEGGdzPathwaysGEO_1.46.0    gplots_3.2.0               
    ##  [51] rhdf5_2.52.1                qvalue_2.40.0              
    ##  [53] SparseArray_1.8.1           grid_4.5.0                 
    ##  [55] blob_1.2.4                  crayon_1.5.3               
    ##  [57] ggtangle_0.0.7              lattice_0.22-7             
    ##  [59] beachmat_2.24.0             msigdbr_25.1.1             
    ##  [61] cowplot_1.2.0               annotate_1.86.1            
    ##  [63] KEGGREST_1.48.1             magick_2.8.7               
    ##  [65] pillar_1.11.0               fgsea_1.34.2               
    ##  [67] hgu133plus2.db_3.13.0       hgu133a.db_3.13.0          
    ##  [69] rjson_0.2.23                widyr_0.1.5                
    ##  [71] codetools_0.2-20            fastmatch_1.1-6            
    ##  [73] mutoss_0.1-13               glue_1.8.0                 
    ##  [75] ggfun_0.2.0                 data.table_1.17.8          
    ##  [77] Rdpack_2.6.4                vctrs_0.6.5                
    ##  [79] png_0.1-8                   treeio_1.32.0              
    ##  [81] org.Mm.eg.db_3.21.0         gtable_0.3.6               
    ##  [83] org.Rn.eg.db_3.21.0         assertthat_0.2.1           
    ##  [85] cachem_1.1.0                xfun_0.52                  
    ##  [87] rbibutils_2.3               S4Arrays_1.8.1             
    ##  [89] survival_3.8-3              SingleCellExperiment_1.30.1
    ##  [91] iterators_1.0.14            statmod_1.5.0              
    ##  [93] TH.data_1.1-3               ellipsis_0.3.2             
    ##  [95] nlme_3.1-168                ggtree_3.16.3              
    ##  [97] bit64_4.6.0-1               rprojroot_2.1.0            
    ##  [99] SnowballC_0.7.1             irlba_2.3.5.1              
    ## [101] KernSmooth_2.23-26          colorspace_2.1-1           
    ## [103] DBI_1.2.3                   mnormt_2.1.1               
    ## [105] tidyselect_1.2.1            bit_4.6.0                  
    ## [107] compiler_4.5.0              curl_6.4.0                 
    ## [109] h5mread_1.0.1               TFisher_0.2.0              
    ## [111] DelayedArray_0.34.1         plotly_4.11.0              
    ## [113] scales_1.4.0                caTools_1.18.3             
    ## [115] SpatialExperiment_1.18.1    digest_0.6.37              
    ## [117] rmarkdown_2.29              XVector_0.48.0             
    ## [119] htmltools_0.5.8.1           pkgconfig_2.0.3            
    ## [121] sparseMatrixStats_1.20.0    fastmap_1.2.0              
    ## [123] rlang_1.1.6                 htmlwidgets_1.6.4          
    ## [125] UCSC.utils_1.4.0            farver_2.1.2               
    ## [127] zoo_1.8-14                  jsonlite_2.0.0             
    ## [129] BiocParallel_1.42.1         GOSemSim_2.34.0            
    ## [131] tokenizers_0.3.0            R.oo_1.27.1                
    ## [133] BiocSingular_1.24.0         RCurl_1.98-1.17            
    ## [135] GenomeInfoDbData_1.2.14     ggplotify_0.1.2            
    ## [137] Rhdf5lib_1.30.0             Rcpp_1.1.0                 
    ## [139] ape_5.8-1                   babelgene_22.9             
    ## [141] stringi_1.8.7               MASS_7.3-65                
    ## [143] globaltest_5.62.0           plyr_1.8.9                 
    ## [145] org.Hs.eg.db_3.21.0         ggstats_0.10.0             
    ## [147] parallel_4.5.0              Biostrings_2.76.0          
    ## [149] splines_4.5.0               multtest_2.64.0            
    ## [151] hms_1.1.3                   qqconf_1.3.2               
    ## [153] locfit_1.5-9.12             igraph_2.1.4               
    ## [155] rngtools_1.5.2              EGSEAdata_1.36.0           
    ## [157] reshape2_1.4.4              ScaledMatrix_1.16.0        
    ## [159] XML_3.99-0.18               GSA_1.03.3                 
    ## [161] evaluate_1.0.4              metap_1.12                 
    ## [163] tidytext_0.4.2              foreach_1.5.2              
    ## [165] tzdb_0.5.0                  rsvd_1.0.5                 
    ## [167] broom_1.0.8                 xtable_1.8-4               
    ## [169] e1071_1.7-16                tidytree_0.4.6             
    ## [171] janeaustenr_1.0.0           viridisLite_0.4.2          
    ## [173] class_7.3-23                clusterProfiler_4.16.0     
    ## [175] aplot_0.2.8                 safe_3.48.0                
    ## [177] memoise_2.0.1               timechange_0.3.0           
    ## [179] sva_3.56.0                  GSEABase_1.70.0

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-alhamdoosh2017egsea" class="csl-entry">

Alhamdoosh, M, M Ng, NJ Wilson, JM Sheridan, H Huynh, MJ Wilson, and ME
Ritchie. 2017. “Combining Multiple Tools Outperforms Individual Methods
in Gene Set Enrichment Analyses.” *Bioinformatics* 33 (3): 414–24.

</div>

<div id="ref-anders2010rle" class="csl-entry">

Anders, Simon, and Wolfgang Huber. 2010. “Differential Expression
Analysis for Sequence Count Data.” *Genome Biology* 11 (10): R106.

</div>

<div id="ref-ashburner2000go" class="csl-entry">

Ashburner, Michael, Catherine A Ball, Judith A Blake, David Botstein,
Heather Butler, J Michael Cherry, J Michael Davis, et al. 2000. “Gene
Ontology: Tool for the Unification of Biology.” *Nature Genetics* 25
(1): 25–29.

</div>

<div id="ref-bullard2010uq" class="csl-entry">

Bullard, James H, Elizabeth Purdom, Kasper D Hansen, and Sandrine
Dudoit. 2010. “Evaluation of Statistical Methods for Normalization and
Differential Expression in mRNA-Seq Experiments.” *BMC Bioinformatics*
11 (1): 94.

</div>

<div id="ref-chen2016edgeR" class="csl-entry">

Chen, Yunshun, Aaron TL Lun, and Gordon K Smyth. 2016. “From Reads to
Genes to Pathways: Differential Expression Analysis of RNA-Seq
Experiments Using Rsubread and the edgeR Quasi-Likelihood Pipeline.”
*F1000Research* 5.

</div>

<div id="ref-hotelling1933pca" class="csl-entry">

Hotelling, Harold. 1933. “Analysis of a Complex of Statistical Variables
into Principal Components.” *Journal of Educational Psychology* 24 (6):
417.

</div>

<div id="ref-kruskal1964mds" class="csl-entry">

Kruskal, Joseph B. 1964. “Multidimensional Scaling by Optimizing
Goodness of Fit to a Nonmetric Hypothesis.” *Psychometrika* 29 (1):
1–27.

</div>

<div id="ref-law2014voom" class="csl-entry">

Law, Charity W, Yunshun Chen, Wei Shi, and Gordon K Smyth. 2014. “Voom:
Precision Weights Unlock Linear Model Analysis Tools for RNA-Seq Read
Counts.” *Genome Biology* 15 (2): R29.

</div>

<div id="ref-liu2015voomweights" class="csl-entry">

Liu, Rui, Anna Z Holik, Shian Su, Nikolas Jansz, Yunshun Chen, Hsiu S
Leong, Marnie E Blewitt, Marie-Liesse Asselin-Labat, Gordon K Smyth, and
Matthew E Ritchie. 2015. “Why Weight? Modelling Sample and Observational
Level Variability Improves Power in RNA-Seq Analyses.” *Nucleic Acids
Research* 43 (15): e97–97.

</div>

<div id="ref-love2014deseq2" class="csl-entry">

Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12): 550.

</div>

<div id="ref-ma2020glmmseq" class="csl-entry">

Ma, Shuang, Bing Zhang, Laura M LaFave, Amanda S Earl, Zhicheng Chiang,
Yuling Hu, Jun Ding, et al. 2020. “A General Linear Mixed Model for
Differential Expression Analysis of Single-Cell RNA-Seq Data.” *Nature
Communications* 11 (1): 1–11.

</div>

<div id="ref-macqueen1967kmeans" class="csl-entry">

MacQueen, J. 1967. “Some Methods for Classification and Analysis of
Multivariate Observations.” In *Proceedings of the Fifth Berkeley
Symposium on Mathematical Statistics and Probability*, 1:281–97. 14.
Oakland, CA, USA.

</div>

<div id="ref-mccarthy2009treat" class="csl-entry">

McCarthy, Davis J, and Gordon K Smyth. 2009. “Testing Significance
Relative to a Fold-Change Threshold Is a TREAT.” *Bioinformatics* 25
(6): 765–71.

</div>

<div id="ref-newman2015cibersort" class="csl-entry">

Newman, Aaron M, Chih Long Liu, Michael R Green, Andrew J Gentles,
Weiguo Feng, Yue Xu, Chuong D Hoang, Maximilian Diehn, and Ash A
Alizadeh. 2015. “Robust Enumeration of Cell Subsets from Tissue
Expression Profiles.” *Nature Methods* 12 (5): 453–57.

</div>

<div id="ref-robinson2010edger" class="csl-entry">

Robinson, Mark D, Davis J McCarthy, and Gordon K Smyth. 2010. “edgeR: A
Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data.” *Bioinformatics* 26 (1): 139–40.

</div>

<div id="ref-subramanian2005gsea" class="csl-entry">

Subramanian, Aravind, Pablo Tamayo, Vamsi K Mootha, Sayan Mukherjee,
Benjamin L Ebert, Michael A Gillette, Amanda Paulovich, et al. 2005.
“Gene Set Enrichment Analysis: A Knowledge-Based Approach for
Interpreting Genome-Wide Expression Profiles.” *Proceedings of the
National Academy of Sciences* 102 (43): 15545–50.

</div>

<div id="ref-zhang2017combatseq" class="csl-entry">

Zhang, Yu, Giovanni Parmigiani, and W Evan Johnson. 2020. “ComBat-Seq:
Batch Effect Adjustment for RNA-Seq Count Data.” *NAR Genomics and
Bioinformatics* 2 (3): lqaa078.

</div>

</div>
