ttBulk - part of tidyTranscriptomics
================

<!-- badges: start --> [![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

<!---

[![Build Status](https://travis-ci.org/stemangiola/ttBulk.svg?branch=master)](https://travis-ci.org/stemangiola/ttBulk) [![Coverage Status](https://coveralls.io/repos/github/stemangiola/ttBulk/badge.svg?branch=master)](https://coveralls.io/github/stemangiola/ttBulk?branch=master)

-->

A modular framework and a tidy data structure for bulk transcriptional
analyses. This package is part of a tidy transcriptomics suite for tidy
and user-friendly grammar of RNA sequencing data exploration and
analysis. A wide range of analyses are included in convenient wrappers,
including: resolving duplicate gene symbol (e.g., isoforms), scaling
(normalisation), identifying variable transcripts, removal of know
unwanted variation, differential transcript abundance analyses
(differential expression), reducing data dimensionality (MDS, PCA and
tSNE), clustering (SNN and kmeans), gene enrichment analyses, remove
redundancy (either samples or transcripts).

# <img src="inst/logo.png" height="139px" width="120px" />

# Installation

Due to Bioconductor submission requiring R \>= 4.0 you need to install
the **dev** branch

``` r
devtools::install_github("stemangiola/ttBulk@dev")
```

# Introduction

ttBulk is a collection of wrapper functions for bulk tanscriptomic
analyses that follows the “tidy” paradigm. The data structure is a
tibble with columns for

  - sample identifier column
  - transcript identifier column
  - count column
  - annotation (and other info) columns

<!-- end list -->

``` r
counts = ttBulk(ttBulk::counts, sample, transcript, count)
counts_tcga = ttBulk(ttBulk::breast_tcga_mini, sample, ens, count)
counts
```

    ## # A tibble: 938,112 x 8
    ##    sample   transcript  `Cell type` count time  condition batch factor_of_inter…
    ##    <fct>    <fct>       <fct>       <dbl> <fct> <lgl>     <int> <lgl>           
    ##  1 SRR1740… DDX11L1     b_cell         17 0 d   TRUE          0 TRUE            
    ##  2 SRR1740… WASH7P      b_cell       3568 0 d   TRUE          0 TRUE            
    ##  3 SRR1740… MIR6859-1   b_cell         57 0 d   TRUE          0 TRUE            
    ##  4 SRR1740… MIR1302-2   b_cell          1 0 d   TRUE          0 TRUE            
    ##  5 SRR1740… FAM138A     b_cell          0 0 d   TRUE          0 TRUE            
    ##  6 SRR1740… OR4F5       b_cell          0 0 d   TRUE          0 TRUE            
    ##  7 SRR1740… LOC729737   b_cell       1764 0 d   TRUE          0 TRUE            
    ##  8 SRR1740… LOC1027251… b_cell         11 0 d   TRUE          0 TRUE            
    ##  9 SRR1740… MIR6859-2   b_cell         40 0 d   TRUE          0 TRUE            
    ## 10 SRR1740… OR4F29      b_cell          0 0 d   TRUE          0 TRUE            
    ## # … with 938,102 more rows

In brief you can: + Going from BAM/SAM to a tidy data frame of counts
(FeatureCounts) + Adding gene symbols from ensembl IDs + Aggregating
duplicated gene symbols + Adding scaled counts + Adding principal
components + Adding MDS components + Rotating principal component or MDS
dimensions + Running differential transcript abunance analyses (edgeR) +
Adding batch adjusted counts (Combat) + Eliminating redunant samples
and/or genes + Clustering samples and/or genes with kmeans + Adding
tissue composition (Cibersort)

# Aggregate duplicated `transcripts`

ttBulk provide the `aggregate_duplicates` function to aggregate
duplicated transcripts (e.g., isoforms, ensembl). For example, we often
have to convert ensembl symbols to gene/transcript symbol, but in doing
so we have to deal with duplicates. `aggregate_duplicates` takes a
tibble and column names (as symbols; for `sample`, `transcript` and
`count`) as arguments and returns a tibble with aggregate transcript
with the same name. All the rest of the column are appended, and factors
and boolean are appended as characters.

<div class="column-left">

TidyTranscriptomics

``` r yellow
counts.aggr = counts %>% aggregate_duplicates()
```

</div>

<div class="column-right">

Standard procedure

``` r
temp = data.frame(
    symbol = dge_list$genes$symbol,
    dge_list$counts
)
dge_list.nr <- by(temp, temp$symbol,
    function(df)
        if(length(df[1,1])>0)
            matrixStats:::colSums(as.matrix(df[,-1]))
)
dge_list.nr <- do.call("rbind", dge_list.nr)
colnames(dge_list.nr) <- colnames(dge_list)
```

</div>

<div style="clear:both;">

</div>

# Scale `counts`

We may want to calculate the scaled counts for library size (e.g., with
TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
`scale_abundance` takes a tibble, column names (as symbols; for
`sample`, `transcript` and `count`) and a method as arguments and
returns a tibble with additional columns with scaled data as `<NAME OF
COUNT COLUMN> scaled`.

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm = counts.aggr %>% scale_abundance()
```

</div>

<div class="column-right">

Standard procedure

``` r
library(edgeR)

myCPM <- cpm(count_m)
keep <- rowSums(myCPM > 0.5) >= 2
count_m.keep <- count_m[keep,]
[...]
dgList <- calcNormFactors(dgList, method="TMM")
dgList <- estimateCommonDisp(dgList)
dgList <- estimateTagwiseDisp(dgList)
norm_counts.table <- t(
    t(dgList$pseudo.counts)*
        (dgList$samples$norm.factors)
)
```

</div>

<div style="clear:both;">

</div>

We can easily plot the scaled density to check the normalisation
outcome. On the x axis we have the log scaled counts, on the y axes we
have the density, data is grouped by sample and coloured by cell type.

``` r
counts.norm %>%
    ggplot(aes(`count scaled` + 1, group=sample, color=`Cell type`)) +
    geom_density() +
    scale_x_log10() +
    my_theme
```

![](README_files/figure-gfm/plot_normalise-1.png)<!-- -->

# Filter `variable transcripts`

We may want to identify and filter variable transcripts.

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.variable = counts.norm %>% filter_variable()
```

    ## Getting the 500 most variable genes

</div>

<div class="column-right">

Standard procedure

``` r
library(edgeR)

x = norm_counts.table

s <- rowMeans((x-rowMeans(x))^2)
o <- order(s,decreasing=TRUE)
x <- x[o[1L:top],,drop=FALSE]

norm_counts.table = norm_counts.table[rownames(x)]

norm_counts.table$cell_type = ttBulk::counts[
    match(
        ttBulk::counts$sample,
        rownames(norm_counts.table)
    ),
    "Cell type"
]
```

</div>

<div style="clear:both;">

</div>

# Reduce `dimensions`

We may want to reduce the dimensions of our data, for example using PCA
or MDS algorithms. `reduce_dimensions` takes a tibble, column names (as
symbols; for `sample`, `transcript` and `count`) and a method (e.g., MDS
or PCA) as arguments and returns a tibble with additional columns for
the reduced dimensions.

**MDS** (Robinson et al., 10.1093/bioinformatics/btp616)

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.MDS =
  counts.norm %>%
  reduce_dimensions(method="MDS", .dims = 6)
```

</div>

<div class="column-right">

Standard procedure

``` r
library(limma)

count_m_log = log(count_m + 1)
cmds1_2 = count_m_log %>% plotMDS(dim.plot = 1:2, plot = FALSE)
cmds3_4 = count_m_log %>% plotMDS(dim.plot = 3:4, plot = FALSE)
cmds5_6 = count_m_log %>% plotMDS(dim.plot = 5:6, plot = FALSE)


cmds = cbind(
    data.frame(cmds1_2$x, cmds1_2$y),
    cbind(
        data.frame(cmds3_4$x, cmds3_4$y),
        data.frame(cmds5_6$x, cmds5_6$y)
    )
) %>%
    setNames(sprintf("Dim%s", 1:6))
cmds$cell_type = ttBulk::counts[
    match(ttBulk::counts$sample, rownames(cmds)),
    "Cell type"
]
```

</div>

<div style="clear:both;">

</div>

On the x and y axes axis we have the reduced dimensions 1 to 3, data is
coloured by cell
type.

``` r
counts.norm.MDS %>% select(sample, contains("Dim"), `Cell type`, time ) %>% distinct()
```

    ## # A tibble: 48 x 9
    ##    sample       Dim1   Dim2    Dim3    Dim4    Dim5  Dim6 `Cell type`      time 
    ##    <chr>       <dbl>  <dbl>   <dbl>   <dbl>   <dbl> <dbl> <chr>            <chr>
    ##  1 SRR1740034  1.84   0.686 -2.79    0.279   0.0862 0.383 b_cell           0 d  
    ##  2 SRR1740035  1.86   0.596 -2.80    0.277   0.109  0.451 b_cell           1 d  
    ##  3 SRR1740036  1.87   0.481 -2.71    0.394   0.109  0.520 b_cell           3 d  
    ##  4 SRR1740037  1.81   0.675 -2.75    0.303   0.0961 0.293 b_cell           7 d  
    ##  5 SRR1740038 -1.05  -2.19  -0.129  -0.108  -1.08   0.172 dendritic_myelo… 0 d  
    ##  6 SRR1740039 -0.989 -2.13  -0.0886 -0.0869 -0.981  0.143 dendritic_myelo… 1 d  
    ##  7 SRR1740040 -0.976 -2.35  -0.153  -0.0691 -1.25   0.221 dendritic_myelo… 3 d  
    ##  8 SRR1740041 -0.943 -2.24  -0.100  -0.0410 -1.05   0.167 dendritic_myelo… 7 d  
    ##  9 SRR1740042 -1.65  -2.28  -0.0460 -0.423   0.977  0.244 monocyte         0 d  
    ## 10 SRR1740043 -1.55  -2.04  -0.0310 -0.541   0.943  0.175 monocyte         1 d  
    ## # … with 38 more rows

``` r
counts.norm.MDS %>%
    select(contains("Dim"), sample, `Cell type`) %>%
  distinct() %>%
  GGally::ggpairs(columns = 1:6, ggplot2::aes(colour=`Cell type`))
```

![](README_files/figure-gfm/plot_mds-1.png)<!-- -->

**PCA**

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.PCA =
  counts.norm %>%
  reduce_dimensions(method="PCA", .dims = 6)
```

</div>

<div class="column-right">

Standard procedure

``` r
count_m_log = log(count_m + 1)
pc = count_m_log %>% prcomp(scale = TRUE)
variance = pc$sdev^2
variance = (variance / sum(variance))[1:6]
pc$cell_type = counts[
    match(counts$sample, rownames(pc)),
    "Cell type"
]
```

</div>

<div style="clear:both;">

</div>

On the x and y axes axis we have the reduced dimensions 1 to 3, data is
coloured by cell
type.

``` r
counts.norm.PCA %>% select(sample, contains("PC"), `Cell type`, time ) %>% distinct()
```

    ## # A tibble: 48 x 9
    ##    sample        PC1     PC2     PC3    PC4     PC5   PC6 `Cell type`      time 
    ##    <chr>       <dbl>   <dbl>   <dbl>  <dbl>   <dbl> <dbl> <chr>            <chr>
    ##  1 SRR1740034  0.130  0.137  -0.144  0.248  -0.0598 0.138 b_cell           0 d  
    ##  2 SRR1740035  0.128  0.137  -0.146  0.253  -0.0536 0.137 b_cell           1 d  
    ##  3 SRR1740036  0.130  0.135  -0.146  0.251  -0.0481 0.138 b_cell           3 d  
    ##  4 SRR1740037  0.128  0.138  -0.151  0.248  -0.0516 0.137 b_cell           7 d  
    ##  5 SRR1740038 -0.177 -0.0647 -0.176  0.0304 -0.0878 0.147 dendritic_myelo… 0 d  
    ##  6 SRR1740039 -0.176 -0.0786 -0.171  0.0328 -0.0874 0.138 dendritic_myelo… 1 d  
    ##  7 SRR1740040 -0.174 -0.0719 -0.173  0.0312 -0.103  0.142 dendritic_myelo… 3 d  
    ##  8 SRR1740041 -0.175 -0.0722 -0.179  0.0306 -0.0856 0.147 dendritic_myelo… 7 d  
    ##  9 SRR1740042 -0.198 -0.0389 -0.102  0.0234 -0.0556 0.118 monocyte         0 d  
    ## 10 SRR1740043 -0.192 -0.0567 -0.0966 0.0315 -0.0976 0.117 monocyte         1 d  
    ## # … with 38 more rows

``` r
counts.norm.PCA %>%
    select(contains("PC"), sample, `Cell type`) %>%
  distinct() %>%
  GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=`Cell type`))
```

![](README_files/figure-gfm/plot_pca-1.png)<!-- -->

**tSNE**

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.tSNE =
    counts_tcga%>%
    scale_abundance() %>%
    reduce_dimensions(
        method = "tSNE",
        perplexity=10,
        pca_scale =TRUE
    )
```

</div>

<div class="column-right">

Standard procedure

``` r
count_m_log = log(count_m + 1)

tsne = Rtsne::Rtsne(
    t(count_m_log),
    perplexity=10,
        pca_scale =TRUE
)$Y
tsne$cell_type = ttBulk::counts[
    match(ttBulk::counts$sample, rownames(tsne)),
    "Cell type"
]
```

</div>

<div style="clear:both;">

</div>

Plot

``` r
counts.norm.tSNE %>%
    select(contains("tSNE", ignore.case = FALSE), sample, Call) %>%
    distinct()
```

    ## # A tibble: 251 x 4
    ##     tSNE1   tSNE2 sample                       Call 
    ##     <dbl>   <dbl> <chr>                        <fct>
    ##  1  -5.36  -7.58  TCGA-A1-A0SD-01A-11R-A115-07 LumA 
    ##  2  10.6    4.69  TCGA-A1-A0SF-01A-11R-A144-07 LumA 
    ##  3  -7.99 -15.0   TCGA-A1-A0SG-01A-11R-A144-07 LumA 
    ##  4  -6.74  -0.591 TCGA-A1-A0SH-01A-11R-A084-07 LumA 
    ##  5  -7.07  -3.48  TCGA-A1-A0SI-01A-11R-A144-07 LumB 
    ##  6  -6.77  -6.33  TCGA-A1-A0SJ-01A-11R-A084-07 LumA 
    ##  7   7.23  32.9   TCGA-A1-A0SK-01A-12R-A084-07 Basal
    ##  8 -13.5    2.09  TCGA-A1-A0SM-01A-11R-A084-07 LumA 
    ##  9 -12.8    3.47  TCGA-A1-A0SN-01A-11R-A144-07 LumB 
    ## 10  -3.28 -23.1   TCGA-A1-A0SQ-01A-21R-A144-07 LumA 
    ## # … with 241 more rows

``` r
counts.norm.tSNE %>%
    select(contains("tSNE", ignore.case = FALSE), sample, Call) %>%
    distinct() %>%
    ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point() + my_theme
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Rotate `dimensions`

We may want to rotate the reduced dimensions (or any two numeric columns
really) of our data, of a set angle. `rotate_dimensions` takes a tibble,
column names (as symbols; for `sample`, `transcript` and `count`) and an
angle as arguments and returns a tibble with additional columns for the
rotated dimensions. The rotated dimensions will be added to the original
data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as
specified in the input arguments.

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.MDS.rotated =
  counts.norm.MDS %>%
    rotate_dimensions(`Dim1`, `Dim2`, rotation_degrees = 45)
```

</div>

<div class="column-right">

Standard procedure

``` r
rotation = function(m, d) {
    r = d * pi / 180
    ((bind_rows(
        c(`1` = cos(r), `2` = -sin(r)),
        c(`1` = sin(r), `2` = cos(r))
    ) %>% as_matrix) %*% m)
}
mds_r = pca %>% rotation(rotation_degrees)
mds_r$cell_type = counts[
    match(counts$sample, rownames(mds_r)),
    "Cell type"
]
```

</div>

<div style="clear:both;">

</div>

**Original** On the x and y axes axis we have the first two reduced
dimensions, data is coloured by cell type.

``` r
counts.norm.MDS.rotated %>%
    distinct(sample, `Dim1`,`Dim2`, `Cell type`) %>%
    ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type` )) +
  geom_point() +
  my_theme
```

![](README_files/figure-gfm/plot_rotate_1-1.png)<!-- -->

**Rotated** On the x and y axes axis we have the first two reduced
dimensions rotated of 45 degrees, data is coloured by cell type.

``` r
counts.norm.MDS.rotated %>%
    distinct(sample, `Dim1 rotated 45`,`Dim2 rotated 45`, `Cell type`) %>%
    ggplot(aes(x=`Dim1 rotated 45`, y=`Dim2 rotated 45`, color=`Cell type` )) +
  geom_point() +
  my_theme
```

![](README_files/figure-gfm/plot_rotate_2-1.png)<!-- -->

# Test `differential abundance`

We may want to test for differential transcription between sample-wise
factors of interest (e.g., with edgeR). `test_differential_abundance`
takes a tibble, column names (as symbols; for `sample`, `transcript` and
`count`) and a formula representing the desired linear model as
arguments and returns a tibble with additional columns for the
statistics from the hypothesis test (e.g., log fold change, p-value and
false discovery rate).

<div class="column-left">

TidyTranscriptomics

``` r
counts.de =
    counts %>%
    test_differential_abundance( ~ condition, action="get")
```

</div>

<div class="column-right">

Standard procedure

``` r
library(edgeR)

design =
        model.matrix(
            object = .formula,
            data = df_for_edgeR
        )

DGEList(counts = counts) %>%
        calcNormFactors(method = "TMM") %>%
        estimateGLMCommonDisp(design) %>%
        estimateGLMTagwiseDisp(design) %>%
        glmFit(design) %>%
        glmLRT(coef = 2) %>%
        topTags(n = 999999) %$%
        table
```

</div>

<div style="clear:both;">

</div>

``` r
counts.de
```

    ## # A tibble: 19,544 x 8
    ##    transcript   logFC  logCPM    LR   PValue      FDR is_de `filter out low cou…
    ##    <chr>        <dbl>   <dbl> <dbl>    <dbl>    <dbl> <lgl> <lgl>               
    ##  1 ANKRD18DP     4.83 -0.436  123.  1.43e-28 1.59e-24 TRUE  FALSE               
    ##  2 IGLL3P        5.35 -0.0631  89.6 2.94e-21 1.63e-17 TRUE  FALSE               
    ##  3 LOC101929322  2.91  1.17    54.7 1.44e-13 5.30e-10 TRUE  FALSE               
    ##  4 GRAMD1B      -5.86  4.75    43.4 4.39e-11 1.16e- 7 TRUE  FALSE               
    ##  5 BTNL9         6.34  4.40    43.1 5.21e-11 1.16e- 7 TRUE  FALSE               
    ##  6 PXDC1         3.91  1.69    40.6 1.89e-10 3.38e- 7 TRUE  FALSE               
    ##  7 MAATS1        2.98 -1.16    40.3 2.22e-10 3.38e- 7 TRUE  FALSE               
    ##  8 SRMS          4.02  0.257   39.8 2.88e-10 3.38e- 7 TRUE  FALSE               
    ##  9 IQGAP2       -3.55  7.72    39.7 3.01e-10 3.38e- 7 TRUE  FALSE               
    ## 10 NT5DC3        2.26  4.56    39.6 3.05e-10 3.38e- 7 TRUE  FALSE               
    ## # … with 19,534 more rows

The functon `test_differential_abundance` operated with contrasts too.
The constrasts hve the name of the design matrix (generally
<NAME_COLUMN_COVARIATE><VALUES_OF_COVARIATE>)

``` r
counts.de =
    counts %>%
    test_differential_abundance(
        ~ 0 + condition,                  
        .contrasts = c( "conditionTRUE - conditionFALSE"), 
        action="get"
    )
```

# Adjust `counts`

We may want to adjust `counts` for (known) unwanted variation.
`adjust_abundance` takes as arguments a tibble, column names (as
symbols; for `sample`, `transcript` and `count`) and a formula
representing the desired linear model where the first covariate is the
factor of interest and the second covariate is the unwanted variation,
and returns a tibble with additional columns for the adjusted counts as
`<COUNT COLUMN> adjusted`. At the moment just an unwanted covariated is
allowed at a time.

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.adj =
    counts.norm %>% adjust_abundance(
        ~ factor_of_interest + batch,
        action = "get"
)
```

</div>

<div class="column-right">

Standard procedure

``` r
library(sva)

count_m_log = log(count_m + 1)

design =
        model.matrix(
            object = ~ factor_of_interest + batch,
            data = annotation
        )

count_m_log.sva =
    ComBat(
            batch = design[,2],
            mod = design,
            ...
        )

count_m_log.sva = ceiling(exp(count_m_log.sva) -1)
count_m_log.sva$cell_type = counts[
    match(counts$sample, rownames(count_m_log.sva)),
    "Cell type"
]
```

</div>

<div style="clear:both;">

</div>

# Deconvolve `Cell type composition`

We may want to infer the cell type composition of our samples (with the
algorithm Cibersort; Newman et al., 10.1038/nmeth.3337).
`deconvolve_cellularity` takes as arguments a tibble, column names (as
symbols; for `sample`, `transcript` and `count`) and returns a tibble
with additional columns for the adjusted cell type proportions.

<div class="column-left">

TidyTranscriptomics

``` r
counts.cibersort =
    counts %>%
    deconvolve_cellularity(action="add", cores=2)
```

</div>

<div class="column-right">

Standard procedure

``` r
source(‘CIBERSORT.R’)
count_m %>% write.table("mixture_file.txt")
results <- CIBERSORT(
    "sig_matrix_file.txt",
    "mixture_file.txt",
    perm=100, QN=TRUE
)
results$cell_type = ttBulk::counts[
    match(ttBulk::counts$sample, rownames(results)),
    "Cell type"
]
```

</div>

<div style="clear:both;">

</div>

With the new annotated data frame, we can plot the distributions of cell
types across samples, and compare them with the nominal cell type labels
to check for the purity of isolation. On the x axis we have the cell
types inferred by Cibersort, on the y axis we have the inferred
proportions. The data is facetted and coloured by nominal cell types
(annotation given by the researcher after FACS sorting).

``` r
counts.cibersort %>%
    select(contains("cibersort:"), everything()) %>%
    gather(`Cell type inferred`, `proportion`, 1:22) %>%
  distinct(sample, `Cell type`, `Cell type inferred`, proportion) %>%
  ggplot(aes(x=`Cell type inferred`, y=proportion, fill=`Cell type`)) +
  geom_boxplot() +
  facet_wrap(~`Cell type`) +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio=1/5)
```

![](README_files/figure-gfm/plot_cibersort-1.png)<!-- -->

# Cluster `samples`

We may want to cluster our data (e.g., using k-means sample-wise).
`cluster_elements` takes as arguments a tibble, column names (as
symbols; for `sample`, `transcript` and `count`) and returns a tibble
with additional columns for the cluster annotation. At the moment only
k-means clustering is supported, the plan is to introduce more
clustering methods.

**k-means**

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.cluster = counts.norm.MDS %>%
  cluster_elements(method="kmeans", centers = 2 )
```

</div>

<div class="column-right">

Standard procedure

``` r
count_m_log = log(count_m + 1)

k = kmeans(count_m_log, iter.max = 1000, ...)
cluster = k$cluster

cluster$cell_type = ttBulk::counts[
    match(ttBulk::counts$sample, rownames(cluster)),
    c("Cell type", "Dim1", "Dim2")
]
```

</div>

<div style="clear:both;">

</div>

We can add cluster annotation to the MDS dimesion reduced data set and
plot.

``` r
 counts.norm.cluster %>%
    distinct(sample, `Dim1`, `Dim2`, `cluster kmeans`) %>%
    ggplot(aes(x=`Dim1`, y=`Dim2`, color=`cluster kmeans`)) +
  geom_point() +
  my_theme
```

![](README_files/figure-gfm/plot_cluster-1.png)<!-- -->

**SNN**

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.SNN =
    counts.norm.tSNE %>%
    cluster_elements(method = "SNN")
```

</div>

<div class="column-right">

Standard procedure

``` r
library(Seurat)

snn = CreateSeuratObject(count_m)
snn = ScaleData(
    snn, display.progress = TRUE,
    num.cores=4, do.par = TRUE
)
snn = FindVariableFeatures(snn, selection.method = "vst")
snn = FindVariableFeatures(snn, selection.method = "vst")
snn = RunPCA(snn, npcs = 30)
snn = FindNeighbors(snn)
snn = FindClusters(snn, method = "igraph", ...)
snn = snn[["seurat_clusters"]]

snn$cell_type = ttBulk::counts[
    match(ttBulk::counts$sample, rownames(snn)),
    c("Cell type", "Dim1", "Dim2")
]
```

</div>

<div style="clear:both;">

</div>

``` r
counts.norm.SNN %>%
    select(contains("tSNE", ignore.case = FALSE), `cluster SNN`, sample) %>%
    distinct()
```

    ## # A tibble: 251 x 4
    ##     tSNE1   tSNE2 `cluster SNN` sample                      
    ##     <dbl>   <dbl> <fct>         <chr>                       
    ##  1  -5.36  -7.58  0             TCGA-A1-A0SD-01A-11R-A115-07
    ##  2  10.6    4.69  2             TCGA-A1-A0SF-01A-11R-A144-07
    ##  3  -7.99 -15.0   1             TCGA-A1-A0SG-01A-11R-A144-07
    ##  4  -6.74  -0.591 0             TCGA-A1-A0SH-01A-11R-A084-07
    ##  5  -7.07  -3.48  0             TCGA-A1-A0SI-01A-11R-A144-07
    ##  6  -6.77  -6.33  1             TCGA-A1-A0SJ-01A-11R-A084-07
    ##  7   7.23  32.9   3             TCGA-A1-A0SK-01A-12R-A084-07
    ##  8 -13.5    2.09  2             TCGA-A1-A0SM-01A-11R-A084-07
    ##  9 -12.8    3.47  2             TCGA-A1-A0SN-01A-11R-A144-07
    ## 10  -3.28 -23.1   1             TCGA-A1-A0SQ-01A-21R-A144-07
    ## # … with 241 more rows

``` r
counts.norm.SNN %>%
    select(contains("tSNE", ignore.case = FALSE), `cluster SNN`, sample, Call) %>%
    gather(source, Call, c("cluster SNN", "Call")) %>%
    distinct() %>%
    ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=Call)) + geom_point() + facet_grid(~source) + my_theme
```

![](README_files/figure-gfm/SNN_plot-1.png)<!-- -->

``` r
# Do differential transcription between clusters
counts.norm.SNN %>%
    mutate(factor_of_interest = `cluster SNN` == 3) %>%
    test_differential_abundance(
    ~ factor_of_interest,
    action="get"
   )
```

    ## # A tibble: 500 x 8
    ##    ens           logFC logCPM    LR   PValue      FDR is_de `filter out low cou…
    ##    <chr>         <dbl>  <dbl> <dbl>    <dbl>    <dbl> <lgl> <lgl>               
    ##  1 ENSG00000186…  6.05   7.91  431. 7.68e-96 3.83e-93 TRUE  FALSE               
    ##  2 ENSG00000111…  2.83   9.58  385. 8.13e-86 2.03e-83 TRUE  FALSE               
    ##  3 ENSG00000181…  7.77   9.01  368. 5.83e-82 9.70e-80 TRUE  FALSE               
    ##  4 ENSG00000140…  2.58   9.50  322. 6.01e-72 7.50e-70 TRUE  FALSE               
    ##  5 ENSG00000065…  1.49  10.2   297. 1.32e-66 1.32e-64 TRUE  FALSE               
    ##  6 ENSG00000137…  3.75   8.21  295. 4.87e-66 4.05e-64 TRUE  FALSE               
    ##  7 ENSG00000124…  4.50   8.51  286. 4.21e-64 3.00e-62 TRUE  FALSE               
    ##  8 ENSG00000196…  4.76   6.92  274. 1.23e-61 7.68e-60 TRUE  FALSE               
    ##  9 ENSG00000091… -5.64  11.1   259. 3.18e-58 1.76e-56 TRUE  FALSE               
    ## 10 ENSG00000092…  2.80   8.33  248. 7.37e-56 3.68e-54 TRUE  FALSE               
    ## # … with 490 more rows

# Drop `redundant` transcripts

We may want to remove redundant elements from the original data set
(e.g., samples or transcripts), for example if we want to define
cell-type specific signatures with low sample redundancy.
`remove_redundancy` takes as arguments a tibble, column names (as
symbols; for `sample`, `transcript` and `count`) and returns a tibble
dropped recundant elements (e.g., samples). Two redundancy estimation
approaches are supported:

  - removal of highly correlated clusters of elements (keeping a
    representative) with method=“correlation”
  - removal of most proximal element pairs in a reduced dimensional
    space.

**Approach 1**

<div class="column-left">

TidyTranscriptomics

``` r
counts.norm.non_redundant =
    counts.norm.MDS %>%
  remove_redundancy(    method = "correlation" )
```

    ## Getting the 19544 most variable genes

</div>

<div class="column-right">

Standard procedure

``` r
library(widyr)

.data.correlated =
    pairwise_cor(
        counts,
        sample,
        transcript,
        rc,
        sort = TRUE,
        diag = FALSE,
        upper = FALSE
    ) %>%
    filter(correlation > correlation_threshold) %>%
    distinct(item1) %>%
    rename(!!.element := item1)

# Return non redudant data frame
counts %>% anti_join(.data.correlated) %>%
    spread(sample, rc, - transcript) %>%
    left_join(annotation)
```

</div>

<div style="clear:both;">

</div>

We can visualise how the reduced redundancy with the reduced dimentions
look like

``` r
counts.norm.non_redundant %>%
    distinct(sample, `Dim1`, `Dim2`, `Cell type`) %>%
    ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type`)) +
  geom_point() +
  my_theme
```

![](README_files/figure-gfm/plot_drop-1.png)<!-- -->

**Approach 2**

``` r
counts.norm.non_redundant =
    counts.norm.MDS %>%
  remove_redundancy(
    method = "reduced_dimensions",
    Dim_a_column = `Dim1`,
    Dim_b_column = `Dim2`
  )
```

We can visualise MDS reduced dimensions of the samples with the closest
pair removed.

``` r
counts.norm.non_redundant %>%
    distinct(sample, `Dim1`, `Dim2`, `Cell type`) %>%
    ggplot(aes(x=`Dim1`, y=`Dim2`, color=`Cell type`)) +
  geom_point() +
  my_theme
```

![](README_files/figure-gfm/plot_drop2-1.png)<!-- -->

# Other useful wrappers

The above wrapper streamline the most common processing of bulk RNA
sequencing data. Other useful wrappers are listed above.

## From BAM/SAM to tibble of gene counts

We can calculate gene counts (using FeatureCounts; Liao Y et al.,
10.1093/nar/gkz114) from a list of BAM/SAM files and format them into a
tidy structure (similar to counts).

``` r
counts = bam_sam_to_featureCounts_tibble(
    file_names,
    genome = "hg38",
    isPairedEnd = TRUE,
    requireBothEndsMapped = TRUE,
    checkFragLength = FALSE,
    useMetaFeatures = TRUE
)
```

## From ensembl IDs to gene symbol IDs

We can add gene symbols from ensembl identifiers. This is useful since
different resources use ensembl IDs while others use gene symbol IDs.

``` r
counts_ensembl %>% annotate_symbol(ens)
```

    ## # A tibble: 119 x 8
    ##    ens   iso   `read count` sample cases_0_project… cases_0_samples… transcript
    ##    <chr> <chr>        <dbl> <chr>  <chr>            <chr>            <fct>     
    ##  1 ENSG… 13             144 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ##  2 ENSG… 13              72 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ##  3 ENSG… 13               0 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ##  4 ENSG… 13            1099 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ##  5 ENSG… 13              11 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ##  6 ENSG… 13               2 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ##  7 ENSG… 13               3 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ##  8 ENSG… 13            2678 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ##  9 ENSG… 13             751 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ## 10 ENSG… 13               1 TARGE… Acute Myeloid L… Primary Blood D… TSPAN6    
    ## # … with 109 more rows, and 1 more variable: hg <fct>

# ADD versus GET modes

Every function takes this structure as input, and outputs either (i) the
new information joint to the original input data frame (default), or
(ii) just the new information, setting action=“add” or action=“get”
respectively. For example, from this data set

``` r
  counts.norm
```

    ## # A tibble: 938,112 x 13
    ##    sample transcript `Cell type` count time  condition batch factor_of_inter…
    ##    <chr>  <chr>      <chr>       <dbl> <chr> <chr>     <dbl> <chr>           
    ##  1 SRR17… A1BG       b_cell        153 0 d   TRUE          0 TRUE            
    ##  2 SRR17… A1BG-AS1   b_cell         83 0 d   TRUE          0 TRUE            
    ##  3 SRR17… A2M-AS1    b_cell          0 0 d   TRUE          0 TRUE            
    ##  4 SRR17… A2ML1      b_cell          3 0 d   TRUE          0 TRUE            
    ##  5 SRR17… A2MP1      b_cell          0 0 d   TRUE          0 TRUE            
    ##  6 SRR17… A3GALT2    b_cell          0 0 d   TRUE          0 TRUE            
    ##  7 SRR17… A4GALT     b_cell          4 0 d   TRUE          0 TRUE            
    ##  8 SRR17… A4GNT      b_cell          0 0 d   TRUE          0 TRUE            
    ##  9 SRR17… AA06       b_cell          0 0 d   TRUE          0 TRUE            
    ## 10 SRR17… AAAS       b_cell        868 0 d   TRUE          0 TRUE            
    ## # … with 938,102 more rows, and 5 more variables: `merged transcripts` <dbl>,
    ## #   `count scaled` <dbl>, TMM <dbl>, multiplier <dbl>, `filter out low
    ## #   counts` <lgl>

**action=“add”** (Default) We can add the MDS dimensions to the original
data set

``` r
  counts.norm %>%
    reduce_dimensions(
        .abundance = `count scaled`,
        method="MDS" ,
        .element = sample,
        .feature = transcript,
        .dims = 3,
        action="add"
    )
```

    ## # A tibble: 938,112 x 16
    ##    sample transcript `Cell type` count time  condition batch factor_of_inter…
    ##    <chr>  <chr>      <chr>       <dbl> <chr> <chr>     <dbl> <chr>           
    ##  1 SRR17… A1BG       b_cell        153 0 d   TRUE          0 TRUE            
    ##  2 SRR17… A1BG-AS1   b_cell         83 0 d   TRUE          0 TRUE            
    ##  3 SRR17… A2M-AS1    b_cell          0 0 d   TRUE          0 TRUE            
    ##  4 SRR17… A2ML1      b_cell          3 0 d   TRUE          0 TRUE            
    ##  5 SRR17… A2MP1      b_cell          0 0 d   TRUE          0 TRUE            
    ##  6 SRR17… A3GALT2    b_cell          0 0 d   TRUE          0 TRUE            
    ##  7 SRR17… A4GALT     b_cell          4 0 d   TRUE          0 TRUE            
    ##  8 SRR17… A4GNT      b_cell          0 0 d   TRUE          0 TRUE            
    ##  9 SRR17… AA06       b_cell          0 0 d   TRUE          0 TRUE            
    ## 10 SRR17… AAAS       b_cell        868 0 d   TRUE          0 TRUE            
    ## # … with 938,102 more rows, and 8 more variables: `merged transcripts` <dbl>,
    ## #   `count scaled` <dbl>, TMM <dbl>, multiplier <dbl>, `filter out low
    ## #   counts` <lgl>, Dim1 <dbl>, Dim2 <dbl>, Dim3 <dbl>

**action=“get”** We can get just the MDS dimensions relative to each
sample

``` r
  counts.norm %>%
    reduce_dimensions(
        .abundance = `count scaled`,
        method="MDS" ,
        .element = sample,
        .feature = transcript,
        .dims = 3,
        action="get"
    )
```

    ## # A tibble: 48 x 4
    ##    sample      Dim1   Dim2    Dim3
    ##    <chr>      <dbl>  <dbl>   <dbl>
    ##  1 SRR1740034 -2.00  0.290 -2.77  
    ##  2 SRR1740035 -1.99  0.259 -2.77  
    ##  3 SRR1740036 -1.97  0.224 -2.69  
    ##  4 SRR1740037 -1.99  0.244 -2.73  
    ##  5 SRR1740038  1.24 -2.02   0.0244
    ##  6 SRR1740039  1.16 -2.04   0.0868
    ##  7 SRR1740040  1.20 -2.01  -0.0334
    ##  8 SRR1740041  1.13 -2.08   0.0615
    ##  9 SRR1740042  1.83 -2.02   0.0885
    ## 10 SRR1740043  1.70 -1.93   0.131 
    ## # … with 38 more rows
