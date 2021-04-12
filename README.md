tidybulk - part of tidyTranscriptomics
================

[![R build
status](https://github.com/stemangiola/tidybulk/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/stemangiola/tidybulk/actions)

**Tidybulk brings transcriptomics to the tidyverse!**

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

The code is released under the version 3 of the GNU General Public
License.

<img src="https://github.com/Bioconductor/BiocStickers/blob/master/tidybulk/tidybulk.png?raw=true" width="120px" />

website:
[stemangiola.github.io/tidybulk/](http://stemangiola.github.io/tidybulk/)

Please have a look also to

-   [tidySummarizedExperiment](https://github.com/stemangiola/tidySummarizedExperiment)
    for bulk data tidy representation
-   [tidySingleCellExperiment](https://github.com/stemangiola/tidySingleCellExperiment)
    for single-cell data tidy representation
-   [tidyseurat](https://github.com/stemangiola/tidyseurat) for
    single-cell data tidy representation
-   [tidyHeatmap](https://github.com/stemangiola/tidyHeatmap) for
    heatmaps produced with tidy principles analysis and manipulation
-   [nanny](https://github.com/stemangiola/nanny) for tidy high-level
    data analysis and manipulation
-   [tidygate](https://github.com/stemangiola/tidygate) for adding
    custom gate information to your tibble

<img src="https://github.com/stemangiola/tidybulk/blob/master/inst/new_SE_usage-01.png?raw=true" width="800px" align="left"/>

<br/> <br/> <br/>

## Functions/utilities available

| Function                          | Description                                                                  |
|-----------------------------------|------------------------------------------------------------------------------|
| `identify_abundant`               | identify the abundant genes                                                  |
| `aggregate_duplicates`            | Aggregate abundance and annotation of duplicated transcripts in a robust way |
| `scale_abundance`                 | Scale (normalise) abundance for RNA sequencing depth                         |
| `reduce_dimensions`               | Perform dimensionality reduction (PCA, MDS, tSNE)                            |
| `cluster_elements`                | Labels elements with cluster identity (kmeans, SNN)                          |
| `remove_redundancy`               | Filter out elements with highly correlated features                          |
| `adjust_abundance`                | Remove known unwanted variation (Combat)                                     |
| `test_differential_abundance`     | Differential transcript abundance testing (DE)                               |
| `deconvolve_cellularity`          | Estimated tissue composition (Cibersort or llsr)                             |
| `test_differential_cellularity`   | Differential cell-type abundance testing                                     |
| `test_stratification_cellularity` | Estimate Kaplan-Meier survival differences                                   |
| `keep_variable`                   | Filter for top variable features                                             |
| `keep_abundant`                   | Filter out lowly abundant transcripts                                        |
| `test_gene_enrichment`            | Gene enrichment analyses (EGSEA)                                             |
| `test_gene_overrepresentation`    | Gene enrichment on list of transcript names (no rank)                        |

| Utilities                  | Description                                                     |
|----------------------------|-----------------------------------------------------------------|
| `get_bibliography`         | Get the bibliography of your workflow                           |
| `tidybulk`                 | add tidybulk attributes to a tibble object                      |
| `tidybulk_SAM_BAM`         | Convert SAM BAM files into tidybulk tibble                      |
| `pivot_sample`             | Select sample-wise columns/information                          |
| `pivot_transcript`         | Select transcript-wise columns/information                      |
| `rotate_dimensions`        | Rotate two dimensions of a degree                               |
| `ensembl_to_symbol`        | Add gene symbol from ensembl IDs                                |
| `symbol_to_entrez`         | Add entrez ID from gene symbol                                  |
| `describe_transcript`      | Add gene description from gene symbol                           |
| `impute_missing_abundance` | Impute abundance for missing data points using sample groupings |
| `fill_missing_abundance`   | Fill abundance for missing data points using an arbitrary value |

All functions are directly compatible with `SummarizedExperiment`
object.
