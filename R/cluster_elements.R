#' Get clusters of elements (e.g., samples or transcripts)
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description cluster_elements() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and identify clusters in the data.
#'
#' @importFrom rlang enquo quo_is_null quo_name
#' @importFrom magrittr not
#' @importFrom dplyr select pull mutate
#' @importFrom stats kmeans
#' @importFrom SummarizedExperiment assays
#'
#'
#' @name cluster_elements
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param method A character string. The cluster algorithm to use, at the moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param ... Further parameters passed to the function kmeans
#'
#' @details identifies clusters in the data, normally of samples.
#' This function returns a tibble with additional columns for the cluster annotation.
#' At the moment only k-means (DOI: 10.2307/2346830) and SNN clustering (DOI:10.1016/j.cell.2019.05.031) is supported, the plan is to introduce more clustering methods.
#'
#' Underlying method for kmeans
#' do.call(kmeans(.data, iter.max = 1000, ...)
#'
#' Underlying method for SNN
#' .data |>
#' Seurat::CreateSeuratObject() |>
#' Seurat::ScaleData(display.progress = TRUE,num.cores = 4, do.par = TRUE) |>
#' Seurat::FindVariableFeatures(selection.method = "vst") |>
#' Seurat::RunPCA(npcs = 30) |>
#' Seurat::FindNeighbors() |>
#' Seurat::FindClusters(method = "igraph", ...)
#'
#'
#' @return A tbl object with additional columns with cluster labels
#'
#'
#' @examples
#' ## Load airway dataset for examples
#'
#'   data('airway', package = 'airway')
#'   # Ensure a 'condition' column exists for examples expecting it
#'
#'     SummarizedExperiment::colData(airway)$condition <- SummarizedExperiment::colData(airway)$dex
#'
#'
#'
#' \dontrun{
#'     cluster_elements(airway, centers = 2, method="kmeans")
#' }
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, 1(14), 281-297. doi:10.1007/978-3-642-05177-7_26
#'
#' Butler, A., Hoffman, P., Smibert, P., Papalexi, E., & Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology, 36(5), 411-420. doi:10.1038/nbt.4096
#'
#' @docType methods
#' @rdname cluster_elements-methods
#' @export
#'
setGeneric("cluster_elements", function(.data,
                                method ,
                                of_samples = TRUE,
                                transform = log1p,
                                ...)
standardGeneric("cluster_elements"))

.cluster_elements_se = function(.data,
                                method ,
                                of_samples = TRUE,
                                transform = log1p,
                                ...) {
  
  # Fix NOTEs
  . = NULL
  
  my_assay =
    .data |>
    # Filter abundant if performed
    filter_if_abundant_were_identified() |>
    assays() |>
    as.list()
  my_assay = my_assay[[get_assay_scaled_if_exists_SE(.data)]]
  
  my_cluster_function  =
    if (method == "kmeans") {
      get_clusters_kmeans_bulk_SE
    } else if (method == "SNN") {
      stop("tidybulk says: Matrix package (v1.3-3) causes an error with Seurat::FindNeighbors used in this method. We are trying to solve this issue. At the moment this option in unaviable.")
    } else {
      stop("tidybulk says: the only supported methods are \"kmeans\" or \"SNN\" ")
    }
  
  my_clusters =
    my_cluster_function(
      my_assay,
      of_samples = of_samples,
      transform = transform,
      ...
    ) |>
    as.character() |>
    as.factor()
  
  my_cluster_column = paste("cluster", method, sep="_")
  
  .data |>
    
    # Add clusters to metadata
    (\(.) {
      .x = (.)
      if (of_samples) {
        colData(.x)[,my_cluster_column] = my_clusters
      } else {
        rowData(.x)[,my_cluster_column] = my_clusters
      }
      .x
    })() |>
    
    # Add bibliography
    (\(.) {
      if (method == "kmeans") {
        memorise_methods_used(., "stats")
      } else if (method == "SNN") {
        memorise_methods_used(., "seurat")
      } else {
        stop("tidybulk says: the only supported methods are \"kmeans\" or \"SNN\" ")
      }
    })()
  
}

#' cluster_elements
#'
#' @docType methods
#' @rdname cluster_elements-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("cluster_elements",
          "SummarizedExperiment",
          .cluster_elements_se)

#' cluster_elements
#'
#' @importFrom rlang inform
#'
#' @docType methods
#' @rdname cluster_elements-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("cluster_elements",
          "RangedSummarizedExperiment",
          .cluster_elements_se)

#' Get K-mean clusters to a tibble
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom stats kmeans
#' @importFrom rlang :=
#' @importFrom purrr map2_dfr
#' @importFrom magrittr "%$%"
#' @importFrom dplyr distinct count desc summarise starts_with n
#' @importFrom purrr map2
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally samples)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
#'
get_clusters_kmeans_bulk_SE <-
  function(.data,
           of_samples = TRUE,
           transform = log1p,
           ...) {
    
    # Fix NOTEs
    cluster = NULL
    seurat_clusters = NULL
    
    # Check if centers is in dots
    dots_args = rlang::dots_list(...)
    if ("centers" %in% names(dots_args) |> not())
      stop("tidybulk says: for kmeans you need to provide the \"centers\" integer argument")
    
    .data = 
      .data |>
      
      # Check if log transform is needed
      transform() 
    
    # Decide if of samples or transcripts
    if (of_samples) 
      .data = .data |> t()
    
    
    # Wrap the do.call because of the centrers check
    
    do.call(kmeans, list(x = .data, iter.max = 1000) |> c(dots_args)) %$%
      cluster
    
  }

#' Get SNN shared nearest neighbour clusters to a tibble
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom rlang :=
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally samples)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
get_clusters_SNN_bulk_SE <-
  function(.data,
           of_samples = TRUE,
           transform = log1p,
           ...) {
    
    
    # Check if package is installed, otherwise install
    check_and_install_packages(c("cluster", "Seurat", "KernSmooth"))
    
    ndims = min(c(nrow(.data), ncol(.data), 30))-1
    
    .data |>
      Seurat::CreateSeuratObject() |>
      Seurat::ScaleData(display.progress = TRUE,
                        num.cores = 4,
                        do.par = TRUE) |>
      Seurat::FindVariableFeatures(selection.method = "vst") |>
      Seurat::RunPCA(npcs = ndims) |>
      Seurat::FindNeighbors(dims = 1:ndims) |>
              Seurat::FindClusters(method = "igraph", ...) |>
        (\(.) .[["seurat_clusters"]])() %$%
      seurat_clusters
    
  }
