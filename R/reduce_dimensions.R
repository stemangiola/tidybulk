#' Dimension reduction of the transcript abundance data
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description reduce_dimensions() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and calculates the reduced dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#'
#'
#' @name reduce_dimensions
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The dimension reduction algorithm to use (PCA, MDS, tSNE).
#' @param top An integer. How many top genes to select for dimensionality reduction
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param .dims An integer. The number of dimensions your are interested in (e.g., 4 for returning the first four principal components).
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean for method="PCA", this will be passed to the `prcomp` function. It is not included in the ... argument because although the default for `prcomp` if FALSE, it is advisable to set it as TRUE.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function prcomp if you choose method="PCA" or Rtsne if you choose method="tSNE", or uwot::tumap if you choose method="umap"
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function reduces the dimensions of the transcript abundances.
#' It can use multi-dimensional scaling (MDS; DOI.org/10.1186/gb-2010-11-3-r25),
#' principal component analysis (PCA), or tSNE (Jesse Krijthe et al. 2018)
#'
#' Underlying method for PCA:
#' prcomp(scale = scale, ...)
#'
#' Underlying method for MDS:
#' limma::plotMDS(ndim = .dims, plot = FALSE, top = top)
#'
#' Underlying method for tSNE:
#' Rtsne::Rtsne(data, ...)
#'
#' Underlying method for UMAP:
#'
#'  df_source =
#' .data |>
#'
#'   # Filter NA symbol
#'   filter(!!.feature |> is.na() |> not()) |>
#'
#'   # Prepare data frame
#'   distinct(!!.feature,!!.element,!!.abundance) |>
#'
#'   # Filter most variable genes
#'   keep_variable_transcripts(top) |>
#'   reduce_dimensions(method="PCA",  .dims = calculate_for_pca_dimensions,  action="get" ) |>
#'   as_matrix(rownames = quo_name(.element)) |>
#'   uwot::tumap(...)
#'
#'
#' @return A tbl object with additional columns for the reduced dimensions
#'
#'
#' @examples
#'
#'
#'
#' counts.MDS =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'  reduce_dimensions( method="MDS", .dims = 3)
#'
#'
#' counts.PCA =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'  reduce_dimensions(method="PCA", .dims = 3)
#'
#'
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#' @export
#'
#'
setGeneric("reduce_dimensions", function(.data,
                                         .element = NULL,
                                         .feature = NULL,
                                         .abundance = NULL,
                                         method,
                                         .dims = 2,
                                         
                                         top = 500,
                                         of_samples = TRUE,
                                         transform = log1p,
                                         scale = TRUE,
                                         action = "add",
                                         ...,
                                         
                                         # DEPRECATED
                                         log_transform = NULL
                                         
)
standardGeneric("reduce_dimensions"))


.reduce_dimensions_se = function(.data,
                                 .abundance = NULL,
                                 
                                 method,
                                 .dims = 2,
                                 top = 500,
                                 of_samples = TRUE,
                                 transform = log1p,
                                 scale = TRUE,
                                 ...) {
  
  # Fix NOTEs
  . = NULL
  
  .abundance = enquo(.abundance)
  
  if(.abundance |> quo_is_symbolic()) my_assay = quo_name(.abundance)
  else my_assay = get_assay_scaled_if_exists_SE(.data)
  
  # adjust top for the max number of features I have
  if(top > nrow(.data)){
    warning(sprintf(
      "tidybulk says: the \"top\" argument %s is higher than the number of features %s",
      top,
      nrow(.data)
    ))
    
    top = min(top, nrow(.data))
  }
  
  my_assay =
    .data %>%
    
    # Filter abundant if performed
    filter_if_abundant_were_identified() %>%
    
    assay(my_assay) %>%
    
    # Filter most variable genes
    keep_variable_transcripts_SE(top = top, transform = transform) %>%
    
    # Check if log transform is needed
    transform()
  
  my_reduction_function  =
    if (tolower(method) == tolower("MDS")) {
      get_reduced_dimensions_MDS_bulk_SE
    } else if (tolower(method) == tolower("PCA")) {
      get_reduced_dimensions_PCA_bulk_SE
    } else if (tolower(method) == tolower("tSNE")) {
      get_reduced_dimensions_TSNE_bulk_SE
    } else if (tolower(method) == tolower("UMAP")) {
      get_reduced_dimensions_UMAP_bulk_SE
    } else {
      stop("tidybulk says: method must be either \"MDS\" or \"PCA\" or \"tSNE\", or \"UMAP\" ")
    }
  
  # Both dataframe and raw result object are returned
  reduced_dimensions =
    my_reduction_function(
      my_assay,
      .dims = .dims,
      top = top,
      of_samples = of_samples,
      transform = transform,
      scale=scale,
      ...
    )
  
  .data %>%
    
    # Add dimensions to metadata
    {
      .x = (.)
      if (of_samples) {
        colData(.x) = colData(.x) %>% cbind(reduced_dimensions$result)
      } else {
        rowData(.x) = rowData(.x) %>% cbind(reduced_dimensions$result)
      }
      .x
    } %>%
    
    # Add bibliography
    {
      if (tolower(method) == tolower("MDS")) {
        memorise_methods_used(., "limma")
      } else if (tolower(method) == tolower("PCA")) {
        memorise_methods_used(., "stats")
      } else if (tolower(method) == tolower("tSNE")) {
        memorise_methods_used(., "rtsne")
      } else if (tolower(method) == tolower("UMAP")) {
        memorise_methods_used(., "uwot")
      } else {
        stop("tidybulk says: method must be either \"MDS\" or \"PCA\" or \"tSNE\", or \"UMAP\" ")
      }
    } %>%
    
    # Attach edgeR for keep variable filtering
    memorise_methods_used(c("edger")) %>%
    
    # Add raw object
    attach_to_internals(reduced_dimensions$raw_result, method) %>%
    
    # Communicate the attribute added
    {
      
      rlang::inform(sprintf("tidybulk says: to access the raw results do `attr(..., \"internals\")$%s`", method), .frequency_id = sprintf("Access %s results", method),  .frequency = "always")
      
      (.)
    }
  
  
}

#' reduce_dimensions
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("reduce_dimensions",
          "SummarizedExperiment",
          .reduce_dimensions_se)

#' reduce_dimensions
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("reduce_dimensions",
          "RangedSummarizedExperiment",
          .reduce_dimensions_se)

