#' Drop redundant elements (e.g., samples) for which feature (e.g., transcript/gene) abundances are correlated
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description remove_redundancy() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) for correlation method or | <DIMENSION 1> | <DIMENSION 2> | <...> | for reduced_dimensions method, and returns a consistent object (to the input) with dropped elements (e.g., samples).
#'
#' @importFrom rlang enquo quo_name
#' @importFrom dplyr filter distinct select mutate arrange
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment assays
#'
#'
#' @name remove_redundancy
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The method to use, correlation and reduced_dimensions are available. The latter eliminates one of the most proximar pairs of samples in PCA reduced dimensions.
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param top An integer. How many top genes to select for correlation based method
#' @param Dim_a_column A character string. For reduced_dimension based calculation. The column of one principal component
#' @param Dim_b_column A character string. For reduced_dimension based calculation. The column of another principal component
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function removes redundant elements from the original data set (e.g., samples or transcripts).
#' For example, if we want to define cell-type specific signatures with low sample redundancy.
#' This function returns a tibble with dropped redundant elements (e.g., samples).
#' Two redundancy estimation approaches are supported:
#' (i) removal of highly correlated clusters of elements (keeping a representative) with method="correlation";
#' (ii) removal of most proximal element pairs in a reduced dimensional space.
#'
#' Underlying method for correlation:
#' widyr::pairwise_cor(sample, transcript,count, sort = TRUE, diag = FALSE, upper = FALSE)
#'
#' Underlying custom method for reduced dimensions:
#' select_closest_pairs = function(df) {
#' 		couples <- df |> head(n = 0)
#'
#' 		while (df |> nrow() > 0) {
#' 			pair <- df |>
#' 			arrange(dist) |>
#' 			head(n = 1)
#' 			couples <- couples |> bind_rows(pair)
#' 			df <- df |>
#' 				filter(
#' 					!`sample 1` %in% (pair |> select(1:2) |> as.character()) &
#' 						!`sample 2` %in% (pair |> select(1:2) |> as.character())
#' 				)
#' 		}
#'
#' 		couples
#'
#' 	}
#'
#'
#'
#' @return A tbl object with with dropped redundant elements (e.g., samples).
#'
#' @examples
#'
#'
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'    remove_redundancy(
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.abundance =  count,
#' 	   	method = "correlation"
#' 	   	)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#' @export
#'
#'
setGeneric("remove_redundancy", function(.data,
                                         .element = NULL,
                                         .feature = NULL,
                                         .abundance = NULL,
                                         method,
                                         
                                         of_samples = TRUE,
                                         
                                         correlation_threshold = 0.9,
                                         top = Inf,
                                         transform = identity,
                                         Dim_a_column,
                                         Dim_b_column,
                                         
                                         # DEPRECATED
                                         log_transform = NULL
)
standardGeneric("remove_redundancy"))


.remove_redundancy_se = function(.data,
                                 .element = NULL,
                                 .feature = NULL,
                                 .abundance = NULL,
                                 method,
                                 of_samples = TRUE,
                                 correlation_threshold = 0.9,
                                 top = Inf,
                                 transform = identity,
                                 
                                 Dim_a_column = NULL,
                                 Dim_b_column = NULL,
                                 
                                 # DEPRECATED
                                 log_transform = NULL) {
  
  
  # Fix NOTEs
  . = NULL
  
  Dim_a_column = enquo(Dim_a_column)
  Dim_b_column = enquo(Dim_b_column)
  
  # Check if .data has more than one element
  if(
    (nrow(.data) <= 1 & of_samples == FALSE) |
    (ncol(.data) <= 1 & of_samples == TRUE)
  )
    stop("tidybulk says: You must have more than one element (trancripts if of_samples == FALSE) to perform remove_redundancy")
  
  redundant_elements =
    if (method == "correlation") {
      # Get counts
      my_assay =
        .data |>
        # Filter abundant if performed
        filter_if_abundant_were_identified() |>
        assays() |>
        as.list()
      my_assay = my_assay[[get_assay_scaled_if_exists_SE(.data)]] |>
        # Filter most variable genes
        keep_variable_transcripts_SE(top = top, transform = transform) |>
        # Check if log transform is needed
        transform()
      # Get correlated elements
      remove_redundancy_elements_through_correlation_SE(
        my_assay,
        correlation_threshold = correlation_threshold,
        of_samples = of_samples
      )
    } else if (method == "reduced_dimensions") {
      # Get dimensions
      my_dims =
        if (of_samples) {
          colData(.data)[,c(quo_name(Dim_a_column), quo_name(Dim_b_column))]
        } else {
          rowData(.data)[,c(quo_name(Dim_a_column), quo_name(Dim_b_column))]
        }
      # Get correlated elements
      remove_redundancy_elements_though_reduced_dimensions_SE(
        my_dims
      )
    } else {
      stop(
        "tidybulk says: method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
      )
    }
  
  .data |>
    (\(.) {
            # Condition on of_samples
      if (of_samples) {
        (.)[,!colnames(.) %in% redundant_elements]
      } else {
        (.)[-!rownames(.) %in% redundant_elements,]
      }
    })() |>
    
    # Add bibliography
    (\(.) {
      if (method == "correlation") {
        memorise_methods_used(., "widyr")
      } else if (method == "reduced_dimensions") {
        (.)
      } else {
        stop("tidybulk says: method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)")
      }
    })()
  
}

#' remove_redundancy
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("remove_redundancy",
          "SummarizedExperiment",
          .remove_redundancy_se)

#' remove_redundancy
#'
#' @importFrom rlang quo
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("remove_redundancy",
          "RangedSummarizedExperiment",
          .remove_redundancy_se)





#' Drop redundant elements (e.g., samples) for which feature (e.g., genes) aboundances are correlated
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom rlang :=
#' @importFrom dplyr mutate_if
#'
#' @param .data A tibble
#' @param correlation_threshold A real number between 0 and 1
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#'
#' @return A tibble with redundant elements removed
#'
#'
remove_redundancy_elements_through_correlation_SE <- function(.data,
                                                              correlation_threshold = 0.9,
                                                              of_samples = TRUE) {
  # Comply with CRAN NOTES
  . = NULL
  cluster = NULL
  
  # Check if package is installed, otherwise install
  check_and_install_packages("widyr")
  
  
  # Get the redundant data frame
  # check that there are non-NA genes for enough samples
  if (nrow(.data) == 0) {
    stop(
      "tidybulk says: In calculating correlation there is no gene that have non NA values is all samples"
    )
  } else if (nrow(.data) < 100) {
    message(
      "tidybulk says: In calculating correlation there is < 100 genes (that have non NA values) is all samples.
The correlation calculation might not be reliable"
    )
  }
  
  .data = 
    .data |>
    as_tibble(rownames="transcript") |>
    # Prepare the data frame
    gather(sample,abundance,-transcript) |>
    (\(.) {
      if (of_samples) {
        dplyr::rename(., rc = abundance,
                      element = sample,
                      feature = transcript)
      } else {
        dplyr::rename(., rc = abundance,
                      element = transcript,
                      feature = sample)
      }
    })() |>
    # Is this necessary?
    mutate_if(is.factor, as.character) |>
    # Run pairwise correlation and return a tibble
    widyr::pairwise_cor(
      element,
      feature,
      rc,
      sort = TRUE,
      diag = FALSE,
      upper = FALSE
    ) |>
    filter(correlation > correlation_threshold) |>
    distinct(item1) |>
    pull(item1)
  
}

#' Identifies the closest pairs in a MDS context and return one of them
#' 
#' @importFrom dplyr rowwise
#' @keywords internal
#' @noRd
#'
#' @importFrom stats setNames
#' @importFrom stats dist
#'
#' @param .data A tibble
#' @param of_samples A boolean
#'
#' @return A tibble with pairs dropped
#'
#'
remove_redundancy_elements_though_reduced_dimensions_SE <-
  function(.data) {
    # This function identifies the closest pairs and return one of them
    
    
    # Calculate distances
    .data |>
      dist() |>
      # Prepare matrix
      as.matrix() |>
      as_tibble(rownames = "sample a") |>
      gather(`sample b`, dist,-`sample a`) |>
      filter(`sample a` != `sample b`) |>
      # Sort the elements of the two columns to avoid eliminating all samples
      rowwise() |>
      mutate(
        `sample 1` = c(`sample a`, `sample b`) |> sort() |> (\(.) .[1])(),
        `sample 2` = c(`sample a`, `sample b`) |> sort() |> (\(.) .[2])()
      ) |>
      ungroup() |>
      select(`sample 1`, `sample 2`, dist) |>
      distinct() |>
      # Select closestpairs
      select_closest_pairs() |>
      # Select pair to keep
      pull(1)
    
  }

