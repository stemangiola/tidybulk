#' Produces the bibliography list of your workflow
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description get_bibliography() takes as input a `tidybulk`
#'
#' @importFrom rlang enquo
#' @importFrom magrittr when
#'
#'
#' @name get_bibliography
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#'
#' @details This methods returns the bibliography list of your workflow from the internals of a tidybulk object (attr(., "internals"))
#'
#'
#' @examples
#'
#'
#' get_bibliography(tidybulk::se_mini)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return NULL. It prints a list of bibliography references for the software used through the workflow.
#' @export
#'
setGeneric("get_bibliography", function(.data)
  standardGeneric("get_bibliography"))

# Set internal
.get_bibliography = 		function(.data)
{
  
  # Fix NOTEs
  . = NULL
  
  default_methods = c("tidybulk", "tidyverse")
  
  # If there is not attributes parameter
  my_methods =
    .data |>
    when(
      !(
        !"internals" %in% (attributes(.) |> names()) &&
          !"methods_used" %in% (attr(., "internals") |> names())
      ) ~       { temp <- attr(., "internals"); temp[["methods_used"]] },
      ~ ""
    )
  
  
      my_bibliography() |>
    (\(.) .[c(default_methods, my_methods)])() |>
    unlist() |>
    writeLines()
  
}



#' get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("get_bibliography",
          "SummarizedExperiment",
          .get_bibliography)

#' get_bibliography
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("get_bibliography",
          "RangedSummarizedExperiment",
          .get_bibliography)


