
#' Fill transcript abundance if missing from sample-transcript pairs
#'
#' \lifecycle{questioning}
#'
#' @description fill_missing_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with new observations
#'
#' @importFrom rlang enquo
#'
#'
#' @name fill_missing_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT>  | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript column
#' @param .abundance The name of the transcript abundance column
#' @param fill_with A numerical abundance with which fill the missing data points
#'
#' @details This function fills the abundance of missing sample-transcript pair using the median of the sample group defined by the formula
#'
#' @return A consistent object (to the input) non-sparse abundance
#'
#'
#'
#'
#' @examples
#'
#' print("Not run for build time.")
#'
#' # tidybulk::se_mini |>  fill_missing_abundance( fill_with = 0)
#'
#'
#' @docType methods
#' @rdname fill_missing_abundance-methods
#'
#' @export
#'
#'
setGeneric("fill_missing_abundance", function(.data,
                                              .sample= NULL,
                                              .transcript= NULL,
                                              .abundance= NULL,
                                              fill_with)
  standardGeneric("fill_missing_abundance"))