#' Keep variable transcripts
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description keep_variable() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#'
#'
#' @name keep_variable
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param top Integer. Number of top transcript to consider
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details At the moment this function uses edgeR \url{https://doi.org/10.1093/bioinformatics/btp616}
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#' Underlying method:
#' 	s <- rowMeans((x - rowMeans(x)) ^ 2)
#'	o <- order(s, decreasing = TRUE)
#'	x <- x[o[1L:top], , drop = FALSE]
#'	variable_trancripts = rownames(x)
#'
#'
#'
#' @examples
#'
#'
#'
#' 	keep_variable(tidybulk::se_mini, top = 500)
#'
#'
#' @docType methods
#' @rdname keep_variable-methods
#' @export
#'
setGeneric("keep_variable", function(.data,
                                     
                                     
                                     .abundance = NULL,
                                     top = 500,
                                     transform = log1p,
                                     
                                     # DEPRECATED
                                     log_transform = TRUE
)
  standardGeneric("keep_variable"))


.keep_variable_se = function(.data,
                             top = 500,
                             transform = log1p)
{
  
  # Fix NOTEs
  . = NULL
  
  
  variable_transcripts =
    .data %>%
    
    # Filter abundant if performed
    filter_if_abundant_were_identified() %>%
    
    assays() %>%
    as.list() %>%
    .[[get_assay_scaled_if_exists_SE(.data)]] %>%
    
    # Filter most variable genes
    keep_variable_transcripts_SE(top = top, transform = transform) %>%
    
    # Take gene names
    rownames()
  
  .data[variable_transcripts]
  
  
}

#' keep_variable
#'
#' @docType methods
#' @rdname keep_variable-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_variable",
          "SummarizedExperiment",
          .keep_variable_se)

#' keep_variable
#'
#' @importFrom purrr map_chr
#' @importFrom tidyr unite
#' @importFrom Matrix colSums
#'
#' @docType methods
#' @rdname keep_variable-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_variable",
          "RangedSummarizedExperiment",
          .keep_variable_se)



#' Identify variable genes for dimensionality reduction
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom Matrix rowMeans
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#' @param top An integer. How many top genes to select
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#'
#' @return A tibble filtered genes
#'
keep_variable_transcripts_SE = function(.data,
                                        top = 500,
                                        transform = log1p) {
  
  
  # Manage Inf
  top = min(top, .data %>% nrow)
  
  message(sprintf("Getting the %s most variable genes", top))
  
  x =
    .data %>%
    
    # Check if log transform is needed
    transform()
  
  
  s <- rowMeans((x - rowMeans(x, na.rm=TRUE)) ^ 2, na.rm=TRUE)
  o <- order(s, decreasing = TRUE)
  x <- x[o[1L:top], , drop = FALSE]
  variable_trancripts = rownames(x)
  
  .data[variable_trancripts,]
  
}

