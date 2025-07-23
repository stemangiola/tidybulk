#' Extract transcript-wise information
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description pivot_transcript() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with only transcript-related columns
#'
#'
#'
#' @name pivot_transcript
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .transcript The name of the transcript column
#'
#'
#' @details This functon extracts only transcript-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to tidybulk function.
#'
#' @return A `tbl` with transcript-related information
#'
#'
#'
#'
#' @examples
#'
#'
#' 	pivot_transcript(tidybulk::se_mini 	)
#'
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#' @export
#'
#'
setGeneric("pivot_transcript", function(.data
)
  standardGeneric("pivot_transcript"))

#' Extract sample-wise information
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description pivot_sample() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with only sample-related columns
#'
#'
#'
#' @name pivot_sample
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#'
#'
#' @details This functon extracts only sample-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to tidybulk function.
#'
#' @return A `tbl` with transcript-related information
#'
#'
#'
#'
#' @examples
#'
#'
#' 	pivot_sample(tidybulk::se_mini )
#'
#'
#' @docType methods
#' @rdname pivot_sample-methods
#' @export
#'
#'
setGeneric("pivot_sample", function(.data
)
  standardGeneric("pivot_sample"))



# Set internal
.pivot_sample = 		function(.data)	{
  
  colData(.data) %>%
    
    # If reserved column names are present add .x
    setNames(
      colnames(.) %>%
        str_replace("^sample$", "sample.x")
    ) %>%
    
    # Convert to tibble
    tibble::as_tibble(rownames=sample__$name)
  
  
  
  
}

#' pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
#' @return A consistent object (to the input)
setMethod("pivot_sample",
          "SummarizedExperiment",
          .pivot_sample)

#' pivot_sample
#'
#' @docType methods
#' @rdname pivot_sample-methods
#'
#' @importFrom stringr str_replace
#'
#'
#' @return A consistent object (to the input)
setMethod("pivot_sample",
          "RangedSummarizedExperiment",
          .pivot_sample)



# Set internal
.pivot_transcript = 		function(.data
)	{
  
  # Fix NOTEs
  . = NULL
  
  range_info <-
    get_special_datasets(.data) %>%
    reduce(left_join, by=feature__$name)
  
  gene_info <-
    rowData(.data) %>%
    
    # If reserved column names are present add .x
    setNames(
      colnames(.) %>%
        str_replace("^feature$", "feature.x")
    ) %>%
    
    # Convert to tibble
    tibble::as_tibble(rownames=feature__$name)
  
  gene_info %>%
    when(
      nrow(range_info) > 0 ~ (.) %>% left_join(range_info, by=feature__$name),
      ~ (.)
    )
}

#' pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
#' @return A consistent object (to the input)
setMethod("pivot_transcript",
          "SummarizedExperiment",
          .pivot_transcript)

#' pivot_transcript
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#'
#' @return A consistent object (to the input)
setMethod("pivot_transcript",
          "RangedSummarizedExperiment",
          .pivot_transcript)

