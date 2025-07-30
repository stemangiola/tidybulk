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
#' library(airway)
#' data(airway)
#' se_mini <- airway[1:100, 1:5]
#' 
#' 	pivot_transcript(se_mini 	)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
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
#' library(airway)
#' data(airway)
#' se_mini <- airway[1:100, 1:5]
#' 
#' 	pivot_sample(se_mini )
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
  
  colData(.data) |>
    
    # If reserved column names are present add .x
    setNames(
      colnames(colData(.data)) |>
        str_replace("^sample$", "sample.x")
    ) |>
    
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
#' @importFrom rlang enquo quo_name
#' @importFrom dplyr select left_join
#' @importFrom SummarizedExperiment colData rowData
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
    get_special_datasets(.data) |>
    reduce(left_join, by=feature__$name)
  
  # Get rowData and fix column names
  row_data <- rowData(.data)
  col_names <- colnames(row_data)
  fixed_col_names <- str_replace(col_names, "^feature$", "feature.x")
  row_data_fixed <- setNames(row_data, fixed_col_names)
  
  gene_info <-
    row_data_fixed |>
    tibble::as_tibble(rownames=feature__$name)
  
  gene_info |>
    (function(gene_data) {
      if (nrow(range_info) > 0) {
        gene_data |> left_join(range_info, by=feature__$name)
      } else {
        gene_data
      }
    })()
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

