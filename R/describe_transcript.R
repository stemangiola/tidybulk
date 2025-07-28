#' Get DESCRIPTION from gene SYMBOL for Human and Mouse
#'
#' @param .data A tt or tbl object.
#' @param .transcript A character. The name of the gene symbol column.
#'
#' @return A tbl
#'
#' @examples
#'
#' describe_transcript(tidybulk::se_mini)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Carlson, M. (2019). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.8.2.
#'
#' Carlson, M. (2019). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.8.2.
#'
#' @docType methods
#' @rdname describe_transcript-methods
#' @export
#'
#'
setGeneric("describe_transcript", function(.data )
  standardGeneric("describe_transcript"))




#' describe_transcript
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom tibble enframe
#' @importFrom dplyr distinct
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr left_join

#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A `SummarizedExperiment` object
#'
.describe_transcript_SE = function(.data ) {
  
  # Fix NOTEs
  . = NULL
  
  # Check if package is installed, otherwise install
  check_and_install_packages(c("org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi"))
  
  
  # .transcript = enquo(.transcript)
  
  # Transcript rownames by default
  my_transcripts = rownames(.data)
  # .transcript %>%
  # when(
  # 	quo_is_null(.) ~ rownames(.data),
  # 	~ rowData(.data)[,quo_name(.transcript)]
  # )
  
  description_df =
    # Human
    tryCatch(suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = my_transcripts,  #ensembl_symbol_mapping$transcript %>% unique,
      column = "GENENAME",
      keytype = "SYMBOL",
      multiVals = "first"
        ))  |>
    (\(.) .[!is.na(.)])(), error = function(x){}) |>
    
    # Mouse
    c(
      tryCatch(suppressMessages(AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        keys = my_transcripts,  #ensembl_symbol_mapping$transcript %>% unique,
        column = "GENENAME",
        keytype = "SYMBOL",
        multiVals = "first"
          )) |> (\(.) .[!is.na(.)])(), error = function(x){})
    
  ) |>
  
  # Parse
  unlist() |>
  #unique() |>
  enframe(name = "transcript", value = "description") |>
    
    # Select just one per transcript
    distinct() |>
    group_by(transcript) |>
    slice(1) |>
    ungroup()
  
  # Create description column for all transcripts
  all_transcripts <- rownames(.data)
  description_matched <- description_df[match(all_transcripts, description_df$transcript), "description", drop = TRUE]
  
  # Add description to rowData
  rowData(.data)$description <- description_matched
  
  .data
}

#' describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A consistent object (to the input) including additional columns for transcript symbol
setMethod("describe_transcript", "SummarizedExperiment", .describe_transcript_SE)

#' describe_transcript
#'
#' @docType methods
#' @rdname describe_transcript-methods
#'
#' @return A consistent object (to the input) including additional columns for transcript symbol
setMethod("describe_transcript", "RangedSummarizedExperiment", .describe_transcript_SE)


