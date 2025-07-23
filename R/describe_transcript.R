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
#' @importFrom magrittr %>%
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
    ))  %>%
      .[!is.na(.)], error = function(x){}) %>%
    
    # Mouse
    c(
      tryCatch(suppressMessages(AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        keys = my_transcripts,  #ensembl_symbol_mapping$transcript %>% unique,
        column = "GENENAME",
        keytype = "SYMBOL",
        multiVals = "first"
      )) %>% .[!is.na(.)], error = function(x){})
      
    ) %>%
    
    # Parse
    unlist() %>%
    #unique() %>%
    enframe(name = "transcript", value = "description") %>%
    
    # Select just one per transcript
    distinct() %>%
    group_by(transcript) %>%
    slice(1) %>%
    ungroup()
  
  rowData(.data) = rowData(.data) %>% cbind(
    tibble(transcript = rownames(!!.data)) %>%
      left_join(description_df, by = "transcript") %>%
      select(description)
  )
  
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


