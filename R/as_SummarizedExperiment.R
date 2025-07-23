
#' as_SummarizedExperiment
#'
#' @description as_SummarizedExperiment() creates a `SummarizedExperiment` object from a `tbl` or `tidybulk` tbl formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#'
#' @importFrom utils data
#' @importFrom tidyr pivot_longer
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @return A `SummarizedExperiment` object
#'
#' @docType methods
#' @rdname as_SummarizedExperiment-methods
#' @export
#'
setGeneric("as_SummarizedExperiment", function(.data,
                                               
                                               
                                               .abundance = NULL)
  standardGeneric("as_SummarizedExperiment"))


.as_SummarizedExperiment = function(.data,
                                    
                                    
                                    .abundance = NULL) {
  
  # Fix NOTEs
  . = NULL
  
  # Get column names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
  .sample = col_names$.sample
  .transcript = col_names$.transcript
  .abundance = col_names$.abundance
  
  # Check if package is installed, otherwise install
  check_and_install_packages(c("SummarizedExperiment", "S4Vectors"))
  
  # If present get the scaled abundance
  if (".abundance_scaled" %in% (get_tt_columns(.data) %>% names) &&
      quo_name(get_tt_columns(.data)$.abundance_scaled) %in% colnames(.data)) {
    .abundance_scaled <- get_tt_columns(.data)$.abundance_scaled
  } else {
    .abundance_scaled <- NULL
  }
  
  # Get which columns are sample wise and which are feature wise
  col_direction = get_x_y_annotation_columns(.data,
                                             !!.sample,
                                             !!.transcript,
                                             !!.abundance,
                                             !!.abundance_scaled)
  sample_cols = col_direction$horizontal_cols
  feature_cols = col_direction$vertical_cols
  counts_cols = col_direction$counts_cols
  
  colData =
    .data %>%
    select(!!.sample, sample_cols) %>%
    distinct() %>%
    
    # Unite if multiple sample columns
    tidyr::unite(!!sample__$name, !!.sample, remove = FALSE, sep = "___") |>
    
    arrange(!!sample__$symbol) %>% {
      S4Vectors::DataFrame(
        (.) %>% select(-!!sample__$symbol),
        row.names = (.) %>% pull(!!sample__$symbol)
      )
    }
  
  rowData =
    .data %>%
    select(!!.transcript, feature_cols) %>%
    distinct() %>%
    
    # Unite if multiple sample columns
    tidyr::unite(!!feature__$name, !!.transcript, remove = FALSE, sep = "___") |>
    
    arrange(!!feature__$symbol) %>% {
      S4Vectors::DataFrame(
        (.) %>% select(-!!feature__$symbol),
        row.names = (.) %>% pull(!!feature__$symbol)
      )
    }
  
  my_assays =
    .data %>%
    
    # Unite if multiple sample columns
    tidyr::unite(!!sample__$name, !!.sample, remove = FALSE, sep = "___") |>
    
    # Unite if multiple sample columns
    tidyr::unite(!!feature__$name, !!.transcript, remove = FALSE, sep = "___") |>
    
    select(!!sample__$symbol,
           !!feature__$symbol,
           !!.abundance,
           !!.abundance_scaled,
           counts_cols) %>%
    distinct() %>%
    
    pivot_longer( cols=-c(!!feature__$symbol,!!sample__$symbol), names_to="assay", values_to= ".a") %>%
    nest(`data` = -`assay`) %>%
    mutate(`data` = `data` %>%  map(
      ~ .x %>%
        spread(!!sample__$symbol, .a) %>%
        
        # arrange sample
        select(!!feature__$symbol, rownames(colData)) |>
        
        # Arrange symbol
        arrange(!!feature__$symbol) |>
        
        # Convert
        as_matrix(rownames = feature__$name)
    ))
  
  # Build the object
  SummarizedExperiment::SummarizedExperiment(
    assays = my_assays %>% pull(`data`) %>% setNames(my_assays$assay),
    rowData = rowData,
    colData = colData
  )
  
}

#' as_SummarizedExperiment
#'
#' @export
#'
#'
#' @docType methods
#' @rdname as_SummarizedExperiment-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("as_SummarizedExperiment", "tbl_df", .as_SummarizedExperiment)