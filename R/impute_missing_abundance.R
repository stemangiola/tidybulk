
#' impute transcript abundance if missing from sample-transcript pairs
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description impute_missing_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional sample-transcript pairs with imputed transcript abundance.
#'
#' @importFrom rlang enquo quo_name
#' @importFrom dplyr select pull mutate
#' @importFrom SummarizedExperiment assays assayNames
#' @importFrom stats setNames
#'
#'
#' @name impute_missing_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_interest + batch)
#' @param suffix A character string. This is added to the imputed count column names. If empty the count column are overwritten
#' @param force_scaling A boolean. In case a abundance-containing column is not scaled (columns with _scale suffix), setting force_scaling = TRUE will result in a scaling by library size, to compensating for a possible difference in sequencing depth.
#' @param ... Further arguments.
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#'
#' @details This function imputes the abundance of missing sample-transcript pair using the median of the sample group defined by the formula
#'
#' @return A consistent object (to the input) non-sparse abundance
#'
#'
#'
#'
#' @examples
#'
#'
#' res =
#' 	impute_missing_abundance(
#' 		tidybulk::se_mini,
#' 	~ condition
#' )
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @export
#'
#'
setGeneric("impute_missing_abundance", function(.data,
                                                .formula,
                                                suffix = "",
                                                force_scaling = FALSE,
                                                ...,
                                                abundance = assayNames(.data)[1],
                                                .abundance = NULL)
  standardGeneric("impute_missing_abundance"))



.impute_missing_abundance_se = function(.data,
                                        .formula,
                                        suffix = "",
                                        force_scaling = FALSE,
                                        ...,
                                        abundance = assayNames(.data)[1],
                                        .abundance = NULL) {
  
  # Fix NOTEs
  . = NULL
  
  # Soft-deprecate .abundance, prefer abundance (character)
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "impute_missing_abundance(.abundance)", "impute_missing_abundance(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }
  
  .abundance = enquo(.abundance)
  
  .assay_to_impute =
    if (quo_is_symbolic(.abundance)) {
      assays(.data)[quo_names(.abundance)]
    } else {
      assays(.data)
    }
  
  
  # Split data by formula and impute
  imputed_dataframe =
    map2(
      
      # Capture assay names as we need to know if scaled is in the name
      as.list(.assay_to_impute), names(.assay_to_impute),
      ~ {
        
        # Pseudo-scale if not scaled
        if(!grepl("_scaled", .y) & force_scaling) {
          library_size = colSums(.x, na.rm = TRUE)
          .x = .x / library_size
        }
        else message(sprintf("tidybulk says: %s appears not to be scaled for sequencing depth (missing _scaled suffix; if you think this column is idependent of sequencing depth ignore this message), therefore the imputation can produce non meaningful results if sequencing depth for samples are highly variable. If you use force_scaling = TRUE library size will be used for eliminatig some sequencig depth effect before imputation", .y))
        
        # Log
        need_log = max(.x, na.rm=TRUE) > 50
        if(need_log) .x = log1p(.x)
        
        # Imputation
        .x = fill_NA_matrix_with_factor_colwise(
          .x,
          # I split according to the formula
          colData(.data)[,parse_formula(.formula)]
        )
        
        # Exp back
        if(need_log) .x = exp(.x)-1
        
        # Scale back if pseudoscaled
        if(!grepl("_scaled", .y) & force_scaling) .x = .x * library_size
        
        # Return
        .x
      }
    ) |>
    
    # Add imputed to the name
    setNames(sprintf("%s%s", names(.assay_to_impute), suffix))
  
  .assays_name_to_port = names(assays(.data)) |> setdiff(names(.assay_to_impute))
  
  assays(.data) =
    as.list(assays(.data))[.assays_name_to_port] |>
    c(imputed_dataframe ) |>
    
    # Add .imputed column
    c(list(.imputed =  which_NA_matrix(.assay_to_impute[[1]] )))

  # Make names unique
  assays_names_unique <- make.unique(names(assays(.data)))
  assays(.data) <- setNames(as.list(assays(.data)), assays_names_unique)
  
  .data |>
    
    # Reattach internals
    reattach_internals(.data)
  
}



#' impute_missing_abundance
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @importFrom stringr str_replace
#'
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("impute_missing_abundance",
          "SummarizedExperiment",
          .impute_missing_abundance_se)

#' impute_missing_abundance
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @importFrom stringr str_replace
#'
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("impute_missing_abundance",
          "RangedSummarizedExperiment",
          .impute_missing_abundance_se)

