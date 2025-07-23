
#' impute transcript abundance if missing from sample-transcript pairs
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description impute_missing_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional sample-transcript pairs with imputed transcript abundance.
#'
#' @importFrom rlang enquo
#'
#'
#' @name impute_missing_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_interest + batch)
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param suffix A character string. This is added to the imputed count column names. If empty the count column are overwritten
#' @param force_scaling A boolean. In case a abundance-containing column is not scaled (columns with _scale suffix), setting force_scaling = TRUE will result in a scaling by library size, to compensating for a possible difference in sequencing depth.
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
#'
#' @docType methods
#' @rdname impute_missing_abundance-methods
#'
#' @export
#'
#'
setGeneric("impute_missing_abundance", function(.data,
                                                .formula,
                                                
                                                
                                                .abundance = NULL,
                                                suffix = "",
                                                force_scaling = FALSE)
  standardGeneric("impute_missing_abundance"))



.impute_missing_abundance_se = function(.data,
                                        .formula,
                                        
                                        
                                        .abundance  = NULL,
                                        suffix = "",
                                        force_scaling = FALSE) {
  
  # Fix NOTEs
  . = NULL
  
  .abundance = enquo(.abundance)
  
  .assay_to_impute =
    .abundance %>%
    when(
      quo_is_symbolic(.) ~ assays(.data)[quo_names(.abundance)],
      ~ assays(.data)
    )
  
  
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
    ) %>%
    
    # Add imputed to the name
    setNames(sprintf("%s%s", names(.assay_to_impute), suffix))
  
  .assays_name_to_port = names(assays(.data)) %>% setdiff(names(.assay_to_impute))
  
  assays(.data) =
    as.list(assays(.data))[.assays_name_to_port] %>%
    c(imputed_dataframe ) %>%
    
    # Add .imputed column
    c(list(.imputed =  which_NA_matrix(.assay_to_impute[[1]] ))) %>%
    
    # Make names unique
    setNames(names(.) %>% make.unique())
  
  .data %>%
    
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

