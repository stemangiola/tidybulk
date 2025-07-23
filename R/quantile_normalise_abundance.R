#' Normalise by quantiles the counts of transcripts/genes
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description quantile_normalise_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and Scales transcript abundance compansating for sequencing depth (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#'
#' @importFrom rlang enquo
#'
#' @importFrom stats median
#' @importFrom dplyr join_by
#'
#' @name quantile_normalise_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A character string. Either "limma_normalize_quantiles" for limma::normalizeQuantiles or "preprocesscore_normalize_quantiles_use_target" for preprocessCore::normalize.quantiles.use.target for large-scale datasets.
#' @param target_distribution A numeric vector. If NULL the target distribution will be calculated by preprocessCore. This argument only affects the "preprocesscore_normalize_quantiles_use_target" method.
#' @param action A character string between "add" (default) and "only". "add" joins the new information to the input tbl (default), "only" return a non-redundant tbl with the just new information.
#'
#'
#' @details Tranform the feature abundance across samples so to have the same quantile distribution (using preprocessCore).
#'
#' Underlying method
#'
#' If `limma_normalize_quantiles` is chosen
#'
#' .data |>limma::normalizeQuantiles()
#'
#'  If `preprocesscore_normalize_quantiles_use_target` is chosen
#'
#' .data |>
#'    preprocessCore::normalize.quantiles.use.target(
#'       target = preprocessCore::normalize.quantiles.determine.target(.data)
#'    )
#'
#'
#' @return A tbl object with additional columns with scaled data as `<NAME OF COUNT COLUMN>_scaled`
#'
#'
#' @examples
#'
#'
#'  tidybulk::se_mini |>
#'    quantile_normalise_abundance()
#'
#'
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#' @export

setGeneric("quantile_normalise_abundance", function(.data,
                                                    
                                                    
                                                    .abundance = NULL,
                                                    method = "limma_normalize_quantiles",
                                                    target_distribution = NULL,
                                                    action = "add")
  standardGeneric("quantile_normalise_abundance"))


#' @importFrom magrittr multiply_by
#' @importFrom magrittr divide_by
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom utils tail
#' @importFrom stats na.omit
#'
.quantile_normalise_abundance_se = function(.data,
                                            
                                            
                                            .abundance = NULL,
                                            method = "limma_normalize_quantiles",
                                            target_distribution = NULL,
                                            action = NULL) {
  
  # Fix NOTEs
  . = NULL
  
  # Set .abundance default if NULL
  if (is.null(.abundance)) {
    .abundance <- rlang::quo(assayNames(.data)[1])
  } else {
    .abundance <- rlang::enquo(.abundance)
  }
  
  # Set column name for value scaled
  
  # If no assay is specified take first
  my_assay = ifelse(
    quo_is_symbol(.abundance),
    quo_name(.abundance),
    .data |>
      assayNames() |>
      extract2(1)
  )
  
  # Set column name for value scaled
  value_scaled = my_assay %>% paste0(scaled_string)
  
  # Check if the matrix is empty and avoid error
  if(.data |> assay(my_assay) |> dim() |> min() == 0)
    .data_norm =
    .data |>
    assay(my_assay) |>
    list() |>
    setNames(value_scaled)
  
  else if(tolower(method) == "limma_normalize_quantiles"){
    
    # Check if package is installed, otherwise install
    check_and_install_packages("limma")
    
    
    .data_norm <-
      .data %>%
      assay(my_assay) |>
      limma::normalizeQuantiles() |>
      list() |>
      setNames(value_scaled)
    
  }
  else if(tolower(method) == "preprocesscore_normalize_quantiles_use_target"){
    
    # Check if package is installed, otherwise install
    check_and_install_packages("preprocessCore")
    
    
    .data_norm =
      .data |>
      assay(my_assay) |>
      as.matrix()
    
    if(is.null(target_distribution)) target_distribution = preprocessCore::normalize.quantiles.determine.target(.data_norm)
    
    .data_norm =
      .data_norm |>
      preprocessCore::normalize.quantiles.use.target(
        target = target_distribution
      )
    
    colnames(.data_norm) = .data |> assay(my_assay) |> colnames()
    rownames(.data_norm) = .data |> assay(my_assay) |> rownames()
    
    .data_norm =
      .data_norm |>
      list() |>
      setNames(value_scaled)
    
  } else stop("tidybulk says: the methods must be limma_normalize_quantiles or preprocesscore")
  
  # Add the assay
  assays(.data) =  assays(.data) %>% c(.data_norm)
  
  .data %>%
    
    # Add methods
    memorise_methods_used(c("quantile")) %>%
    
    # Attach column internals
    add_tt_columns(.abundance_scaled = !!(((function(x, v)	enquo(v))(x,!!as.symbol(value_scaled))) |> drop_enquo_env()) )
  
}

#' quantile_normalise_abundance
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("quantile_normalise_abundance",
          "SummarizedExperiment",
          .quantile_normalise_abundance_se)

#' quantile_normalise_abundance
#'
#' @docType methods
#' @rdname quantile_normalise_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("quantile_normalise_abundance",
          "RangedSummarizedExperiment",
          .quantile_normalise_abundance_se)

