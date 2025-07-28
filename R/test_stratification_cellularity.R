#' Test of stratification of biological replicates based on tissue composition, one cell-type at the time, using Kaplan-meier curves.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @description 
#' **DEPRECATED**: This function is deprecated and will be removed in a future version.
#' test_stratification_cellularity() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo quo_name
#' @importFrom stringr str_replace
#' @importFrom dplyr mutate
#'
#' @name test_stratification_cellularity
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula representing the desired linear model. The formula can be of two forms: multivariable (recommended) or univariable Respectively: \"factor_of_interest ~ .\" or \". ~ factor_of_interest\". The dot represents cell-type proportions, and it is mandatory. If censored regression is desired (coxph) the formula should be of the form \"survival::Surv\(y, dead\) ~ .\"
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. Either \"cibersort\", \"epic\" or \"llsr\". The regression method will be chosen based on being multivariable: lm or cox-regression (both on logit-transformed proportions); or univariable: beta or cox-regression (on logit-transformed proportions). See .formula for multi- or univariable choice.
#' @param reference A data frame. The transcript/cell_type data frame of integer transcript abundance
#' @param ... Further parameters passed to the method deconvolve_cellularity
#'
#' @details This routine applies a deconvolution method (e.g., Cibersort; DOI: 10.1038/nmeth.3337)
#' and passes the proportions inferred into a generalised linear model (DOI:dx.doi.org/10.1007/s11749-010-0189-z)
#' or a cox regression model (ISBN: 978-1-4757-3294-8)
#'
#'
#' Underlying method for the test:
#' data |>
#' deconvolve_cellularity(
#' 	!!.sample, !!.transcript, !!.abundance,
#' 	method=method,
#' 	reference = reference,
#' 	...
#' )  %>%
#' 	[..] |>
#' 	mutate(.high_cellularity = .proportion > median(.proportion)) |>
#' 	survival::survdiff(data = data, .my_formula)
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#' # Examples commented out due to deprecation
#' # # Example with survival analysis
#' # tidybulk::se_mini |>
#' # test_stratification_cellularity(
#' # 	survival::Surv(days, dead) ~ .,
#' # 	method = "cibersort",
#' # 	cores = 1
#' # )
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453-457. doi:10.1038/nmeth.3337
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#' @noRd
#' @export
#'
setGeneric("test_stratification_cellularity", function(.data,
                                                       .formula,
                                                       
                                                       
                                                       .abundance = NULL,
                                                       method = c("cibersort", "llsr", "epic", "mcp_counter", "quantiseq", "xcell"),
                                                       reference = X_cibersort,
                                                       ...)
  standardGeneric("test_stratification_cellularity"))




# Set internal
#' @importFrom stringr str_replace
.test_stratification_cellularity_SE = 		function(.data,
                                                 .formula,
                                                 
                                                 
                                                 .abundance = NULL,
                                                 method = c("cibersort", "llsr", "epic", "mcp_counter", "quantiseq", "xcell"),
                                                 reference = X_cibersort,
                                                 ...)
{
  
  # Deprecation warning
  lifecycle::deprecate_warn(
    when = "2.0.0",
    what = "test_stratification_cellularity()",
    details = "This function is deprecated and will be removed in a future version."
  )
  
  # Fix NOTEs
  . = NULL
  
  # Validate method parameter - ensure only one method is selected
  if (length(method) > 1) {
    stop("Multiple methods provided. Please select only one method from: ", paste(method, collapse = ", "))
  }
  
  # Validate method is one of the supported methods
  valid_methods <- c("cibersort", "llsr", "epic", "mcp_counter", "quantiseq", "xcell")
  if (!method %in% valid_methods) {
    stop(paste("Invalid method. Please choose from:", paste(valid_methods, collapse = ", ")))
  }
  
  # Validate formula
  formula_str = format(.formula)
  if(is.na(formula_str) || !grepl(" \\.|\\. ", formula_str))
    stop("tidybulk says: in the formula a dot must be present in either these forms \". ~\" or \"~ .\" with a white-space after or before respectively")
  
  deconvoluted =
    .data |>
    
    # Deconvolution
    deconvolve_cellularity(
      method=method,
      prefix = sprintf("%s:", method),
      reference = reference,
      ...
    )
  
  # Check if test is univaiable or multivariable
  result = .formula |>
    (\(.) {
      # Parse formula
      .my_formula =
        format(.) |>
        str_replace("([~ ])(\\.)", "\\1.high_cellularity") |>
        as.formula()
      # Test
      univariable_differential_tissue_stratification_SE(deconvoluted,
                                                        method,
                                                        .my_formula) %>%
        # Attach attributes
        reattach_metadata(deconvoluted) %>%
        memorise_methods_used("test_stratification_cellularity")
    })()
  
  # Eliminate prefix if .cell_type column exists
  if (".cell_type" %in% colnames(result)) {
    result = result |> mutate(.cell_type = str_remove(.cell_type, sprintf("%s:", method)))
  }
  
  result
  
}

#' test_stratification_cellularity
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_stratification_cellularity",
          "SummarizedExperiment",
          .test_stratification_cellularity_SE)

#' test_stratification_cellularity
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
setMethod("test_stratification_cellularity",
          "RangedSummarizedExperiment",
          .test_stratification_cellularity_SE)



