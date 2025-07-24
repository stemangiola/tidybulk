#' Add differential tissue composition information to a tbl
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_differential_cellularity() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom stringr str_detect
#'
#' @name test_differential_cellularity
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula representing the desired linear model. The formula can be of two forms: multivariable (recommended) or univariable Respectively: \"factor_of_interest ~ .\" or \". ~ factor_of_interest\". The dot represents cell-type proportions, and it is mandatory. If censored regression is desired (coxph) the formula should be of the form \"survival::Surv\(y, dead\) ~ .\"
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. Either \"cibersort\", \"epic\" or \"llsr\". The regression method will be chosen based on being multivariable: lm or cox-regression (both on logit-transformed proportions); or univariable: beta or cox-regression (on logit-transformed proportions). See .formula for multi- or univariable choice.
#' @param reference A data frame. The transcript/cell_type data frame of integer transcript abundance
#' @param significance_threshold A real between 0 and 1 (usually 0.05).
#' @param ... Further parameters passed to the method deconvolve_cellularity
#'
#' @details This routine applies a deconvolution method (e.g., Cibersort; DOI: 10.1038/nmeth.3337)
#' and passes the proportions inferred into a generalised linear model (DOI:dx.doi.org/10.1007/s11749-010-0189-z)
#' or a cox regression model (ISBN: 978-1-4757-3294-8)
#'
#' Underlying method for the generalised linear model:
#' data |>
#' deconvolve_cellularity(
#' 	!!.sample, !!.transcript, !!.abundance,
#' 	method=method,
#' 	reference = reference,
#' 	action="get",
#' 	...
#' )  %>%
#' 	[..] %>%
#' 	betareg::betareg(.my_formula, .)
#'
#' Underlying method for the cox regression:
#' data |>
#' deconvolve_cellularity(
#' 	!!.sample, !!.transcript, !!.abundance,
#' 	method=method,
#' 	reference = reference,
#' 	action="get",
#' 	...
#' )  %>%
#' 	[..] %>%
#' 	mutate(.proportion_0_corrected = .proportion_0_corrected  |> boot::logit()) %>%
#' 	survival::coxph(.my_formula, .)
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#'
#'  # Regular regression
#' 	test_differential_cellularity(
#' 	 tidybulk::se_mini ,
#' 	    . ~ condition,
#' 	    cores = 1
#' 	)
#'
#' 	# Cox regression - multiple
#'
#'	tidybulk::se_mini |>
#'
#'		# Test
#'		test_differential_cellularity(
#'		    survival::Surv(days, dead) ~ .,
#'		    cores = 1
#'		)
#'
#'
#'
#' @docType methods
#' @rdname test_differential_cellularity-methods
#' @export
#'
setGeneric("test_differential_cellularity", function(.data,
                                                     .formula,
                                                     
                                                     
                                                     .abundance = NULL,
                                                     method = "cibersort",
                                                     reference = X_cibersort,
                                                     significance_threshold = 0.05,
                                                     ...)
  standardGeneric("test_differential_cellularity"))



#' Differential cellularity test with standard-error propagation
#'
#' Deconvolves bulk expression into cell–type proportions, then fits either a
#' univariable or multivariable model to compare cellularity across conditions,
#' returning estimates and standard errors.  The helper chooses the appropriate
#' engine (beta-regression, survival, or bootstrap) according to the input
#' formula.
#'
#' @param .data      A `SummarizedExperiment` (or derivative) that holds bulk
#'                   expression and sample-level metadata.
#' @param .formula   A two-sided formula specifying the outcome and covariates.
#' @param method     Character, the deconvolution method name (default
#'                   `"cibersort"`).
#' @param reference  Matrix or data frame with reference profiles.  Defaults to
#'                   the internal object `X_cibersort`.
#' @param ...        Additional arguments forwarded to
#'                   `deconvolve_cellularity()`.
#'
#' @return A tibble with one row per cell type and columns for log-fold change,
#'         standard error, confidence interval, and \(P\) value.
#'
#' @family differential-cellularity helpers
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_lgl equals when
#' @importFrom tibble as_tibble
#' @importFrom dplyr select starts_with filter pull mutate
#' @importFrom tidyr gather
#' @importFrom stringr str_split str_replace str_replace_all str_remove
#' @keywords internal
.test_differential_cellularity_se = function(.data,
                                             .formula,
                                             method = "cibersort",
                                             reference = X_cibersort,
                                             ...)
{
  
  # Fix NOTEs
  . = NULL
  
  #
  # 	if (find.package("broom", quiet = TRUE) %>% length %>% equals(0)) {
  # 		message("Installing broom needed for analyses")
  # 		install.packages("broom", repos = "https://cloud.r-project.org")
  # 	}
  
  deconvoluted =
    .data %>%
    
    # Deconvolution
    deconvolve_cellularity(
      method=method,
      prefix = sprintf("%s:", method),
      reference = reference,
      ...
    ) %>%
    colData()
  
  min_detected_proportion =
    deconvoluted %>%
    as_tibble(rownames = "sample") %>%
    select(starts_with(method)) %>%
    gather(cell_type, prop) %>%
    filter(prop > 0) %>%
    pull(prop) %>%
    min()
  
  
  # Check if test is univaiable or multivariable
  .formula %>%
    when(
      
      # Univariable
      format(.) %>%
        str_split("~") %>%
        .[[1]] %>%
        map_lgl( ~ grepl("\\. | \\.", .x)) %>%
        which	== 1 ~ {
          # Parse formula
          .my_formula =
            .formula %>%
            when(
              # If I have the dot, needed definitely for censored
              format(.) %>% grepl("\\.", .) %>% any ~ format(.) %>% str_replace("([-\\+\\*~ ]?)(\\.)", "\\1.proportion_0_corrected"),
              
              # If normal formula
              ~ sprintf(".proportion_0_corrected%s", format(.))
            ) %>%
            
            as.formula
          
          # Test
          univariable_differential_tissue_composition_SE(deconvoluted,
                                                         method,
                                                         .my_formula,
                                                         min_detected_proportion) %>%
            
            # Attach attributes
            reattach_internals(.data) %>%
            
            # Add methods used
            when(
              grepl("Surv", .my_formula) ~ (.) %>% memorise_methods_used(c("survival", "boot")),
              ~ (.) %>% memorise_methods_used("betareg")
            )
        },
      
      # Multivariable
      ~ {
        # Parse formula
        covariates =
          deconvoluted %>%
          as_tibble(rownames = "sample") %>%
          select(starts_with(method)) %>%
          colnames() %>%
          gsub(sprintf("%s:", method), "", .) %>%
          str_replace_all("[ \\(\\)]", "___")
        
        # Parse formula
        .my_formula =
          .formula %>%
          when(
            # If I have the dot, needed definitely for censored
            format(.) %>% grepl("\\.", .) %>% any ~
              format(.formula) %>%
              str_replace("([-\\+\\*~ ])(\\.)",
                          sprintf(
                            "\\1%s", paste(covariates, collapse = " + ")
                          )),
            
            # If normal formula
            ~ sprintf(".proportion_0_corrected%s", format(.))
          ) %>%
          
          as.formula
        
        # Test
        multivariable_differential_tissue_composition_SE(deconvoluted,
                                                         method,
                                                         .my_formula,
                                                         min_detected_proportion) %>%
          
          # Attach attributes
          reattach_internals(.data) %>%
          
          # Add methods used
          when(grepl("Surv", .my_formula) ~ (.) %>% memorise_methods_used(c("survival", "boot")),
               ~ (.))
        
      }) %>%
    
    # Eliminate prefix
    mutate(.cell_type = str_remove(.cell_type, sprintf("%s:", method)))
  
  
}

#' test_differential_cellularity
#'
#' @docType methods
#' @rdname test_differential_cellularity-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "test_differential_cellularity",
  "SummarizedExperiment",
  .test_differential_cellularity_se
)

#' test_differential_cellularity
#'
#' @docType methods
#' @rdname test_differential_cellularity-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "test_differential_cellularity",
  "RangedSummarizedExperiment",
  .test_differential_cellularity_se
)



#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stringr str_remove
#' @importFrom stringr str_replace_all
#'
multivariable_differential_tissue_composition_SE = function(
    deconvoluted,
    method,
    .my_formula,
    min_detected_proportion
){
  # Fix NOTEs
  . = NULL
  constrast = NULL
  results_regression =
    deconvoluted %>%
    as_tibble(rownames = "sample") %>%
    
    # Replace 0s - before
    mutate(across(starts_with(method), function(.x) if_else(.x==0, min_detected_proportion, .x))) %>%
    mutate(across(starts_with(method), boot::logit)) %>%
    
    # Rename columns - after
    setNames(
      str_remove(colnames(.), sprintf("%s:", method)) %>%
        str_replace_all("[ \\(\\)]", "___")
    ) %>%
    
    # Beta or Cox
    when(
      grepl("Surv", .my_formula) %>% any ~ {
        # Check if package is installed, otherwise install
        check_and_install_packages(c("survival", "boot"))
        
        
        (.) %>%
          survival::coxph(.my_formula, .)	%>%
          broom::tidy()
      } ,
      ~ {
        (.) %>%
          lm(.my_formula, .) %>%
          broom::tidy() %>%
          filter(term != "(Intercept)")
      }
    )
  
  # Join results
  deconvoluted %>%
    as_tibble(rownames = "sample") %>%
    pivot_longer(
      names_prefix = sprintf("%s: ", method),
      cols = starts_with(method),
      names_to = ".cell_type",
      values_to = ".proportion"
    ) %>%
    tidyr::nest(cell_type_proportions = -.cell_type) %>%
    bind_cols(
      results_regression %>%
        select(-term)
    )
}

univariable_differential_tissue_composition_SE = function(
    deconvoluted,
    method,
    .my_formula,
    min_detected_proportion
){
  # Fix NOTEs
  . = NULL
  .cell_type = NULL
  .proportion = NULL
  surv_test = NULL
  deconvoluted %>%
    as_tibble(rownames = "sample") %>%
    
    # Test
    pivot_longer(
      names_prefix = sprintf("%s: ", method),
      cols = starts_with(method),
      names_to = ".cell_type",
      values_to = ".proportion"
    ) %>%
    
    # Replace 0s
    mutate(.proportion_0_corrected = if_else(.proportion==0, min_detected_proportion, .proportion)) %>%
    
    # Test survival
    tidyr::nest(cell_type_proportions = -.cell_type) %>%
    mutate(surv_test = map(
      cell_type_proportions,
      ~ {
        if(pull(., .proportion_0_corrected) %>% unique %>% length %>%  `<=` (3)) return(NULL)
        
        # See if regression if censored or not
        .x %>%
          when(
            grepl("Surv", .my_formula) %>% any ~ {
              # Check if package is installed, otherwise install
              check_and_install_packages(c("survival", "boot"))
              
              
              (.) %>%
                mutate(.proportion_0_corrected = .proportion_0_corrected  %>% boot::logit()) %>%
                survival::coxph(.my_formula, .)	%>%
                broom::tidy() %>%
                select(-term)
            } ,
            ~ {
              # Check if package is installed, otherwise install
              check_and_install_packages("betareg")
              
              (.) %>%
                betareg::betareg(.my_formula, .) %>%
                broom::tidy() %>%
                filter(component != "precision") %>%
                pivot_wider(names_from = term, values_from = c(estimate, std.error, statistic,   p.value)) %>%
                select(-c(`std.error_(Intercept)`, `statistic_(Intercept)`, `p.value_(Intercept)`)) %>%
                select(-component)
            }
          )
      }
    )) %>%
    
    unnest(surv_test, keep_empty = TRUE)
}
