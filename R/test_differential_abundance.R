#' Perform differential transcription testing using edgeR quasi-likelihood (QLT), edgeR likelihood-ratio (LR), limma-voom, limma-voom-with-quality-weights or DESeq2
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_differential_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#'
#'
#' @name test_differential_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula representing the desired linear model. If there is more than one factor, they should be in the order factor of interest + additional factors.
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param contrasts This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT), "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb"
#' @param test_above_log2_fold_change A positive real value. This works for edgeR and limma_voom methods. It uses the `treat` function, which tests that the difference in abundance is bigger than this threshold rather than zero \url{https://pubmed.ncbi.nlm.nih.gov/19176553}.
#' @param scaling_method A character string. The scaling method passed to the back-end functions: edgeR and limma-voom (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile"). Setting the parameter to \"none\" will skip the compensation for sequencing-depth for the method edgeR or limma-voom.
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param prefix A character string. The prefix you would like to add to the result columns. It is useful if you want to compare several methods.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param significance_threshold DEPRECATED - A real between 0 and 1 (usually 0.05).
#' @param fill_missing_values DEPRECATED - A boolean. Whether to fill missing sample/transcript values with the median of the transcript. This is rarely needed.
#' @param .contrasts DEPRECATED - This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param ... Further arguments passed to some of the internal experimental functions. For example for glmmSeq, it is possible to pass .dispersion, and .scaling_factor column tidyeval to skip the caluclation of dispersion and scaling and use precalculated values. This is helpful is you want to calculate those quantities on many genes and do DE testing on fewer genes. .scaling_factor is the TMM value that can be obtained with tidybulk::scale_abundance.
#'
#'
#' @details This function provides the option to use edgeR \url{https://doi.org/10.1093/bioinformatics/btp616}, limma-voom \url{https://doi.org/10.1186/gb-2014-15-2-r29}, limma_voom_sample_weights \url{https://doi.org/10.1093/nar/gkv412} or  DESeq2 \url{https://doi.org/10.1186/s13059-014-0550-8} to perform the testing.
#' All methods use raw counts, irrespective of if scale_abundance or adjust_abundance have been calculated, therefore it is essential to add covariates such as batch effects (if applicable) in the formula.
#'
#' Underlying method for edgeR framework:
#'
#' 	.data |>
#'
#' 	# Filter
#'	keep_abundant(
#'			factor_of_interest = !!(as.symbol(parse_formula(.formula)[1])),
#'			minimum_counts = minimum_counts,
#'			minimum_proportion = minimum_proportion
#'		) |>
#'
#'			# Format
#'			select(!!.transcript,!!.sample,!!.abundance) |>
#'			spread(!!.sample,!!.abundance) |>
#'			as_matrix(rownames = !!.transcript) %>%
#'
#'			# edgeR
#'			edgeR::DGEList(counts = .) |>
#'			edgeR::calcNormFactors(method = scaling_method) |>
#'			edgeR::estimateDisp(design) |>
#'
#'			# Fit
#'			edgeR::glmQLFit(design) |> // or glmFit according to choice
#'			edgeR::glmQLFTest(coef = 2, contrast = my_contrasts) // or glmLRT according to choice
#'
#'
#'	Underlying method for DESeq2 framework:
#'
#'	keep_abundant(
#'			factor_of_interest = !!as.symbol(parse_formula(.formula)[[1]]),
#'			minimum_counts = minimum_counts,
#'			minimum_proportion = minimum_proportion
#'	) |>
#'
#'	# DESeq2
#'	DESeq2::DESeqDataSet(design = .formula) |>
#'	DESeq2::DESeq() |>
#'	DESeq2::results()
#'
#'
#'
#' Underlying method for glmmSeq framework:
#'
#' counts =
#' .data %>%
#'   assay(my_assay)
#'
#' # Create design matrix for dispersion, removing random effects
#' design =
#'   model.matrix(
#'     object = .formula |> lme4::nobars(),
#'     data = metadata
#'   )
#'
#' dispersion = counts |> edgeR::estimateDisp(design = design) %$% tagwise.dispersion |> setNames(rownames(counts))
#'
#'   glmmSeq( .formula,
#'            countdata = counts ,
#'            metadata =   metadata |> as.data.frame(),
#'            dispersion = dispersion,
#'            progress = TRUE,
#'            method = method |> str_remove("(?i)^glmmSeq_" ),
#'   )
#'
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#'
#'  # edgeR
#'
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#' 	test_differential_abundance( ~ condition )
#'
#' 	# The function `test_differential_abundance` operates with contrasts too
#'
#'  tidybulk::se_mini |>
#'  identify_abundant(factor_of_interest = condition) |>
#'  test_differential_abundance(
#' 	    ~ 0 + condition,
#' 	    contrasts = c( "conditionTRUE - conditionFALSE")
#'  )
#'
#'  # DESeq2 - equivalent for limma-voom
#'
#' my_se_mini = tidybulk::se_mini
#' my_se_mini$condition  = factor(my_se_mini$condition)
#'
#' # demontrating with `fitType` that you can access any arguments to DESeq()
#' my_se_mini  |>
#'    identify_abundant(factor_of_interest = condition) |>
#'        test_differential_abundance( ~ condition, method="deseq2", fitType="local")
#'
#' # testing above a log2 threshold, passes along value to lfcThreshold of results()
#' res <- my_se_mini  |>
#'    identify_abundant(factor_of_interest = condition) |>
#'         test_differential_abundance( ~ condition, method="deseq2",
#'             fitType="local",
#'             test_above_log2_fold_change=4 )
#'
#' # Use random intercept and random effect models
#'
#'  se_mini[1:50,] |>
#'   identify_abundant(factor_of_interest = condition) |>
#'   test_differential_abundance(
#'     ~ condition + (1 + condition | time),
#'     method = "glmmseq_lme4", cores = 1
#'   )
#'
#' # confirm that lfcThreshold was used
#' \dontrun{
#'     res |>
#'         mcols() |>
#'         DESeq2::DESeqResults() |>
#'         DESeq2::plotMA()
#' }
#'
#' # The function `test_differential_abundance` operates with contrasts too
#'
#'  my_se_mini |>
#'  identify_abundant() |>
#'  test_differential_abundance(
#' 	    ~ 0 + condition,
#' 	    contrasts = list(c("condition", "TRUE", "FALSE")),
#' 	    method="deseq2",
#'          fitType="local"
#'  )
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#' @export
#'
setGeneric("test_differential_abundance", function(.data,
                                                   .formula,
                                                   
                                                   
                                                   .abundance = NULL,
                                                   contrasts = NULL,
                                                   method = "edgeR_quasi_likelihood",
                                                   test_above_log2_fold_change = NULL,
                                                   scaling_method = "TMM",
                                                   omit_contrast_in_colnames = FALSE,
                                                   prefix = "",
                                                   action = "add",
                                                   ...,
                                                   
                                                   # DEPRECATED
                                                   significance_threshold = NULL,
                                                   fill_missing_values = NULL,
                                                   .contrasts = NULL
)
standardGeneric("test_differential_abundance"))



#' @importFrom rlang inform
.test_differential_abundance_se = function(.data,
                                           .formula,
                                           .abundance = NULL,
                                           contrasts = NULL,
                                           method = "edgeR_quasi_likelihood",
                                           test_above_log2_fold_change = NULL,
                                           scaling_method = "TMM",
                                           omit_contrast_in_colnames = FALSE,
                                           prefix = "",
                                           ...)
{
  
  .abundance = enquo(.abundance)
  
  # Fix NOTEs
  . = NULL
  
  # DEPRECATION OF .constrasts
  if (is_present(.contrasts) & !is.null(.contrasts)) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.7.4", "tidybulk::test_differential_abundance(.contrasts = )", details = "The argument .contrasts is now deprecated please use contrasts (without the dot).")
    
    contrasts = .contrasts
  }
  
  # Clearly state what counts are used
  rlang::inform("=====================================
tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance
or adjust_abundance have been calculated. Therefore, it is essential to add covariates
such as batch effects (if applicable) in the formula.
=====================================", .frequency_id = "All testing methods use raw counts",  .frequency = "once")
  
  
  # Test test_above_log2_fold_change
  if(!is.null(test_above_log2_fold_change) && test_above_log2_fold_change < 0)
    stop("tidybulk says: test_above_log2_fold_change should be a positive real or NULL")
  
  # Filter abundant if performed
  .data = filter_if_abundant_were_identified(.data)
  
  if(tolower(method) %in% c("edger_quasi_likelihood", "edger_likelihood_ratio", "edger_robust_likelihood_ratio"))
    my_differential_abundance =
    get_differential_transcript_abundance_bulk_SE(
      .data,
      .formula,
      .abundance = !!.abundance,
      .contrasts = contrasts,
      sample_annotation = colData(.data),
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      ...
    )
  
  else if (grepl("voom", method))
    my_differential_abundance =
    get_differential_transcript_abundance_bulk_voom_SE(
      .data,
      .formula,
      .abundance = !!.abundance,
      .contrasts = contrasts,
      sample_annotation = colData(.data),
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      ...
    )
  
  else if(tolower(method)=="deseq2")
    my_differential_abundance =
    get_differential_transcript_abundance_deseq2_SE(
      .data,
      .formula,
      .abundance = !!.abundance,
      .contrasts = contrasts,
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      ...
    )
  
  
  else if(	tolower(method) %in% c("glmmseq_lme4", "glmmseq_glmmtmb"))
    my_differential_abundance =
    get_differential_transcript_abundance_glmmSeq_SE(
      .data,
      .formula,
      .abundance = !!.abundance,
      .contrasts = contrasts,
      sample_annotation = colData(.data),
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      ...
    )
  else
    stop("tidybulk says: the only methods supported at the moment are \"edgeR_quasi_likelihood\" (i.e., QLF), \"edgeR_likelihood_ratio\" (i.e., LRT), \"limma_voom\", \"limma_voom_sample_weights\", \"DESeq2\", \"glmmseq_lme4\", \"glmmseq_glmmTMB\"")
  
  # If action is get just return the statistics
  if(action == "get") return(my_differential_abundance$result)
  
  # Add results
  rowData(.data) = rowData(.data) %>% cbind(
    
    # Parse the statistics
    my_differential_abundance$result %>%
      as_matrix(rownames = "transcript") %>%
      .[match(rownames(rowData(.data)), rownames(.)),,drop=FALSE]
  )
  
  
  .data %>%
    
    # Add bibliography
    when(
      tolower(method) ==  "edger_likelihood_ratio" ~ (.) %>% memorise_methods_used(c("edger", "edgeR_likelihood_ratio")),
      tolower(method) ==  "edger_quasi_likelihood" ~ (.) %>% memorise_methods_used(c("edger", "edgeR_quasi_likelihood")),
      tolower(method) ==  "edger_robust_likelihood_ratio" ~ (.) %>% memorise_methods_used(c("edger", "edger_robust_likelihood_ratio")),
      tolower(method) == "limma_voom" ~ (.) %>% memorise_methods_used("voom"),
      tolower(method) == "limma_voom_sample_weights" ~ (.) %>% memorise_methods_used("voom_sample_weights"),
      tolower(method) == "deseq2" ~ (.) %>% memorise_methods_used("deseq2"),
      tolower(method) %in% c("glmmseq_lme4", "glmmseq_glmmtmb") ~ (.) %>% memorise_methods_used("glmmseq"),
      ~ stop("tidybulk says: method not supported")
    ) %>%
    
    when(
      !is.null(test_above_log2_fold_change) ~ (.) %>% memorise_methods_used("treat"),
      ~ (.)
    ) %>%
    
    attach_to_internals(my_differential_abundance$result_raw, method) %>%
    
    # Communicate the attribute added
    {
      rlang::inform(sprintf("tidybulk says: to access the raw results (fitted GLM) do `attr(..., \"internals\")$%s`", method), .frequency_id = sprintf("Access DE results %s", method),  .frequency = "always")
      
      (.)
    }
  
  
  
}

#' test_differential_abundance
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "test_differential_abundance",
  "SummarizedExperiment",
  .test_differential_abundance_se
)

#' test_differential_abundance
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "test_differential_abundance",
  "RangedSummarizedExperiment",
  .test_differential_abundance_se
)


