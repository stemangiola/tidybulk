#' Resolve Complete Confounders of Non-Interest
#'
#' This function identifies and resolves complete confounders among specified factors of non-interest within a `SummarizedExperiment` object.
#' Complete confounders occur when the levels of one factor are entirely predictable based on the levels of another factor.
#' Such relationships can interfere with downstream analyses by introducing redundancy or collinearity.
#'
#' The function systematically examines pairs of specified factors and determines whether they are completely confounded.
#' If a pair of factors is found to be confounded, one of the factors is adjusted or removed to resolve the issue.
#' The adjusted `SummarizedExperiment` object is returned, preserving all assays and metadata except the resolved factors.
#'
#' @param se A `SummarizedExperiment` object. This object contains assay data, row data (e.g., gene annotations), and column data (e.g., sample annotations).
#' @param ... Factors of non-interest (column names from `colData(se)`) to examine for complete confounders.
#'
#' @details
#' Complete confounders of non-interest can create dependencies between variables that may bias statistical models or violate their assumptions.
#' This function systematically addresses this by:
#' 1. Creating new columns with the suffix "___altered" for each specified factor to preserve original values
#' 2. Identifying pairs of factors in the specified columns that are fully confounded
#' 3. Resolving confounding by adjusting one of the factors in the "___altered" columns
#' 
#' The function creates new columns with the "___altered" suffix to store the modified values while preserving the original data. This allows users to compare the original and adjusted values if needed.
#'
#' The resolution strategy depends on the analysis context and can be modified in the helper function
#' `resolve_complete_confounders_of_non_interest_pair_SE()`. By default, the function adjusts one of the confounded factors in the "___altered" columns.
#'
#' @return
#' A `SummarizedExperiment` object with resolved confounders. The object retains its structure, including assays and metadata,
#' but the column data (`colData`) is updated with new "___altered" columns containing the resolved factors.
#'
#' @examples
#' # Load necessary libraries
#' library(SummarizedExperiment)
#' library(dplyr)
#'
#' # Sample annotations
#' sample_annotations <- data.frame(
#'   sample_id = paste0("Sample", seq(1, 9)),
#'   factor_of_interest = c(rep("treated", 4), rep("untreated", 5)),
#'   A = c("a1", "a2", "a1", "a2", "a1", "a2", "a1", "a2", "a3"),
#'   B = c("b1", "b1", "b2", "b1", "b1", "b1", "b2", "b1", "b3"),
#'   C = c("c1", "c1", "c1", "c1", "c1", "c1", "c1", "c1", "c3"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Simulated assay data
#' assay_data <- matrix(rnorm(100 * 9), nrow = 100, ncol = 9)
#'
#' # Row data (e.g., gene annotations)
#' row_data <- data.frame(gene_id = paste0("Gene", seq_len(100)))
#'
#' # Create SummarizedExperiment object
#' se <- SummarizedExperiment(
#'   assays = list(counts = assay_data),
#'   rowData = row_data,
#'   colData = DataFrame(sample_annotations)
#' )
#'
#' # Apply the function to resolve confounders
#' se_resolved <- resolve_complete_confounders_of_non_interest(se, A, B, C)
#'
#' # View the updated column data
#' colData(se_resolved)
#'
#' @seealso
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} for creating and handling `SummarizedExperiment` objects.
#'
#' @importFrom dplyr select
#' @importFrom rlang set_names
#' @importFrom tibble as_tibble
#' @export
setGeneric("resolve_complete_confounders_of_non_interest", function(se, ...) {
  standardGeneric("resolve_complete_confounders_of_non_interest")
})

#' Get matrix from tibble
#'
#'
#'
#'
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames The column name of the input tibble that will become the rownames of the output matrix
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#'
#' tibble(.feature = "CD3G", count=1) |> as_matrix(rownames=.feature)
#'
#' @export
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  
  # Fix NOTEs
  . = NULL
  
  rownames = enquo(rownames)
  tbl %>%
    
    # Through warning if data frame is not numerical beside the rownames column (if present)
    ifelse_pipe(
      do_check &&
        tbl %>%
        # If rownames defined eliminate it from the data frame
        ifelse_pipe(!quo_is_null(rownames), ~ .x[,-1], ~ .x) |>
        dplyr::summarise_all(class) |>
        tidyr::gather(variable, class) |>
        pull(class) |>
        unique() %>%
        `%in%`(c("numeric", "integer")) |> not() |> any(),
      ~ {
        warning("tidybulk says: there are NON-numerical columns, the matrix will NOT be numerical")
        .x
      }
    ) |>
    as.data.frame() |>
    
    # Deal with rownames column if present
    ifelse_pipe(
      !quo_is_null(rownames),
      ~ .x |>
        magrittr::set_rownames(tbl |> pull(!!rownames)) |>
        select(-1)
    ) |>
    
    # Convert to matrix
    as.matrix()
}

#' Produces the bibliography list of your workflow
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description get_bibliography() takes as input a `tidybulk`
#'
#' @importFrom rlang enquo
#'
#'
#' @name get_bibliography
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#'
#' @details This methods returns the bibliography list of your workflow from the internals of a tidybulk object (attr(., "internals"))
#'
#'
#' @examples
#'
#'
#' get_bibliography(tidybulk::se_mini)
#'
#'
#'
#' @docType methods
#' @rdname get_bibliography-methods
#'
#' @return NULL. It prints a list of bibliography references for the software used through the workflow.
#' @export
#'
setGeneric("get_bibliography", function(.data)
  standardGeneric("get_bibliography"))

# Set internal
.get_bibliography = 		function(.data)
{
  
  # Fix NOTEs
  . = NULL
  
  default_methods = c("tidybulk", "tidyverse")
  
  # If there is not attributes parameter
  my_methods =
    .data %>%
    when(
      !(
        !"internals" %in% (attributes(.) |> names()) &&
          !"methods_used" %in% (attr(., "internals") |> names())
      ) ~ 	attr(., "internals") %>% .[["methods_used"]],
      ~ ""
    )
  
  
  my_bibliography() %>%
    .[c(default_methods, my_methods)] |>
    unlist() |>
    writeLines()
  
}


#' Test of stratification of biological replicates based on tissue composition, one cell-type at the time, using Kaplan-meier curves.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_stratification_cellularity() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#'
#' @importFrom stringr str_detect
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
#' 	action="get",
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
#'
#'
#'	tidybulk::se_mini |>
#'	test_stratification_cellularity(
#'		survival::Surv(days, dead) ~ .,
#'		cores = 1
#'	)
#'
#'
#'
#' @docType methods
#' @rdname test_stratification_cellularity-methods
#' @export
#'
setGeneric("test_stratification_cellularity", function(.data,
                                                       .formula,
                                                       .sample = NULL,
                                                       .transcript = NULL,
                                                       .abundance = NULL,
                                                       method = "cibersort",
                                                       reference = X_cibersort,
                                                       ...)
  standardGeneric("test_stratification_cellularity"))


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
                                                     .sample = NULL,
                                                     .transcript = NULL,
                                                     .abundance = NULL,
                                                     method = "cibersort",
                                                     reference = X_cibersort,
                                                     significance_threshold = 0.05,
                                                     ...)
  standardGeneric("test_differential_cellularity"))


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
                                                .sample = NULL,
                                                .transcript = NULL,
                                                .abundance = NULL,
                                                suffix = "",
                                                force_scaling = FALSE)
  standardGeneric("impute_missing_abundance"))

# Set internal
.impute_missing_abundance = 	function(.data,
                                      .formula,
                                      .sample = NULL,
                                      .transcript = NULL,
                                      .abundance = NULL,
                                      suffix = "",
                                      force_scaling = FALSE)
{
  
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
  
  # Get scaled abundance if present, otherwise get abundance
  .abundance_scaled = NULL
  if(
    .data |> get_tt_columns() |> is.null() |> not() &&
    ".abundance_scaled" %in% (.data |> get_tt_columns() |> names()) &&
    quo_name(.data |> get_tt_columns() %$% .abundance_scaled) %in% (.data |> colnames()) &&
    quo_name(.data |> get_tt_columns() %$% .abundance_scaled) != quo_name(.abundance)
  )
    .abundance_scaled = get_tt_columns(.data)$.abundance_scaled
  
  # Validate data frame
  if(do_validate())  validation(.data, !!.sample, !!.transcript, !!.abundance)
  
  fill_NA_using_formula(
    .data,
    .formula,
    .sample = !!.sample,
    .transcript = !!.transcript,
    .abundance = !!.abundance,
    .abundance_scaled = !!.abundance_scaled,
    suffix = suffix,
    force_scaling = force_scaling) |>
    
    # Reattach internals
    reattach_internals(.data)
  
}


#' Fill transcript abundance if missing from sample-transcript pairs
#'
#' \lifecycle{questioning}
#'
#' @description fill_missing_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with new observations
#'
#' @importFrom rlang enquo
#'
#'
#' @name fill_missing_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT>  | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript column
#' @param .abundance The name of the transcript abundance column
#' @param fill_with A numerical abundance with which fill the missing data points
#'
#' @details This function fills the abundance of missing sample-transcript pair using the median of the sample group defined by the formula
#'
#' @return A consistent object (to the input) non-sparse abundance
#'
#'
#'
#'
#' @examples
#'
#' print("Not run for build time.")
#'
#' # tidybulk::se_mini |>  fill_missing_abundance( fill_with = 0)
#'
#'
#' @docType methods
#' @rdname fill_missing_abundance-methods
#'
#' @export
#'
#'
setGeneric("fill_missing_abundance", function(.data,
                                              .sample= NULL,
                                              .transcript= NULL,
                                              .abundance= NULL,
                                              fill_with)
  standardGeneric("fill_missing_abundance"))


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
#' @param .transcript The name of the transcript column
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
#'
#'
#' 	pivot_transcript(tidybulk::se_mini 	)
#'
#'
#' @docType methods
#' @rdname pivot_transcript-methods
#' @export
#'
#'
setGeneric("pivot_transcript", function(.data,
                                        .transcript = NULL)
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
#' @param .sample The name of the sample column
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
#'
#'
#' 	pivot_sample(tidybulk::se_mini )
#'
#'
#' @docType methods
#' @rdname pivot_sample-methods
#' @export
#'
#'
setGeneric("pivot_sample", function(.data,
                                    .sample = NULL)
  standardGeneric("pivot_sample"))


#' analyse gene rank with GSEA
#'
#' \lifecycle{maturing}
#'
#' @description test_gene_rank() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_missing
#'
#'
#' @name test_gene_rank
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .arrange_desc A column name of the column to arrange in decreasing order
#' @param species A character. For example, human or mouse. MSigDB uses the latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")
#' @param gene_sets A character vector or a list. It can take one or more of the following built-in collections as a character vector: c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"), to be used with EGSEA buildIdx. c1 is human specific. Alternatively, a list of user-supplied gene sets can be provided, to be used with EGSEA buildCustomIdx. In that case, each gene set is a character vector of Entrez IDs and the names of the list are the gene set names.
#'
#' @param gene_set DEPRECATED. Use gene_sets instead.
#'
#' @details This wrapper execute gene enrichment analyses of the dataset using a list of transcripts and GSEA.
#' This wrapper uses clusterProfiler (DOI: doi.org/10.1089/omi.2011.0118) on the back-end.
#'
#' Undelying method:
## Get gene sets signatures
#' msigdbr::msigdbr(species = species) %>%
#'
#'	# Filter specific gene_sets  if specified. This was introduced to speed up examples executionS
#'	when(
#'		!is.null(gene_sets ) ~ filter(., gs_collection %in% gene_sets ),
#'		~ (.)
#'	) |>
#'
#'	# Execute calculation
#'	nest(data = -gs_collection) |>
#'	mutate(fit =
#'				 	map(
#'				 		data,
#'				 		~ 	clusterProfiler::GSEA(
#'				 			my_entrez_rank,
#'				 			TERM2GENE=.x |> select(gs_name, ncbi_gene),
#'				 			pvalueCutoff = 1
#'				 		)
#'
#'				 	))
#'
#' @return A consistent object (to the input)
#'
#'
#'
#'
#' @examples
#'
#' print("Not run for build time.")
#'
#' \dontrun{
#'
#' df_entrez = tidybulk::se_mini
#' df_entrez = mutate(df_entrez, do_test = .feature %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
#' df_entrez  = df_entrez |> test_differential_abundance(~ condition)
#'
#'
#'	test_gene_rank(
#'		df_entrez,
#' 		.sample = .sample,
#'		.entrez = entrez,
#' 		species="Homo sapiens",
#'    gene_sets =c("C2"),
#'  .arrange_desc = logFC
#' 	)
#' }
#'
#' @docType methods
#' @rdname test_gene_rank-methods
#' @export
#'
#'
setGeneric("test_gene_rank", function(.data,
                                      .entrez,
                                      .arrange_desc,
                                      species,
                                      .sample = NULL,
                                      gene_sets  = NULL,
                                      gene_set = NULL  # DEPRECATED
)
  standardGeneric("test_gene_rank"))


#' analyse gene over-representation with GSEA
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_gene_overrepresentation() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` with the GSEA statistics
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_missing
#'
#'
#' @name test_gene_overrepresentation
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .do_test A boolean column name symbol. It indicates the transcript to check
#' @param species A character. For example, human or mouse. MSigDB uses the latin species names (e.g., \"Mus musculus\", \"Homo sapiens\")
#' @param gene_sets  A character vector. The subset of MSigDB datasets you want to test against (e.g. \"C2\"). If NULL all gene sets are used (suggested). This argument was added to avoid time overflow of the examples.
#'
#' @param gene_set DEPRECATED. Use gene_sets instead.
#'
#' @details This wrapper execute gene enrichment analyses of the dataset using a list of transcripts and GSEA.
#' This wrapper uses clusterProfiler (DOI: doi.org/10.1089/omi.2011.0118) on the back-end.
#'
#' # Get MSigDB data
#' msigdb_data = msigdbr::msigdbr(species = species)
#' 
#' # Filter for specific gene collections if provided
#' if (!is.null(gene_collections)) {
#'   msigdb_data = filter(msigdb_data, gs_collection %in% gene_collections)
#' }
#' 
#' # Process the data
#' msigdb_data |>
#'   nest(data = -gs_collection) |>
#'   mutate(test =
#'            map(
#'              data,
#'              ~ clusterProfiler::enricher(
#'                my_entrez_rank,
#'                TERM2GENE=.x |> select(gs_name, ncbi_gene),
#'                pvalueCutoff = 1
#'              ) |>
#'                as_tibble()
#'            ))
#'
#' @return A consistent object (to the input)
#'
#'
#'
#'
#' @examples
#'
#' print("Not run for build time.")
#'
#' #se_mini = aggregate_duplicates(tidybulk::se_mini, .transcript = entrez)
#' #df_entrez = mutate(df_entrez, do_test = feature %in% c("TNFRSF4", "PLCH2", "PADI4", "PAX7"))
#'
#' \dontrun{
#' 	test_gene_overrepresentation(
#' 		df_entrez,
#' 		.sample = sample,
#' 		.entrez = entrez,
#' 		.do_test = do_test,
#' 		species="Homo sapiens",
#'    gene_sets =c("C2")
#' 	)
#' }
#'
#' @docType methods
#' @rdname test_gene_overrepresentation-methods
#' @export
#'
#'
setGeneric("test_gene_overrepresentation", function(.data,
                                                    .entrez,
                                                    .do_test,
                                                    species,
                                                    .sample = NULL,
                                                    gene_sets  = NULL,
                                                    gene_set = NULL # DEPRECATED
)
  standardGeneric("test_gene_overrepresentation"))


#' analyse gene enrichment with EGSEA
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_gene_enrichment() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` of gene set information
#'
#' @importFrom rlang enquo
#'
#'
#' @name test_gene_enrichment
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .abundance The name of the transcript/gene abundance column
#' @param contrasts This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param methods A character vector. One or 3 or more methods to use in the testing (currently EGSEA errors if 2 are used). Type EGSEA::egsea.base() to see the supported GSE methods.
#' @param gene_sets A character vector or a list. It can take one or more of the following built-in collections as a character vector: c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"), to be used with EGSEA buildIdx. c1 is human specific. Alternatively, a list of user-supplied gene sets can be provided, to be used with EGSEA buildCustomIdx. In that case, each gene set is a character vector of Entrez IDs and the names of the list are the gene set names.
#' @param species A character. It can be human, mouse or rat.
#' @param cores An integer. The number of cores available
#'
#' @param method DEPRECATED. Please use methods.
#' @param .contrasts DEPRECATED - This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#'
#' @details This wrapper executes ensemble gene enrichment analyses of the dataset using EGSEA (DOI:0.12688/f1000research.12544.1)
#'
#'
#' dge =
#' 	data |>
#' 	keep_abundant(
#' 		factor_of_interest = !!as.symbol(parse_formula(.formula)[[1]]),
#' 		!!.sample, !!.entrez, !!.abundance
#' 	) %>%
#'
#' 	# Make sure transcript names are adjacent
#' 	[...] %>%
#' 	as_matrix(rownames = !!.entrez) %>%
#' 	edgeR::DGEList(counts = .)
#'
#' idx =  buildIdx(entrezIDs = rownames(dge), species = species, msigdb.gsets = msigdb.gsets,
#'	               kegg.exclude = kegg.exclude)
#'
#' dge |>
#'
#' 	# Calculate weights
#' 	limma::voom(design, plot = FALSE) |>
#'
#' 	# Execute EGSEA
#' 	egsea(
#' 		contrasts = my_contrasts,
#' 		baseGSEAs = methods,
#' 		gs.annots = idx,
#' 		sort.by = "med.rank",
#' 		num.threads = cores,
#' 		report = FALSE
#' 	)
#'
#' @return A consistent object (to the input)
#'
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' library(SummarizedExperiment)
#' se = tidybulk::se_mini
#' rowData( se)$entrez = rownames(se )
#' df_entrez = aggregate_duplicates(se,.transcript = entrez )
#'
#' library("EGSEA")
#'
#' 	test_gene_enrichment(
#'			df_entrez,
#'			~ condition,
#'			.sample = sample,
#'			.entrez = entrez,
#'			.abundance = count,
#'          methods = c("roast" , "safe", "gage"  ,  "padog" , "globaltest", "ora" ),
#'          gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
#'			species="human",
#'			cores = 2
#'		)
#'
#'}
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#' @export
#'
#'
setGeneric("test_gene_enrichment", function(.data,
                                            .formula,
                                            .sample = NULL,
                                            .entrez,
                                            .abundance = NULL,
                                            contrasts = NULL,
                                            methods = c("camera" ,    "roast" ,     "safe",       "gage"  ,     "padog" ,     "globaltest",  "ora" ),
                                            gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
                                            species,
                                            cores = 10,
                                            
                                            # DEPRECATED
                                            method = NULL,
                                            .contrasts = NULL
)
standardGeneric("test_gene_enrichment"))


#' Filter to keep only abundant transcripts/genes
#'
#' \lifecycle{questioning}
#'
#' @description 
#' Filters the data to keep only transcripts/genes that are consistently expressed above 
#' a threshold across samples. This is a filtering version of identify_abundant() that 
#' removes low-abundance features instead of just marking them.
#'
#' @param .data A `tbl` or `SummarizedExperiment` object containing transcript/gene abundance data
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param factor_of_interest The name of the column containing groups/conditions for filtering. 
#'        Used by edgeR's filterByExpr to define sample groups.
#' @param design A design matrix for more complex experimental designs. If provided, 
#'        this is passed to filterByExpr instead of factor_of_interest.
#' @param minimum_counts A positive number specifying the minimum counts per million (CPM) threshold
#'        for a transcript to be kept (default = 10)
#' @param minimum_proportion A number between 0 and 1 specifying the minimum proportion of samples
#'        that must exceed the minimum_counts threshold (default = 0.7)
#'
#' @details 
#' This function uses edgeR's filterByExpr() function to identify and keep consistently expressed features.
#' A feature is kept if it has CPM > minimum_counts in at least minimum_proportion of samples
#' in at least one experimental group (defined by factor_of_interest or design).
#' 
#' This function is similar to identify_abundant() but instead of adding an .abundant column,
#' it filters out the low-abundance features directly.
#'
#' @return 
#' Returns a filtered version of the input object containing only the features that passed
#' the abundance threshold criteria.
#'
#' @examples
#' # Basic usage
#' se_mini |> keep_abundant()
#'
#' # With custom thresholds
#' se_mini |> keep_abundant(
#'   minimum_counts = 5,
#'   minimum_proportion = 0.5
#' )
#'
#' # Using a factor of interest
#' se_mini |> keep_abundant(factor_of_interest = condition)
#'
#' @references
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of 
#' multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 
#' 40(10), 4288-4297. DOI: 10.1093/bioinformatics/btp616
#'
#' @importFrom rlang enquo
#' @importFrom dplyr filter
#'
#' @docType methods
#' @rdname keep_abundant-methods
#' @export
setGeneric("keep_abundant", function(.data,
                                     .sample = NULL,
                                     .transcript = NULL,
                                     .abundance = NULL,
                                     factor_of_interest = NULL,
                                     design = NULL,
                                     minimum_counts = 10,
                                     minimum_proportion = 0.7)
  standardGeneric("keep_abundant"))

#' Identify abundant transcripts/genes
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description 
#' Identifies transcripts/genes that are consistently expressed above a threshold across samples. 
#' This function adds a logical column `.abundant` to indicate which features pass the filtering criteria.
#'
#' @param .data A `tbl` or `SummarizedExperiment` object containing transcript/gene abundance data
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param factor_of_interest The name of the column containing groups/conditions for filtering. 
#'        Used by edgeR's filterByExpr to define sample groups.
#' @param design A design matrix for more complex experimental designs. If provided, 
#'        this is passed to filterByExpr instead of factor_of_interest.
#' @param minimum_counts A positive number specifying the minimum counts per million (CPM) threshold
#'        for a transcript to be considered abundant (default = 10)
#' @param minimum_proportion A number between 0 and 1 specifying the minimum proportion of samples
#'        that must exceed the minimum_counts threshold (default = 0.7)
#'
#' @details 
#' This function uses edgeR's filterByExpr() function to identify consistently expressed features.
#' A feature is considered abundant if it has CPM > minimum_counts in at least minimum_proportion 
#' of samples in at least one experimental group (defined by factor_of_interest or design).
#'
#' @return 
#' Returns the input object with an additional logical column `.abundant` indicating which 
#' features passed the abundance threshold criteria.
#'
#' @examples
#' # Basic usage
#' se_mini |> identify_abundant()
#'
#' # With custom thresholds
#' se_mini |> identify_abundant(
#'   minimum_counts = 5,
#'   minimum_proportion = 0.5
#' )
#'
#' # Using a factor of interest
#' se_mini |> identify_abundant(factor_of_interest = condition)
#'
#' @references
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of 
#' multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 
#' 40(10), 4288-4297. DOI: 10.1093/bioinformatics/btp616
#'
#' @importFrom rlang enquo
#' @importFrom dplyr filter
#' @importFrom tidyr drop_na
#'
#' @docType methods
#' @rdname identify_abundant-methods
#' @export
setGeneric("identify_abundant", function(.data,
                                         .sample = NULL,
                                         .transcript = NULL,
                                         .abundance = NULL,
                                         factor_of_interest = NULL,
                                         design = NULL,
                                         minimum_counts = 10,
                                         minimum_proportion = 0.7)
  standardGeneric("identify_abundant"))


#' Keep variable transcripts
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description keep_variable() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#'
#'
#' @name keep_variable
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param top Integer. Number of top transcript to consider
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details At the moment this function uses edgeR \url{https://doi.org/10.1093/bioinformatics/btp616}
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#' Underlying method:
#' 	s <- rowMeans((x - rowMeans(x)) ^ 2)
#'	o <- order(s, decreasing = TRUE)
#'	x <- x[o[1L:top], , drop = FALSE]
#'	variable_trancripts = rownames(x)
#'
#'
#'
#' @examples
#'
#'
#'
#' 	keep_variable(tidybulk::se_mini, top = 500)
#'
#'
#' @docType methods
#' @rdname keep_variable-methods
#' @export
#'
setGeneric("keep_variable", function(.data,
                                     .sample = NULL,
                                     .transcript = NULL,
                                     .abundance = NULL,
                                     top = 500,
                                     transform = log1p,
                                     
                                     # DEPRECATED
                                     log_transform = TRUE
)
  standardGeneric("keep_variable"))


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
                                                   .sample = NULL,
                                                   .transcript = NULL,
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
setGeneric("describe_transcript", function(.data,
                                           .transcript = NULL)
  standardGeneric("describe_transcript"))


#' Get cell type proportions from samples
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description deconvolve_cellularity() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with the estimated cell type abundance for each sample
#'
#' @importFrom rlang enquo
#'
#'
#' @name deconvolve_cellularity
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param reference A data frame. The methods cibersort and llsr can accept a custom rectangular dataframe with genes as rows names, cell types as column names and gene-transcript abundance as values. For exampler tidybulk::X_cibersort. The transcript/cell_type data frame of integer transcript abundance. If NULL, the default reference for each algorithm will be used. For llsr will be LM22.
#' @param method A character string. The method to be used. At the moment Cibersort (default, can accept custom reference), epic (can accept custom reference) and llsr (linear least squares regression, can accept custom reference), mcp_counter, quantiseq, xcell are available.
#' @param prefix A character string. The prefix you would like to add to the result columns. It is useful if you want to reshape data.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function Cibersort
#'
#' @details This function infers the cell type composition of our samples
#' (with the algorithm Cibersort; Newman et al., 10.1038/nmeth.3337).
#'
#' Underlying method:
#' CIBERSORT(Y = data, X = reference, ...)
#'
#' @return A consistent object (to the input) including additional columns for each cell type estimated
#'
#'
#'
#'
#' @examples
#'
#'
#' # Subsetting for time efficiency
#' tidybulk::se_mini |> deconvolve_cellularity(cores = 1)
#'
#'
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#' @export
#'
setGeneric("deconvolve_cellularity", function(.data,
                                              .sample = NULL,
                                              .transcript = NULL,
                                              .abundance = NULL,
                                              reference = NULL,
                                              method = "cibersort",
                                              prefix = "",
                                              action = "add",
                                              ...)
  standardGeneric("deconvolve_cellularity"))


#' Aggregates multiple counts from the same samples (e.g., from isoforms), concatenates other character columns, and averages other numeric columns
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description aggregate_duplicates() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with aggregated transcripts that were duplicated.
#'
#' @importFrom rlang enquo
#'
#'
#' @name aggregate_duplicates
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param aggregation_function A function for counts aggregation (e.g., sum,  median, or mean)
#' @param keep_integer A boolean. Whether to force the aggregated counts to integer
#'
#' @details This function aggregates duplicated transcripts (e.g., isoforms, ensembl).
#' For example, we often have to convert ensembl symbols to gene/transcript symbol,
#'  but in doing so we have to deal with duplicates. `aggregate_duplicates` takes a tibble
#'  and column names (as symbols; for `sample`, `transcript` and `count`) as arguments and
#'  returns a tibble with aggregate transcript with the same name. All the rest of the column
#'  are appended, and factors and boolean are appended as characters.
#'
#'  Underlying custom method:
#'  data |>
#' 		filter(n_aggr > 1) |>
#' 		group_by(!!.sample,!!.transcript) |>
#' 		dplyr::mutate(!!.abundance := !!.abundance |> aggregation_function())
#'
#' @return A consistent object (to the input) with aggregated transcript abundance and annotation
#'
#'
#'
#'
#' @examples
#'
#' # Create a aggregation column
#' se_mini = tidybulk::se_mini
#' SummarizedExperiment::rowData(se_mini )$gene_name = rownames(se_mini )
#'
#'    aggregate_duplicates(
#'      se_mini,
#'    .transcript = gene_name
#'    )
#'
#'
#' @docType methods
#' @rdname aggregate_duplicates-methods
#' @export
#'
#'
setGeneric("aggregate_duplicates", function(.data,
                                            
                                            .sample = NULL,
                                            .transcript = NULL,
                                            .abundance = NULL,
                                            aggregation_function = sum,
                                            keep_integer = TRUE)
  standardGeneric("aggregate_duplicates"))


#' Adjust transcript abundance for unwanted variation
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description adjust_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with an additional adjusted abundance column. This method uses scaled counts if present.
#'
#' @importFrom rlang enquo
#'
#'
#' @name adjust_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .factor_unwanted A tidy select, e.g. column names without double quotation. c(batch, country) These are the factor that we want to adjust for, including unwanted batcheffect, and unwanted biological effects.
#' @param .factor_of_interest A tidy select, e.g. column names without double quotation. c(treatment) These are the factor that we want to preserve.
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A character string. Methods include combat_seq (default), combat and limma_remove_batch_effect.
#'
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function sva::ComBat
#'
#' @param .formula DEPRECATED - A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_interest + batch)
#' @param transform DEPRECATED - A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param inverse_transform DEPRECATED - A function that is the inverse of transform (e.g. expm1 is inverse of log1p). This is needed to tranform back the counts after analysis.
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function adjusts the abundance for (known) unwanted variation.
#' At the moment just an unwanted covariate is allowed at a time using Combat (DOI: 10.1093/bioinformatics/bts034)
#'
#' Underlying method:
#' 	sva::ComBat(data, batch = my_batch,	mod = design,	prior.plots = FALSE, ...)
#'
#' @return A consistent object (to the input) with additional columns for the adjusted counts as `<COUNT COLUMN>_adjusted`
#'
#'
#'
#'
#' @examples
#'
#'
#'
#' cm = tidybulk::se_mini
#' cm$batch = 0
#' cm$batch[colnames(cm) %in% c("SRR1740035", "SRR1740043")] = 1
#'
#' cm |>
#' identify_abundant() |>
#'	adjust_abundance(	.factor_unwanted = batch, .factor_of_interest =  condition, method="combat"	)
#'
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#' @export
#'
#'
setGeneric("adjust_abundance", function(.data,
                                        
                                        # DEPRECATED
                                        .formula = NULL,
                                        .factor_unwanted =NULL,
                                        .factor_of_interest = NULL,
                                        .sample = NULL,
                                        .transcript = NULL,
                                        .abundance = NULL,
                                        method = "combat_seq",
                                        
                                        action = "add",
                                        ...,
                                        
                                        # DEPRECATED
                                        log_transform = NULL,
                                        transform = NULL,
                                        inverse_transform = NULL
                                        
)
standardGeneric("adjust_abundance"))


#' Drop redundant elements (e.g., samples) for which feature (e.g., transcript/gene) abundances are correlated
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description remove_redundancy() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) for correlation method or | <DIMENSION 1> | <DIMENSION 2> | <...> | for reduced_dimensions method, and returns a consistent object (to the input) with dropped elements (e.g., samples).
#'
#' @importFrom rlang enquo
#'
#'
#' @name remove_redundancy
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The method to use, correlation and reduced_dimensions are available. The latter eliminates one of the most proximar pairs of samples in PCA reduced dimensions.
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param top An integer. How many top genes to select for correlation based method
#' @param Dim_a_column A character string. For reduced_dimension based calculation. The column of one principal component
#' @param Dim_b_column A character string. For reduced_dimension based calculation. The column of another principal component
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function removes redundant elements from the original data set (e.g., samples or transcripts).
#' For example, if we want to define cell-type specific signatures with low sample redundancy.
#' This function returns a tibble with dropped redundant elements (e.g., samples).
#' Two redundancy estimation approaches are supported:
#' (i) removal of highly correlated clusters of elements (keeping a representative) with method="correlation";
#' (ii) removal of most proximal element pairs in a reduced dimensional space.
#'
#' Underlying method for correlation:
#' widyr::pairwise_cor(sample, transcript,count, sort = TRUE, diag = FALSE, upper = FALSE)
#'
#' Underlying custom method for reduced dimensions:
#' select_closest_pairs = function(df) {
#' 		couples <- df |> head(n = 0)
#'
#' 		while (df |> nrow() > 0) {
#' 			pair <- df |>
#' 			arrange(dist) |>
#' 			head(n = 1)
#' 			couples <- couples |> bind_rows(pair)
#' 			df <- df |>
#' 				filter(
#' 					!`sample 1` %in% (pair |> select(1:2) |> as.character()) &
#' 						!`sample 2` %in% (pair |> select(1:2) |> as.character())
#' 				)
#' 		}
#'
#' 		couples
#'
#' 	}
#'
#'
#'
#' @return A tbl object with with dropped redundant elements (e.g., samples).
#'
#' @examples
#'
#'
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'    remove_redundancy(
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.abundance =  count,
#' 	   	method = "correlation"
#' 	   	)
#'
#' counts.MDS =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'   reduce_dimensions( method="MDS", .dims = 3)
#'
#' remove_redundancy(
#' 	counts.MDS,
#' 	Dim_a_column = `Dim1`,
#' 	Dim_b_column = `Dim2`,
#' 	.element = sample,
#'   method = "reduced_dimensions"
#' )
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#' @export
#'
#'
setGeneric("remove_redundancy", function(.data,
                                         .element = NULL,
                                         .feature = NULL,
                                         .abundance = NULL,
                                         method,
                                         
                                         of_samples = TRUE,
                                         
                                         correlation_threshold = 0.9,
                                         top = Inf,
                                         transform = identity,
                                         Dim_a_column,
                                         Dim_b_column,
                                         
                                         # DEPRECATED
                                         log_transform = NULL
)
standardGeneric("remove_redundancy"))


#' Rotate two dimensions (e.g., principal components) of an arbitrary angle
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description rotate_dimensions() takes as input a `tbl` formatted as | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#'
#'
#' @name rotate_dimensions
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#'
#' @param dimension_1_column A character string. The column of the dimension 1
#' @param dimension_2_column  A character string. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param dimension_1_column_rotated A character string. The column of the rotated dimension 1 (optional)
#' @param dimension_2_column_rotated A character string. The column of the rotated dimension 2 (optional)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details This function to rotate two dimensions such as the reduced dimensions.
#'
#' Underlying custom method:
#' 	rotation = function(m, d) {
#' 		// r = the angle
#' 		// m data matrix
#'    r = d * pi / 180
#'    ((bind_rows(
#' 	  c(`1` = cos(r), `2` = -sin(r)),
#' 	  c(`1` = sin(r), `2` = cos(r))
#'   ) |> as_matrix()) %*% m)
#'  }
#'
#'
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
#'
#'
#' @examples
#'
#' counts.MDS =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'  reduce_dimensions( method="MDS", .dims = 3)
#'
#' counts.MDS.rotated =  rotate_dimensions(counts.MDS, `Dim1`, `Dim2`, rotation_degrees = 45, .element = sample)
#'
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#' @export
#'
setGeneric("rotate_dimensions", function(.data,
                                         dimension_1_column,
                                         dimension_2_column,
                                         rotation_degrees,
                                         .element = NULL,
                                         of_samples = TRUE,
                                         dimension_1_column_rotated = NULL,
                                         dimension_2_column_rotated = NULL,
                                         action = "add")
  standardGeneric("rotate_dimensions"))


#' Dimension reduction of the transcript abundance data
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description reduce_dimensions() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and calculates the reduced dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#'
#'
#' @name reduce_dimensions
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The dimension reduction algorithm to use (PCA, MDS, tSNE).
#' @param top An integer. How many top genes to select for dimensionality reduction
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param .dims An integer. The number of dimensions your are interested in (e.g., 4 for returning the first four principal components).
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean for method="PCA", this will be passed to the `prcomp` function. It is not included in the ... argument because although the default for `prcomp` if FALSE, it is advisable to set it as TRUE.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function prcomp if you choose method="PCA" or Rtsne if you choose method="tSNE", or uwot::tumap if you choose method="umap"
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function reduces the dimensions of the transcript abundances.
#' It can use multi-dimensional scaling (MDS; DOI.org/10.1186/gb-2010-11-3-r25),
#' principal component analysis (PCA), or tSNE (Jesse Krijthe et al. 2018)
#'
#' Underlying method for PCA:
#' prcomp(scale = scale, ...)
#'
#' Underlying method for MDS:
#' limma::plotMDS(ndim = .dims, plot = FALSE, top = top)
#'
#' Underlying method for tSNE:
#' Rtsne::Rtsne(data, ...)
#'
#' Underlying method for UMAP:
#'
#'  df_source =
#' .data |>
#'
#'   # Filter NA symbol
#'   filter(!!.feature |> is.na() |> not()) |>
#'
#'   # Prepare data frame
#'   distinct(!!.feature,!!.element,!!.abundance) |>
#'
#'   # Filter most variable genes
#'   keep_variable_transcripts(top) |>
#'   reduce_dimensions(method="PCA",  .dims = calculate_for_pca_dimensions,  action="get" ) |>
#'   as_matrix(rownames = quo_name(.element)) |>
#'   uwot::tumap(...)
#'
#'
#' @return A tbl object with additional columns for the reduced dimensions
#'
#'
#' @examples
#'
#'
#'
#' counts.MDS =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'  reduce_dimensions( method="MDS", .dims = 3)
#'
#'
#' counts.PCA =
#'  tidybulk::se_mini |>
#'  identify_abundant() |>
#'  reduce_dimensions(method="PCA", .dims = 3)
#'
#'
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#' @export
#'
#'
setGeneric("reduce_dimensions", function(.data,
                                         .element = NULL,
                                         .feature = NULL,
                                         .abundance = NULL,
                                         method,
                                         .dims = 2,
                                         
                                         top = 500,
                                         of_samples = TRUE,
                                         transform = log1p,
                                         scale = TRUE,
                                         action = "add",
                                         ...,
                                         
                                         # DEPRECATED
                                         log_transform = NULL
                                         
)
standardGeneric("reduce_dimensions"))



#' Get clusters of elements (e.g., samples or transcripts)
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description cluster_elements() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and identify clusters in the data.
#'
#' @importFrom rlang enquo
#'
#'
#' @name cluster_elements
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .abundance The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The cluster algorithm to use, at the moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function kmeans
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details identifies clusters in the data, normally of samples.
#' This function returns a tibble with additional columns for the cluster annotation.
#' At the moment only k-means (DOI: 10.2307/2346830) and SNN clustering (DOI:10.1016/j.cell.2019.05.031) is supported, the plan is to introduce more clustering methods.
#'
#' Underlying method for kmeans
#' do.call(kmeans(.data, iter.max = 1000, ...)
#'
#' Underlying method for SNN
#' .data %>%
#' Seurat::CreateSeuratObject() %>%
#' Seurat::ScaleData(display.progress = TRUE,num.cores = 4, do.par = TRUE) %>%
#' Seurat::FindVariableFeatures(selection.method = "vst") %>%
#' Seurat::RunPCA(npcs = 30) %>%
#' Seurat::FindNeighbors() %>%
#' Seurat::FindClusters(method = "igraph", ...)
#'
#'
#' @return A tbl object with additional columns with cluster labels
#'
#'
#' @examples
#'
#'
#'     cluster_elements(tidybulk::se_mini,	centers = 2, method="kmeans")
#'
#' @docType methods
#' @rdname cluster_elements-methods
#' @export
#'
setGeneric("cluster_elements", function(.data,
                                        .element = NULL,
                                        .feature = NULL,
                                        .abundance = NULL,
                                        method,
                                        of_samples = TRUE,
                                        transform = log1p,
                                        
                                        action = "add",
                                        ...,
                                        
                                        # DEPRECATED
                                        log_transform = NULL
)
standardGeneric("cluster_elements"))


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
                                                    .sample = NULL,
                                                    .transcript = NULL,
                                                    .abundance = NULL,
                                                    method = "limma_normalize_quantiles",
                                                    target_distribution = NULL,
                                                    action = "add")
  standardGeneric("quantile_normalise_abundance"))


#' Scale the counts of transcripts/genes
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description scale_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and Scales transcript abundance compansating for sequencing depth (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#'
#' @importFrom rlang enquo
#'
#' @importFrom stats median
#'
#' @name scale_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A character string. The scaling method passed to the back-end function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param reference_sample A character string. The name of the reference sample. If NULL the sample with highest total read count will be selected as reference.
#' @param .subset_for_scaling A gene-wise quosure condition. This will be used to filter rows (features/genes) of the dataset. For example
#' @param action A character string between "add" (default) and "only". "add" joins the new information to the input tbl (default), "only" return a non-redundant tbl with the just new information.
#'
#' @param reference_selection_function DEPRECATED. please use reference_sample.
#'
#' @details Scales transcript abundance compensating for sequencing depth
#' (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#' Lowly transcribed transcripts/genes (defined with minimum_counts and minimum_proportion parameters)
#' are filtered out from the scaling procedure.
#' The scaling inference is then applied back to all unfiltered data.
#'
#' Underlying method
#' edgeR::calcNormFactors(.data, method = c("TMM","TMMwsp","RLE","upperquartile"))
#'
#'
#'
#' @return A tbl object with additional columns with scaled data as `<NAME OF COUNT COLUMN>_scaled`
#'
#'
#' @examples
#'
#'
#'  tidybulk::se_mini |>
#'    identify_abundant() |>
#'    scale_abundance()
#'
#'
#'
#' @docType methods
#' @rdname scale_abundance-methods
#' @export

setGeneric("scale_abundance", function(.data,
                                       .sample = NULL,
                                       .transcript = NULL,
                                       .abundance = NULL,
                                       method = "TMM",
                                       reference_sample = NULL,
                                       .subset_for_scaling = NULL,
                                       action = "add",
                                       
                                       # DEPRECATED
                                       reference_selection_function = NULL)
  standardGeneric("scale_abundance"))


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
                                               .sample = NULL,
                                               .transcript = NULL,
                                               .abundance = NULL)
  standardGeneric("as_SummarizedExperiment"))


.as_SummarizedExperiment = function(.data,
                                    .sample = NULL,
                                    .transcript = NULL,
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
  .abundance_scaled =
    .data %>%
    ifelse_pipe(
      ".abundance_scaled" %in% ((.) %>% get_tt_columns() %>% names) &&
        # .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% not() &&
        quo_name((.) %>% get_tt_columns() %$% .abundance_scaled) %in% ((.) %>% colnames),
      ~ .x %>% get_tt_columns() %$% .abundance_scaled,
      ~ NULL
    )
  
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
setMethod("as_SummarizedExperiment", "spec_tbl_df", .as_SummarizedExperiment)

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