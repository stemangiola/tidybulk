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
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param method A character string. The scaling method passed to the back-end function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param reference_sample A character string. The name of the reference sample. If NULL the sample with highest total read count will be selected as reference.
#' @param .subset_for_scaling A gene-wise quosure condition. This will be used to filter rows (features/genes) of the dataset. For example
#' @param suffix A character string to append to the scaled abundance column name. Default is "_scaled".
#'
#' @param reference_selection_function DEPRECATED. please use reference_sample.
#' @param ... Further arguments.
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
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
                                       
                                       
                                       abundance = assayNames(.data)[1],
                                       method = "TMM",
                                       reference_sample = NULL,
                                       .subset_for_scaling = NULL,
                                       suffix = "_scaled",
                                       # DEPRECATED
                                       reference_selection_function = NULL,
                                       ...,
                                       .abundance = NULL)
  standardGeneric("scale_abundance"))

#' @importFrom magrittr multiply_by
#' @importFrom magrittr divide_by
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom utils tail
#' @importFrom stats na.omit
#' @importFrom stringr str_c
#' @importFrom dplyr mutate select pull arrange slice n_distinct across
#' @importFrom tidyr nest unnest pivot_longer pivot_wider drop_na
#' @importFrom purrr map when
#' @importFrom stringr str_subset str_remove str_replace str_replace_all
#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @importFrom purrr map2_dfr
#' @importFrom SummarizedExperiment "colData<-"
#' @importFrom SummarizedExperiment "rowData<-"
#' @importFrom magrittr "%$%"
#' @importFrom SummarizedExperiment "assays<-"
#'
.scale_abundance_se = function(.data,
                               
                               
                               abundance = assayNames(.data)[1],
                               method = "TMM",
                               reference_sample = NULL,
                               .subset_for_scaling = NULL,
                               suffix = "_scaled",
                               # DEPRECATED
                               reference_selection_function = NULL,
                               ...,
                               .abundance = NULL) {
  
  # Fix NOTEs
  . = NULL
  
  # Soft-deprecate .abundance, prefer abundance (character)
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "scale_abundance(.abundance)", "scale_abundance(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }
  my_assay <- abundance
  
  # Check if package is installed, otherwise install
  check_and_install_packages("edgeR")
  
  
  # DEPRECATION OF reference function
  if (is_present(reference_selection_function) & !is.null(reference_selection_function)) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.1.8", "tidybulk::scale_abundance(reference_selection_function = )", details = "The argument reference_selection_function is now deprecated please use reference_sample. By default the reference selection function is max()")
    
  }
  
  # Check that reference sample exists
  if(!is.null(reference_sample) && !reference_sample %in% (.data %>% colnames))
    stop("tidybulk says: your reference sample is not among the samples in your data frame")
  
  .subset_for_scaling = enquo(.subset_for_scaling)
  
  
  .data_filtered =
    filter_if_abundant_were_identified(.data)
  
  
  if (!quo_is_null(.subset_for_scaling))
    .data_filtered = filter_genes_on_condition(.data_filtered, !!.subset_for_scaling)
  
  # Filter based on user condition
  
  # Check I have genes left
  if (nrow(.data_filtered) == 0)
    stop("tidybulk says: there are 0 genes that passes the filters (.abundant and/or .subset_for_scaling). Please check your filtering or your data.")
  
  # Determine the correct assay name
  my_counts_filtered = assays(.data_filtered)[[my_assay]] %>% na.omit()
  library_size_filtered = my_counts_filtered %>% colSums(na.rm  = TRUE)
  
  # If not enough genes, warning
  if(nrow(my_counts_filtered)<100) warning(warning_for_scaling_with_few_genes)
  
  # Set column name for value scaled
  value_scaled = paste0(my_assay, suffix)
  
  # Get reference
  reference <-
    reference_sample
  
  if (is.null(reference))
    reference = library_size_filtered %>%
    sort() %>%
    tail(1) %>%
    names()
  
  
  # Communicate the reference if chosen by default
  if(is.null(reference_sample)) message(sprintf("tidybulk says: the sample with largest library size %s was chosen as reference for scaling", reference))
  
  # Calculate TMM
  nf <-
    edgeR::calcNormFactors(
      my_counts_filtered,
      refColumn = reference,
      method = method
    )
  
  # Calculate multiplier
  multiplier =
    # Relecting the ratio of effective library size of the reference sample to the effective library size of each sample
    (library_size_filtered[reference] * nf[reference]) %>%
    divide_by(library_size_filtered * nf) 
  
  # NOT HELPING - Put everything to the reference sample scale
  # multiply_by(library_size_filtered[reference])
  
  # At the moment no because would be different from TIBBLE behaviour
  # %>%
  #
  # # Make reference == 1
  # divide_by(.[reference])
  
  # Add to sample info
  colData(.data)$TMM = nf
  colData(.data)$multiplier = multiplier
  
  my_counts_scaled =
    list(
      assay(.data) %*%
        diag(multiplier)
      
    ) %>%
    setNames(value_scaled)
  colnames(my_counts_scaled[[1]]) = assay(.data)  |> colnames()
  
  
  # Add the assay
  assays(.data, withDimnames=FALSE) =  assays(.data) %>% c(my_counts_scaled)
  
  .data %>%
    
    # Add methods
    memorise_methods_used(c("edger", "tmm")) %>%
    
    # Attach column internals
    add_tt_columns(.abundance_scaled = !!(((function(x, v)	enquo(v))(x,!!as.symbol(value_scaled))) |> drop_enquo_env()) )
  
}

#' scale_abundance
#'
#' @docType methods
#' @rdname scale_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("scale_abundance",
          "SummarizedExperiment",
          .scale_abundance_se)

#' scale_abundance
#'
#' @docType methods
#' @rdname scale_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("scale_abundance",
          "RangedSummarizedExperiment",
          .scale_abundance_se)
