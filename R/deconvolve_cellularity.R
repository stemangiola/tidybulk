#' Get cell type proportions from samples
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description deconvolve_cellularity() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with the estimated cell type abundance for each sample
#'
#' @importFrom rlang enquo quo_name
#' @importFrom dplyr select pull mutate as_tibble
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment colData assays
#' @importFrom magrittr %$%
#'
#'
#' @name deconvolve_cellularity
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .abundance The name of the transcript/gene abundance column
#' @param reference A data frame. The methods cibersort and llsr can accept a custom rectangular dataframe with genes as rows names, cell types as column names and gene-transcript abundance as values. For exampler tidybulk::X_cibersort. The transcript/cell_type data frame of integer transcript abundance. If NULL, the default reference for each algorithm will be used. For llsr will be LM22.
#' @param method A character string. The method to be used. Available methods: "cibersort", "llsr", "epic", "mcp_counter", "quantiseq", "xcell". If a vector is provided, an error will be thrown. Default is all available methods.
#' @param prefix A character string. The prefix you would like to add to the result columns. It is useful if you want to reshape data.
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
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453-457. doi:10.1038/nmeth.3337
#'
#' Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife, 6, e26476. doi:10.7554/eLife.26476
#'
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#' @export
#'
setGeneric("deconvolve_cellularity", function(.data,
                                              
                                              
                                              .abundance = NULL,
                                              reference = NULL,
                                              method = "cibersort",
                                              prefix = "",
                                              
                                              ...)
  standardGeneric("deconvolve_cellularity"))

#' Deconvolve bulk data to sample-level cell–type proportions
#'
#' Wrapper that applies a chosen deconvolution engine (CIBERSORT, LLSR, EPIC,
#' MCPcounter, quanTIseq, or xCell via immunedeconv) to the active assay in a
#' `SummarizedExperiment`, then appends the estimated proportions to `colData`.
#' You may supply a custom reference matrix; otherwise sensible defaults are
#' used per method.
#'
#' @param .data    A `SummarizedExperiment`.
#' @param reference Reference matrix or method-specific handle (see Details).
#' @param method   Character string naming the deconvolution method.
#' @param prefix   Optional prefix to prepend to output column names.
#' @param ...      Additional arguments passed through to the underlying
#'                 deconvolution function.
#'
#' @return The input `SummarizedExperiment` with additional proportion columns
#'         in `colData`.
#'
#' @details
#' * `"cibersort"` uses `my_CIBERSORT()` and requires the packages *class*,
#'   *e1071*, and *preprocessCore*.
#' * `"llsr"` uses `run_llsr()`.
#' * `"epic"` uses `run_epic()`.
#' * `"mcp_counter"`, `"quantiseq"`, and `"xcell"` use `immunedeconv::deconvolute()`
#'   and require that *immunedeconv* is attached.
#'
#' @family deconvolution helpers
#'

#' @importFrom purrr map
#' @importFrom rlang quo_is_symbolic quo_name
#' @importFrom tibble as_tibble
#' @importFrom dplyr select
#' @importFrom tidyr gather spread
#' @importFrom SummarizedExperiment assays colData
#' @keywords internal
.deconvolve_cellularity_se = function(.data,
                                      reference = X_cibersort,
                                      method = c("cibersort", "llsr", "epic", "mcp_counter", "quantiseq", "xcell"),
                                      prefix = "",
                                      ...) {
  
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
  
  .sample = s_(.data)$symbol
  
  my_assay =
    .data |>
    assays() |>
    as.list()
  my_assay = my_assay[[get_assay_scaled_if_exists_SE(.data)]]
  
  # 	  # Change row names
  # 	  if (quo_is_symbolic(.transcript)) {
  #   	    .x = (.)
  #   	    rownames(.x) = .data %>% pivot_transcript() %>% pull(!!.transcript)
  #   	    .x
  #   	  } else {
  #   	    (.)
  #   	  }
  
  # Get the dots arguments
  dots_args = rlang::dots_list(...)
  
  my_proportions = 
    if (method |> tolower() |> equals("cibersort")) {
      
      # Check if package is installed, otherwise install
      check_and_install_packages(c("class", "e1071", "preprocessCore"))
      
      # Choose reference
      reference = if (is.null(reference)) X_cibersort else reference
      
      # Validate reference
      validate_signature_SE(my_assay, reference)
      
      do.call(my_CIBERSORT, list(Y = my_assay, X = reference, QN=FALSE) |> c(dots_args)) %$%
        proportions |>
        as_tibble(rownames = quo_name(.sample)) |>
        select(-`P-value`,-Correlation,-RMSE)
    } else if (method |> tolower() |> equals("llsr")) {
      
      # Choose reference
      reference = if (is.null(reference)) X_cibersort else reference
      
      # Validate reference
      validate_signature_SE(my_assay, reference)
      
      my_assay |>
        run_llsr(reference, ...) |>
        as_tibble(rownames = quo_name(.sample))
    } else if (method |> tolower() |> equals("epic")) {
      
      # Choose reference
      reference = if (is.null(reference)) "BRef" else reference
      
      my_assay |>
        run_epic(reference) |>
        as_tibble(rownames = quo_name(.sample))
    } else if (method |> tolower() %in% c("mcp_counter", "quantiseq", "xcell")) {
      
      # Check if package is installed, otherwise install
      check_and_install_packages("immunedeconv")
      
      if(method %in% c("mcp_counter", "quantiseq", "xcell") & !"immunedeconv" %in% (.packages()))
        stop("tidybulk says: for xcell, mcp_counter, or quantiseq deconvolution you should have the package immunedeconv attached. Please execute library(immunedeconv)")
      
      my_assay |>
        deconvolute(method |> tolower(), tumor = FALSE) |>
        gather(!!.sample, .proportion, -cell_type) |>
        spread(cell_type,  .proportion)
    } else {
      stop(
        "tidybulk says: please choose between cibersort, llsr and epic methods"
      )
    }
    
    # Parse results and return
    if (is.null(my_proportions)) {
      stop("tidybulk says: The deconvolution method returned NULL. Please check your input data and method parameters.")
    }
    col_names = colnames(my_proportions)
    sample_col = col_names[1]
    other_cols = col_names[-1]
    new_names = sprintf("%s%s", prefix, other_cols)
    my_proportions = setNames(my_proportions, c(sample_col, new_names))
  
  # Att proportions
    colData(.data) = colData(.data) |> cbind(
    my_proportions |>
    as_matrix(rownames = .sample) |>
      (function(matrix_data) {
        matrix_data[match(rownames(colData(.data)), rownames(matrix_data)),]
      })()
  )
  
  .data |>
    
    # Attach attributes
    memorise_methods_used(tolower(method))
  
}

#' deconvolve_cellularity
#'
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("deconvolve_cellularity",
          "SummarizedExperiment",
          .deconvolve_cellularity_se)

#' deconvolve_cellularity
#'
#' @importFrom rlang inform
#'
#' @docType methods
#' @rdname deconvolve_cellularity-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "deconvolve_cellularity",
  "RangedSummarizedExperiment",
  .deconvolve_cellularity_se
)


