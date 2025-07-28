#' Aggregates multiple counts from the same samples (e.g., from isoforms), concatenates other character columns, and averages other numeric columns
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description aggregate_duplicates() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with aggregated transcripts that were duplicated.
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_null quo_name
#' @importFrom dplyr group_by summarise arrange select left_join
#' @importFrom dplyr across select_if as_tibble pull
#' @importFrom methods is
#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment rowData colData SummarizedExperiment
#' @importFrom SummarizedExperiment rowRanges assays
#' @importFrom SummarizedExperiment "rowRanges<-"
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom S4Vectors as.data.frame
#' @importFrom utils capture.output
#' @importFrom crayon blue
#'
#'
#' @name aggregate_duplicates
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param aggregation_function A function for counts aggregation (e.g., sum,  median, or mean)
#' @param keep_integer A boolean. Whether to force the aggregated counts to integer
#' @param ... Additional arguments passed to the aggregation function
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
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Lawrence, M., Huber, W., Pagès, H., Aboyoun, P., Carlson, M., Gentleman, R., Morgan, M. T., & Carey, V. J. (2013). Software for computing and annotating genomic ranges. PLoS Computational Biology, 9(8), e1003118. doi:10.1371/journal.pcbi.1003118
#'
#' Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 1.1.0. https://CRAN.R-project.org/package=dplyr
#'
#' @docType methods
#' @rdname aggregate_duplicates-methods
#' @export
#'
#'
setGeneric("aggregate_duplicates", function(.data,
                                            
                                            
                                            .transcript = NULL,
                                            .abundance = NULL,
                                            aggregation_function = sum,
                                            keep_integer = TRUE,
                                            ...)
  standardGeneric("aggregate_duplicates"))



#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom dplyr setdiff
#' @importFrom dplyr select_if
.aggregate_duplicates_se = function(.data,
                                    
                                    
                                    
                                    .abundance = NULL,
                                    aggregation_function = sum,
                                    keep_integer = TRUE,
                                    ...) {
  
  # Fix NOTEs
  . = NULL
  
  # Make col names
  .transcript = enquo(.transcript)
  
  
  if(quo_is_null(.transcript)) stop("tidybulk says: using SummarizedExperiment with aggregate_duplicates, you need to specify .transcript parameter. It should be a feature-wise column (e.g. gene symbol) that you want to collapse he features with (e.g. ensembl). It cannot be the representation of rownames(SummarizedExperiment), as those are unique by definition, and not part of rowData per-se.")
  
  if(!quo_name(.transcript) %in% colnames( .data |> rowData()))
    stop("tidybulk says: the .transcript argument must be a feature-wise column names. The feature-wise information can be found with rowData()")
  if( !is.null(.abundance))
    warning("tidybulk says: for SummarizedExperiment objects only the argument .transcript (feature ID to collapse) is considered")
  
  collapse_function = function(x){ x |> unique() |> paste(collapse = "___")	}
  
  
  # Non standard column classes
  non_standard_columns =
    .data |>
    rowData() |>
    as_tibble() |>
    select_if(select_non_standard_column_class) |>
    colnames()
  
  # GRanges
  columns_to_collapse =
    .data |>
    rowData() |>
    colnames() |>
    outersect(non_standard_columns) |>
    setdiff(quo_name(.transcript)) |>
    c(feature__$name)
  # when(
  #   !is.null(rownames(.data)) ~ c(., feature__$name),
  #   ~ (.)
  # )
  
  # Row data
  new_row_data =
    .data |>
    rowData() |>
    S4Vectors::as.data.frame(optional = TRUE) |>
    tibble::as_tibble(rownames = feature__$name, .name_repair = "minimal") |> 
    group_by(!!.transcript) |>
    summarise(
      across(columns_to_collapse, ~ .x |> collapse_function()),
      across(non_standard_columns, ~ .x[1]),
      merged_transcripts = n()
    ) |>
    
    arrange(!!as.symbol(feature__$name)) |>
    as.data.frame() |>
    (\(.) {
      .x = (.)
      rownames(.x) = .x[,feature__$name]
      .x = .x |> select(-feature__$name)
      .x
    })()
  
  # If no duplicate exit
  if(!nrow(new_row_data)<nrow(rowData(.data))){
    message(sprintf("tidybulk says: your object does not have duplicates along the %s column. The input dataset is returned.", quo_name(.transcript)))
    return(.data)
  }
  
  # Counts
  new_count_data =
    .data |>
    assays() |>
    as.list() |>
    map(
      ~ {
        is_data_frame = .x |> is("data.frame")
        if(is_data_frame) .x = .x |> as.matrix()
        
        # Gove duplicated rownames
        rownames(.x) = rowData(.data)[,quo_name(.transcript)]
        
        # Combine
        if(rownames(.x) |> is.na() |> which() |> length() |> gt(0))
          stop(sprintf("tidybulk says: you have some %s that are NAs", quo_name(.transcript)))
        
        .x =  combineByRow(.x, aggregation_function)
        .x = .x[match(new_row_data[,quo_name(.transcript)], rownames(.x)),,drop=FALSE]
        rownames(.x) = rownames(new_row_data)
        
        if(is_data_frame) .x = .x |> as.data.frame()
        .x
      }
    )
  
  if(!is.null(rowRanges(.data))){
    
    new_range_data = rowRanges(.data) |> as_tibble()
    
    # If GRangesList & and .transcript is not there add .transcript
    if(is(rowRanges(.data), "CompressedGRangesList") & !quo_name(.transcript) %in% colnames(new_range_data)){
      
      new_range_data =
        new_range_data |> left_join(
          rowData(.data)[,quo_name(.transcript),drop=FALSE] |>
            as_tibble(rownames = feature__$name) ,
          by=c("group_name" = feature__$name)
        ) |>
        select(-group_name, -group)
    }
    
    # Through warning if there are logicals of factor in the data frame
    # because they cannot be merged if they are not unique
    if (length(non_standard_columns)>0 & new_range_data |>  pull(!!.transcript) |> duplicated() |> which() |> length() |> gt(0) ) {
      warning(paste(capture.output({
        cat(crayon::blue("tidybulk says: If duplicates exist from the following columns, only the first instance was taken (lossy behaviour), as aggregating those classes with concatenation is not possible.\n"))
        print(rowData(.data)[1,non_standard_columns,drop=FALSE])
      }), collapse = "\n"))
    }
    
    new_range_data = new_range_data |>
      
      # I have to use this trick because rowRanges() and rowData() share @elementMetadata
      select(-any_of(colnames(new_row_data) |> outersect(quo_name(.transcript)))) |>
      suppressWarnings()
    
    
    #if(is(rr, "CompressedGRangesList") | nrow(new_row_data)<nrow(rowData(.data))) {
    new_range_data = makeGRangesListFromDataFrame(
      new_range_data,
      split.field = quo_name(.transcript),
      keep.extra.columns = TRUE
    )
    
    # Give back rownames
    new_range_data = new_range_data[match(new_row_data[,quo_name(.transcript)], names(new_range_data))]
    #names(new_range_data) = rownames(new_row_data)
    #}
    # else if(is(rr, "GRanges")) new_range_data = makeGRangesFromDataFrame(new_range_data, keep.extra.columns = TRUE)
    # else stop("tidybulk says: riowRanges should be either GRanges or CompressedGRangesList. Or am I missing something?")
    
  }
  
  # Build the object
  .data_collapsed =
    SummarizedExperiment(
      assays = new_count_data,
      colData = colData(.data)
    )
  
  if(!is.null(rowRanges(.data))) rowRanges(.data_collapsed) = new_range_data
  
  rowData(.data_collapsed) = new_row_data
  
  .data_collapsed
  
}

#' aggregate_duplicates
#'
#' @docType methods
#' @rdname aggregate_duplicates-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("aggregate_duplicates",
          "SummarizedExperiment",
          .aggregate_duplicates_se)

#' aggregate_duplicates
#'
#' @docType methods
#' @rdname aggregate_duplicates-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("aggregate_duplicates",
          "RangedSummarizedExperiment",
          .aggregate_duplicates_se)


