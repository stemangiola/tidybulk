eliminate_GRanges_metadata_columns_also_present_in_Rowdata = function(.my_data, se){
  .my_data |>
    select(-any_of(colnames(rowData(se)))) |>

    # In case there is not metadata column
    suppressWarnings()
}


#' @importFrom dplyr select
#' @importFrom tidyselect any_of
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom tibble rowid_to_column
#'
#' @noRd
get_special_datasets <- function(se) {

  rr =  se |>
    rowRanges()

  rr |>
    when(

      # If no ranges
      as.data.frame(.) |>
        nrow() |>
        equals(0) ~ tibble(),

      # If it is a range list (multiple rows per feature)
      is(., "CompressedGRangesList") ~ {

        # If GRanges does not have row names
        if(is.null(rr@partitioning@NAMES)) rr@partitioning@NAMES = as.character(1:nrow(se))

        tibble::as_tibble(rr) |>
          eliminate_GRanges_metadata_columns_also_present_in_Rowdata(se) |>
          nest(GRangesList = -group_name) |>
          rename(!!f_(se)$symbol := group_name)

      },

      # If standard GRanges (one feature per line)
      ~ {

        # If GRanges does not have row names
        if(is.null(rr@ranges@NAMES)) rr@ranges@NAMES = as.character(1:nrow(se))

        tibble::as_tibble(rr) |>
          eliminate_GRanges_metadata_columns_also_present_in_Rowdata(se) |>
          mutate(!!f_(se)$symbol := rr@ranges@NAMES)
      }

    ) %>%
    list()

}

#' @importFrom stringr str_replace
change_reserved_column_names = function(col_data, .data ){

  col_data %>%

    setNames(
      colnames(.) |>
        sapply(function(x) if(x==f_(.data)$name) sprintf("%s.x", f_(.data)$name) else x) %>%
        sapply(function(x) if(x==s_(.data)$name) sprintf("%s.x", s_(.data)$name) else x) %>%
        str_replace("^coordinate$", "coordinate.x")
    )

}

#' @importFrom tidyr gather
#' @importFrom dplyr rename
#' @importFrom dplyr left_join
#' @importFrom tibble as_tibble
#' @importFrom purrr reduce
#' @importFrom SummarizedExperiment assays
#'
#' @noRd
get_count_datasets <- function(se) {
  map2(
    assays(se) |> as.list(),
    names(assays(se)),
    ~ {

      # If the counts are in a sparse matrix convert to a matrix
      # This might happen because the user loaded tidySummarizedExperiment and is
      # print a SingleCellExperiment
      if(is(.x, "dgCMatrix")) {
        .x = as.matrix(.x)
      }

      .x |>
        # matrix() %>%
        # as.data.frame() %>%
        tibble::as_tibble(rownames = f_(se)$name, .name_repair = "minimal") |>

        # If the matrix does not have sample names, fix column names
        when(colnames(.x) |> is.null() ~ setNames(., c(
          f_(se)$name,  seq_len(ncol(.x))
        )),
        ~ (.)
        ) %>%

        gather(!!s_(se)$symbol, !!.y,-!!f_(se)$symbol)

      #%>%
      #  rename(!!.y := count)
    }) %>%
    when( 
      length(.)>0 ~ bind_cols(.,  .name_repair = c("minimal")) %>% .[!duplicated(colnames(.))], # reduce(., left_join, by = c(f_(se)$name, s_(se)$name)),
      ~ expand.grid(
        rownames(se), colnames(se)
      ) %>%
        setNames(c(f_(se)$name, s_(se)$name)) %>%
        tibble::as_tibble()
    ) %>%

    # Add dummy sample or feature if we have empty assay.
    # This is needed for a correct isualisation of the tibble form
    when(
      f_(se)$name %in% colnames(.) |> not() ~ mutate(., !!f_(se)$symbol := as.character(NA)),
      s_(se)$name %in% colnames(.) |> not() ~ mutate(., !!s_(se)$symbol := as.character(NA)),
      ~ (.)
    )
}

subset_tibble_output = function(.data, count_info, sample_info, gene_info, range_info, .subset){
  # This function outputs a tibble after subsetting the columns
  .subset = enquo(.subset)

  # Build template of the output
  output_colnames =
    slice(count_info, 0) %>%
    left_join(slice(sample_info, 0), by=s_(.data)$name) %>%
    left_join(slice(gene_info, 0), by = f_(.data)$name) %>%
    when(nrow(range_info) > 0 ~ (.) %>% left_join(range_info, by=f_(.data)$name), ~ (.)) %>%
    select(!!.subset) %>%
    colnames()


  # Sample table
  sample_info =
    sample_info %>%
    when(
      colnames(.) |> intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., any_of(s_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )

  # Ranges table
  range_info =
    range_info %>%
    when(
      colnames(.) |> intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., any_of(f_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )

  # Ranges table
  gene_info =
    gene_info %>%
    when(
      colnames(.) |> intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., any_of(f_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )

  # Ranges table
  count_info =
    count_info %>%
    when(
      colnames(.) |> intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., any_of(f_(.data)$name, s_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )


  if(
    !is.null(count_info) &
    (
      !is.null(sample_info) & !is.null(gene_info) |

      # Make exception for weirs cases (e.g. c(sample, counts))
      (colnames(count_info) %>% outersect(c(f_(.data)$name, s_(.data)$name)) %>% length() %>% gt(0))
    )
  ) {
    output_df =
      count_info %>%
      when(!is.null(sample_info) ~ (.) %>% left_join(sample_info, by=s_(.data)$name), ~ (.)) %>%
      when(!is.null(gene_info) ~ (.) %>% left_join(gene_info, by=f_(.data)$name), ~ (.)) %>%
      when(!is.null(range_info) ~ (.) %>% left_join(range_info, by=f_(.data)$name), ~ (.))
  }
  else if(!is.null(sample_info) ){
    output_df = sample_info
  }
  else if(!is.null(gene_info)){
    output_df = gene_info %>%

      # If present join GRanges
      when(!is.null(range_info) ~ (.) %>% left_join(range_info, by=f_(.data)$name), ~ (.))
  }

  output_df %>%

    # Cleanup
    select(any_of(output_colnames)) %>%
    suppressWarnings()

}


.as_tibble_optimised = function(x, skip_GRanges = FALSE, .subset = NULL,
                                .name_repair=c("check_unique", "unique", "universal", "minimal"),
                                rownames=pkgconfig::get_config("tibble::rownames", NULL)){

  .subset = enquo(.subset)

  sample_info <-
    colData(x)  %>%

    # If reserved column names are present add .x
    change_reserved_column_names(x) %>%

    # Convert to tibble
    tibble::as_tibble(rownames=s_(x)$name) %>%
    setNames(c(s_(x)$name, colnames(colData(x))))

  range_info <-
    skip_GRanges %>%
    when(
      (.) ~ tibble() %>% list,
      ~  get_special_datasets(x)
    ) %>%
    reduce(left_join, by="coordinate")

  gene_info <-
    rowData(x)  %>%

    # If reserved column names are present add .x
    change_reserved_column_names(x)%>%

    # Convert to tibble
    tibble::as_tibble(rownames=f_(x)$name) %>%
    setNames(c(f_(x)$name, colnames(rowData(x))))


  count_info <- get_count_datasets(x)

  # Return
  if(quo_is_null(.subset))

    # If I want to return all columns
    count_info %>%
    inner_join(sample_info, by=s_(x)$name) %>%
    inner_join(gene_info, by=f_(x)$name) %>%
    when(nrow(range_info) > 0 ~ (.) %>% left_join(range_info) %>% suppressMessages(), ~ (.))

  # This function outputs a tibble after subsetting the columns
  else subset_tibble_output(x, count_info, sample_info, gene_info, range_info, !!.subset)


}

#' @importFrom S4Vectors metadata
f_ =  function(x){
  # Check if old deprecated columns are used
  if("feature__" %in% names(metadata(x))) feature__ = metadata(x)$feature__
  return(feature__)
}

#' @importFrom S4Vectors metadata
s_ = function(x){
  if("sample__" %in% names(metadata(x))) sample__ = metadata(x)$sample__
  return(sample__)
}
