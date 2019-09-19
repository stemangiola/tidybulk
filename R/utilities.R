

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import dplyr
#' @import tidyr
#' @importFrom purrr as_mapper
#'
#' @param .x A tibble
#' @param .p A boolean
#' @param .f1 A function
#' @param .f2 A function
#'
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
  switch(.p %>% `!` %>% sum(1),
         as_mapper(.f1)(.x),
         if (.f2 %>% is.null %>% `!`)
           as_mapper(.f2)(.x)
         else
           .x)

}

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import dplyr
#' @import tidyr
#'
#' @param .x A tibble
#' @param .p1 A boolean
#' @param .p2 ELSE IF condition
#' @param .f1 A function
#' @param .f2 A function
#' @param .f3 A function
#'
#' @return A tibble
ifelse2_pipe = function(.x, .p1, .p2, .f1, .f2, .f3 = NULL) {
  # Nested switch
  switch(# First condition
    .p1 %>% `!` %>% sum(1),

    # First outcome
    as_mapper(.f1)(.x),
    switch(
      # Second condition
      .p2 %>% `!` %>% sum(1),

      # Second outcome
      as_mapper(.f2)(.x),

      # Third outcome - if there is not .f3 just return the original data frame
      if (.f3 %>% is.null %>% `!`)
        as_mapper(.f3)(.x)
      else
        .x
    ))
}

#' Get matrix from tibble
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @export
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  rownames = enquo(rownames)
  tbl %>%

    # Through warning if data frame is not numerical beside the rownames column (if present)
    ifelse_pipe(
      do_check &&
        tbl %>%
        # If rownames defined eliminate it from the data frame
        ifelse_pipe(!quo_is_null(rownames), ~ .x[, -1], ~ .x) %>%
        dplyr::summarise_all(class) %>%
        tidyr::gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
      ~ {
        warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
        .x
      }
    ) %>%
    as.data.frame() %>%

    # Deal with rownames column if present
    ifelse_pipe(
      !quo_is_null(rownames),
      ~ .x %>%
        magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
        select(-1)
    ) %>%

    # Convert to matrix
    as.matrix()
}

#' Check whether a numeric vector has been log transformed
#'
#' @param x A numeric vector
#' @param counts_column A character name of the count column
#'
error_if_log_transformed <- function(x, counts_column) {
  counts_column = enquo(counts_column)

  if (x %>% nrow %>% `>` (0))
    if (x %>% summarise(m = !!counts_column %>% max) %>% pull(m) < 50)
      stop(
        "The input was log transformed, this algorithm requires raw (un-normalised) read counts"
      )
}

#' Check whether there are duplicated genes/transcripts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param input.df A tibble of read counts
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
error_if_duplicated_genes <- function(input.df,
                                      sample_column = `sample`,
                                      transcript_column = `transcript`,
                                      counts_column = `read count`) {
  sample_column = enquo(sample_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  duplicates <-
    input.df %>%
    select(!!sample_column,!!transcript_column,!!counts_column) %>%
    distinct() %>%
    count(!!sample_column,!!transcript_column) %>%
    filter(n > 1) %>%
    arrange(n %>% desc())

  if (duplicates %>% nrow() > 0) {
    writeLines("Those are the duplicated genes")
    duplicates %>% print()
    stop(
      "Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding."
    )
  }

  input.df

}

#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param input.df A tibble of read counts
#' @param counts_column A character name of the read count column
error_if_counts_is_na = function(input.df, counts_column) {
  counts_column = enquo(counts_column)

  # Do the check
  if (input.df %>% filter(!!counts_column %>% is.na) %>% nrow %>% `>` (0))
    stop("You have NA values in your counts")

  # If all good return original data frame
  input.df
}

#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param input.df A tibble of read counts
#' @param list_input A list
#' @param expected_type A character string
#'
error_if_wrong_input = function(input.df, list_input, expected_type) {




  # Do the check
  if (
    list_input %>%
    map(~ .x %>% class() %>% `[` (1)) %>%
    unlist %>%
    equals(expected_type) %>%
    `!`
  )
    stop("You have passed the wrong argument to the function. Please check again.")

  # If all good return original data frame
  input.df
}


#' Formula parser
#'
#' @param fm A formula
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
  if (attr(terms(fm), "response") == 1)
    stop("The formula must be of the kind \"~ covariates\" ")
  else
    as.character(attr(terms(fm), "variables"))[-1]
}

#' Scale design matrix
#'
#' @param df A tibble
#' @param formula A formula
#'
#' @return A tibble
#'
#'
scale_design = function(df, formula) {
  df %>%
    setNames(c("sample_idx", "(Intercept)", parse_formula(formula))) %>%
    gather(cov, value,-sample_idx) %>%
    group_by(cov) %>%
    mutate(value = ifelse(
      !grepl("Intercept", cov) &
        length(union(c(0, 1), value)) != 2,
      scale(value),
      value
    )) %>%
    ungroup() %>%
    spread(cov, value) %>%
    arrange(as.integer(sample_idx)) %>%
    select(`(Intercept)`, one_of(parse_formula(formula)))
}

#' Add attribute to abject
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}

#' From rlang deprecated
#'
#' @param x An array
#' @param values An array
#' @param before A boolean
#'
prepend = function (x, values, before = 1)
{
  n <- length(x)
  stopifnot(before > 0 && before <= n)
  if (before == 1) {
    c(values, x)
  }
  else {
    c(x[1:(before - 1)], values, x[before:n])
  }
}

#' Add class to abject
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class = function(var, name) {

  class(var) <- prepend(class(var),name)

  var
}

#' Get column names either from user or from attributes
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param input.df A tibble
#' @param sample_column A character name of the sample column
#' @param transcript_column A character name of the gene/transcript name column
#' @param counts_column A character name of the read count column
#'
#' @return A list of column enquo or error
get_sample_transcript_counts = function(input.df, sample_column, transcript_column, counts_column){

  # If setted by the user, enquo those
  if(
    sample_column %>% quo_is_symbol() &
    transcript_column %>% quo_is_symbol() &
    counts_column %>% quo_is_symbol()
  )
    return(list(
      sample_column = sample_column,
      transcript_column = transcript_column,
      counts_column = counts_column
    ))

  # Otherwise check if attribute exists
  else {

    # If so, take them from the attribute
    if(input.df %>% attr("parameters") %>% is.null %>% `!`)
      return(list(
        sample_column = attr(input.df, "parameters")$sample_column,
        transcript_column = attr(input.df, "parameters")$transcript_column,
        counts_column = attr(input.df, "parameters")$counts_column
      ))
    # Else through error
    else
      stop("
        The fucntion does not know what your sample, transcript and counts columns are.\n
        You have to either enter those as symbols (e.g., `sample`), \n
        or use the funtion create_tt_from_tibble() to pass your column names that will be remembered.
      ")
  }
}

#' Get column names either from user or from attributes
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param input.df A tibble
#' @param elements_column A character name of the sample column
#' @param feature_column A character name of the gene/transcript name column
#' @param of_samples A boolean
#'
#' @return A list of column enquo or error
#'
get_elements_features = function(input.df, elements_column, feature_column, of_samples){

  # If setted by the user, enquo those
  if(
    elements_column %>% quo_is_symbol() &
    feature_column %>% quo_is_symbol()
  )
    return(list(
      elements_column = elements_column,
      feature_column = feature_column
    ))

  # Otherwise check if attribute exists
  else {

    # If so, take them from the attribute
    if(input.df %>% attr("parameters") %>% is.null %>% `!`)

      return(list(
        elements_column =  switch(
          of_samples %>% `!` %>% sum(1),
          attr(input.df, "parameters")$sample_column,
          attr(input.df, "parameters")$transcript_column
        ),
        feature_column = switch(
          of_samples %>% `!` %>% sum(1),
          attr(input.df, "parameters")$transcript_column,
          attr(input.df, "parameters")$sample_column
        )
      ))
    # Else through error
    else
      stop("
        The fucntion does not know what your elements (e.g., sample) and features (e.g., transcripts) are.\n
        You have to either enter those as symbols (e.g., `sample`), \n
        or use the funtion create_tt_from_tibble() to pass your column names that will be remembered.
      ")
  }
}

#' Get column names either from user or from attributes
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param input.df A tibble
#' @param elements_column A character name of the sample column
#'
#' @return A list of column enquo or error
get_elements = function(input.df, elements_column){

  # If setted by the user, enquo those
  if(
    elements_column %>% quo_is_symbol()
  )
    return(list(
      elements_column = elements_column
    ))

  # Otherwise check if attribute exists
  else {

    # If so, take them from the attribute
    if(input.df %>% attr("parameters") %>% is.null %>% `!`)

      return(list(
        elements_column =  switch(
          of_samples %>% `!` %>% sum(1),
          attr(input.df, "parameters")$sample_column,
          attr(input.df, "parameters")$transcript_column
        )
      ))
    # Else through error
    else
      stop("
        The fucntion does not know what your elements (e.g., sample) are.\n
        You have to either enter those as symbols (e.g., `sample`), \n
        or use the funtion create_tt_from_tibble() to pass your column names that will be remembered.
      ")
  }
}

#' Sub function of drop_redundant_elements_though_reduced_dimensions
#'
#' @param df A tibble
#'
#'
#' @return A tibble with pairs to drop
select_closest_pairs = function(df) {
  couples <- df %>% head(n = 0)

  while (df %>% nrow() > 0) {
    pair <- df %>%
      arrange(dist) %>%
      head(n = 1)
    couples <- couples %>% bind_rows(pair)
    df <- df %>%
      filter(
        !`sample 1` %in% (pair %>% select(1:2) %>% as.character()) &
          !`sample 2` %in% (pair %>% select(1:2) %>% as.character())
      )
  }

  couples

}
