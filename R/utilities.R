my_stop = function() {
  stop("
        You should call tidybulk library *after* tidyverse libraries.
        tidybulk says: The function does not know what your sample, transcript and counts columns are.
        You have to either enter those as arguments, or use the function tidybulk() to pass your column names that will be remembered.
      ")
}

#' This is a generalisation of ifelse that accepts an object and return an objects
#'
#' @keywords internal
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
  switch(.p %>% not() %>% sum(1),
         as_mapper(.f1)(.x),
         if (.f2 %>% is.null %>% not())
           as_mapper(.f2)(.x)
         else
           .x)

}

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @keywords internal
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
    .p1 %>% not() %>% sum(1),

    # First outcome
    as_mapper(.f1)(.x),
    switch(
      # Second condition
      .p2 %>% not() %>% sum(1),

      # Second outcome
      as_mapper(.f2)(.x),

      # Third outcome - if there is not .f3 just return the original data frame
      if (.f3 %>% is.null %>% not())
        as_mapper(.f3)(.x)
      else
        .x
    ))
}



#' Check whether a numeric vector has been log transformed
#'
#' @keywords internal
#'
#' @param x A numeric vector
#' @param .abundance A character name of the transcript/gene abundance column
#'
#' @return NA
error_if_log_transformed <- function(x, .abundance) {
  .abundance = enquo(.abundance)

  if (x %>% nrow %>% gt(0))
    if (x %>% summarise(m = !!.abundance %>% max) %>% pull(m) < 50)
      stop(
        "tidybulk says: The input was log transformed, this algorithm requires raw (un-scaled) read counts"
      )
}

#' Check whether there are duplicated genes/transcripts
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom utils capture.output
#'
#'
#' @param .data A tibble of read counts
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
error_if_duplicated_genes <- function(.data,
                                      .sample = `sample`,
                                      .transcript = `transcript`,
                                      .abundance = `read count`) {
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)

  duplicates <-
    distinct( .data, !!.sample,!!.transcript,!!.abundance) %>%
    count(!!.sample,!!.transcript) %>%
    filter(n > 1) %>%
    arrange(n %>% desc())

  if (duplicates %>% nrow() > 0) {
    message("Those are the duplicated genes")
    duplicates %>% capture.output() %>% paste0(collapse = "\n") %>% message()
    stop(
      "tidybulk says: Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding."
    )
  }

  .data

}

#' Check whether there are NA counts
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @param .data A tibble of read counts
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
#'
error_if_counts_is_na = function(.data, .abundance) {
  .abundance = enquo(.abundance)

  # Do the check
  if (.data %>% filter(!!.abundance %>% is.na) %>% nrow %>% gt(0))
    stop("tidybulk says: You have NA values in your counts")

  # If all good return original data frame
  .data
}

#' Check whether there are NA counts
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map
#'
#'
#' @param .data A tibble of read counts
#' @param list_input A list
#' @param expected_type A character string
#'
#' @return A tbl
#'
error_if_wrong_input = function(.data, list_input, expected_type) {




  # Do the check
  if (
    list_input %>%
    map(~ .x %>% class() %>% `[` (1)) %>%
    unlist %>%
    equals(expected_type) %>%
    not()
  )
    stop("tidybulk says: You have passed the wrong argument to the function. Please check again.")

  # If all good return original data frame
  .data
}


#' .formula parser
#'
#' @keywords internal
#'
#' @importFrom stats terms
#'
#' @param fm a formula
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
  if (attr(terms(fm), "response") == 1)
    stop("tidybulk says: The .formula must be of the kind \"~ covariates\" ")
  else
    as.character(attr(terms(fm), "variables"))[-1]
}

#' Formula parser with survival
#'
#' @param fm A formula
#'
#' @return A character vector
#'
#'
parse_formula_survival <- function(fm) {
  pars = as.character(attr(terms(fm), "variables"))[-1]

  response = NULL
  if(attr(terms(fm), "response") == 1) response = pars[1]
  covariates = ifelse(attr(terms(fm), "response") == 1, pars[-1], pars)

  list(
    response = response,
    covariates = covariates
  )
}

#' Scale design matrix
#'
#' @keywords internal
#'
#' @importFrom stats setNames
#' @importFrom stats cov
#'
#' @param df A tibble
#' @param .formula a formula
#'
#' @return A tibble
#'
#'
scale_design = function(df, .formula) {
  df %>%
    setNames(c("sample_idx", "(Intercept)", parse_formula(.formula))) %>%
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
    select(`(Intercept)`, one_of(parse_formula(.formula)))
}

get_tt_columns = function(.data){
  if(.data %>% attr("internals") %>% is.list() & "internals" %in% (.data %>% attr("internals") %>% names()))
    .data %>% attr("internals") %$% tt_columns
  else NULL
}

add_tt_columns = function(.data,
                          .sample,
                          .transcript,
                          .abundance,
                          .abundance_scaled = NULL,
													.abundance_adjusted = NULL){

  # Make col names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  .abundance_scaled = enquo(.abundance_scaled)
  .abundance_adjusted = enquo(.abundance_adjusted)

  # Add tt_columns
  .data %>% attach_to_internals(
     list(
      .sample = .sample,
      .transcript = .transcript,
      .abundance = .abundance
    ) %>%

    # If .abundance_scaled is not NULL add it to tt_columns
    ifelse_pipe(
      .abundance_scaled %>% quo_is_symbol,
      ~ .x %>% c(		list(.abundance_scaled = .abundance_scaled))
    ) %>%

  	ifelse_pipe(
  		.abundance_adjusted %>% quo_is_symbol,
  		~ .x %>% c(		list(.abundance_adjusted = .abundance_adjusted))
  	),
    "tt_columns"
  )

}

initialise_tt_internals = function(.data){
  .data %>%
    ifelse_pipe(
      "internals" %in% ((.) %>% attributes %>% names) %>% not(),
      ~ .x %>% add_attr(list(), "internals")
    )
}

reattach_internals = function(.data, .data_internals_from = NULL){
  if(.data_internals_from %>% is.null)
    .data_internals_from = .data

  .data %>% add_attr(.data_internals_from %>% attr("internals"), "internals")
}

attach_to_internals = function(.data, .object, .name){

  internals =
    .data %>%
    initialise_tt_internals() %>%
    attr("internals")

  # Add tt_bolumns
  internals[[.name]] = .object

  .data %>% add_attr(internals, "internals")
}

drop_internals = function(.data){

  .data %>% drop_attr("internals")
}

memorise_methods_used = function(.data, .method, object_containing_methods = .data){
  # object_containing_methods is used in case for example of test gene rank where the output is a tibble that has lost its attributes

	.data %>%
		attach_to_internals(
		  object_containing_methods %>%
				attr("internals") %>%
				.[["methods_used"]] %>%
				c(.method) %>%	unique(),
			"methods_used"
		)

}

#' Add attribute to abject
#'
#' @keywords internal
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

#' Drop attribute to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
drop_attr = function(var, name) {
  attr(var, name) <- NULL
  var
}

#' Remove class to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
#' @keywords internal
drop_class = function(var, name) {
  class(var) <- class(var)[!class(var)%in%name]
  var
}

#' From rlang deprecated
#'
#' @keywords internal
#'
#' @param x An array
#' @param values An array
#' @param before A boolean
#'
#' @return An array
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
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class = function(var, name) {

  if(!name %in% class(var)) class(var) <- prepend(class(var),name)

  var
}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#'
#' @return A list of column enquo or error
get_sample_transcript_counts = function(.data, .sample, .transcript, .abundance){

    if( .sample %>% quo_is_symbol() ) .sample = .sample
    else if(".sample" %in% (.data %>% get_tt_columns() %>% names))
      .sample =  get_tt_columns(.data)$.sample
    else my_stop()

    if( .transcript %>% quo_is_symbol() ) .transcript = .transcript
    else if(".transcript" %in% (.data %>% get_tt_columns() %>% names))
      .transcript =  get_tt_columns(.data)$.transcript
    else my_stop()

    if( .abundance %>% quo_is_symbol() ) .abundance = .abundance
    else if(".abundance" %in% (.data %>% get_tt_columns() %>% names))
      .abundance = get_tt_columns(.data)$.abundance
    else my_stop()

    list(.sample = .sample, .transcript = .transcript, .abundance = .abundance)

}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .abundance A character name of the read count column
#'
#' @return A list of column enquo or error
get_sample_counts = function(.data, .sample, .abundance){

  if( .sample %>% quo_is_symbol() ) .sample = .sample
  else if(".sample" %in% (.data %>% get_tt_columns() %>% names))
    .sample =  get_tt_columns(.data)$.sample
  else my_stop()

  if( .abundance %>% quo_is_symbol() ) .abundance = .abundance
  else if(".abundance" %in% (.data %>% get_tt_columns() %>% names))
    .abundance = get_tt_columns(.data)$.abundance
  else my_stop()

  list(.sample = .sample, .abundance = .abundance)

}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#'
#' @return A list of column enquo or error
get_sample = function(.data, .sample){

  if( .sample %>% quo_is_symbol() ) .sample = .sample
  else if(".sample" %in% (.data %>% get_tt_columns() %>% names))
    .sample =  get_tt_columns(.data)$.sample
  else my_stop()

  list(.sample = .sample)

}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .transcript A character name of the transcript column
#'
#' @return A list of column enquo or error
get_transcript = function(.data, .transcript){


  my_stop = function() {
    stop("
        tidybulk says: The function does not know what your transcript, transcript and counts columns are.\n
        You have to either enter those as symbols (e.g., `transcript`), \n
        or use the function create_tt_from_tibble() to pass your column names that will be remembered.
      ")
  }

  if( .transcript %>% quo_is_symbol() ) .transcript = .transcript
  else if(".transcript" %in% (.data %>% get_tt_columns() %>% names))
    .transcript =  get_tt_columns(.data)$.transcript
  else my_stop()

  list(.transcript = .transcript)

}


#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#'
#' @return A list of column enquo or error
get_sample_transcript = function(.data, .sample, .transcript){

  if( .sample %>% quo_is_symbol() ) .sample = .sample
  else if(".sample" %in% (.data %>% get_tt_columns() %>% names))
    .sample =  get_tt_columns(.data)$.sample
  else my_stop()

  if( .transcript %>% quo_is_symbol() ) .transcript = .transcript
  else if(".transcript" %in% (.data %>% get_tt_columns() %>% names))
    .transcript =  get_tt_columns(.data)$.transcript
  else my_stop()


  list(.sample = .sample, .transcript = .transcript)

}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#'
#' @return A list of column enquo or error
get_sample = function(.data, .sample){

  if( .sample %>% quo_is_symbol() ) .sample = .sample
  else if(".sample" %in% (.data %>% get_tt_columns() %>% names))
    .sample =  get_tt_columns(.data)$.sample
  else my_stop()

  list(.sample = .sample)

}



#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .element A character name of the sample column
#' @param .feature A character name of the transcript/gene column
#' @param of_samples A boolean
#'
#' @return A list of column enquo or error
#'
get_elements_features = function(.data, .element, .feature, of_samples = TRUE){

  # If setted by the user, enquo those
  if(
    .element %>% quo_is_symbol() &
    .feature %>% quo_is_symbol()
  )
    return(list(
      .element = .element,
      .feature = .feature
    ))

  # Otherwise check if attribute exists
  else {

    # If so, take them from the attribute
    if(.data %>% get_tt_columns() %>% is.null %>% not())

      return(list(
        .element =  switch(
          of_samples %>% not() %>% sum(1),
          get_tt_columns(.data)$.sample,
          get_tt_columns(.data)$.transcript
        ),
        .feature = switch(
          of_samples %>% not() %>% sum(1),
          get_tt_columns(.data)$.transcript,
          get_tt_columns(.data)$.sample
        )
      ))
    # Else through error
    else
      stop("
        tidybulk says: The function does not know what your elements (e.g., sample) and features (e.g., transcripts) are.\n
        You have to either enter those as symbols (e.g., `sample`), \n
        or use the function create_tt_from_tibble() to pass your column names that will be remembered.
      ")
  }
}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .element A character name of the sample column
#' @param .feature A character name of the transcript/gene column
#' @param .abundance A character name of the read count column

#' @param of_samples A boolean
#'
#' @return A list of column enquo or error
#'
get_elements_features_abundance = function(.data, .element, .feature, .abundance, of_samples = TRUE){

  my_stop = function() {
    stop("
        tidybulk says: The function does not know what your elements (e.g., sample) and features (e.g., transcripts) are.\n
        You have to either enter those as symbols (e.g., `sample`), \n
        or use the function create_tt_from_tibble() to pass your column names that will be remembered.
      ")
  }

  if( .element %>% quo_is_symbol() ) .element = .element
  else if(of_samples & ".sample" %in% (.data %>% get_tt_columns() %>% names))
    .element =  get_tt_columns(.data)$.sample
  else if((!of_samples) & ".transcript" %in% (.data %>% get_tt_columns() %>% names))
     .element =  get_tt_columns(.data)$.transcript
  else my_stop()

  if( .feature %>% quo_is_symbol() ) .feature = .feature
  else if(of_samples & ".transcript" %in% (.data %>% get_tt_columns() %>% names))
    .feature =  get_tt_columns(.data)$.transcript
  else if((!of_samples) & ".sample" %in% (.data %>% get_tt_columns() %>% names))
    .feature =  get_tt_columns(.data)$.sample
  else my_stop()

  if( .abundance %>% quo_is_symbol() ) .abundance = .abundance
  else if(".abundance" %in% (.data %>% get_tt_columns() %>% names))
    .abundance = get_tt_columns(.data)$.abundance
  else my_stop()

  list(.element = .element, .feature = .feature, .abundance = .abundance)
}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .element A character name of the sample column
#' @param of_samples A boolean
#'
#' @return A list of column enquo or error
get_elements = function(.data, .element, of_samples = TRUE){

  # If setted by the user, enquo those
  if(
    .element %>% quo_is_symbol()
  )
    return(list(
      .element = .element
    ))

  # Otherwise check if attribute exists
  else {

    # If so, take them from the attribute
    if(.data %>% get_tt_columns() %>% is.null %>% not())

      return(list(
        .element =  switch(
          of_samples %>% not() %>% sum(1),
          get_tt_columns(.data)$.sample,
          get_tt_columns(.data)$.transcript
        )
      ))
    # Else through error
    else
      stop("
        tidybulk says: The function does not know what your elements (e.g., sample) are.\n
        You have to either enter those as symbols (e.g., `sample`), \n
        or use the function create_tt_from_tibble() to pass your column names that will be remembered.
      ")
  }
}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .abundance A character name of the abundance column
#'
#' @return A list of column enquo or error
get_abundance_norm_if_exists = function(.data, .abundance){

  # If setted by the user, enquo those
  if(
    .abundance %>% quo_is_symbol()
  )
    return(list(
      .abundance = .abundance
    ))

  # Otherwise check if attribute exists
  else {

    # If so, take them from the attribute
    if(.data %>% get_tt_columns() %>% is.null %>% not())

      return(list(
        .abundance =  switch(
          (".abundance_scaled" %in% (.data %>% get_tt_columns() %>% names) &&
             # .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% not() &&
             quo_name(.data %>% get_tt_columns() %$% .abundance_scaled) %in% (.data %>% colnames)
           ) %>% not() %>% sum(1),
          get_tt_columns(.data)$.abundance_scaled,
          get_tt_columns(.data)$.abundance
        )
      ))
    # Else through error
    else
      stop("
        tidybulk says: The function does not know what your elements (e.g., sample) are.\n
        You have to either enter those as symbols (e.g., `sample`), \n
        or use the function create_tt_from_tibble() to pass your column names that will be remembered.
      ")
  }
}

#' Sub function of remove_redundancy_elements_though_reduced_dimensions
#'
#' @keywords internal
#'
#' @importFrom stats dist
#' @importFrom utils head
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

#' This function is needed for DE in case the matrix is not rectangular, but includes NA
#'
#' @keywords internal
#'
#' @param .matrix A matrix
#'
#' @return A matrix
fill_NA_with_row_median = function(.matrix){

  if(length(which(rowSums(is.na(.matrix)) > 0)) > 0)
    rbind(
      .matrix[rowSums(is.na(.matrix)) == 0,],
      apply(.matrix[rowSums(is.na(.matrix)) > 0,], 1, FUN = function(.x) { .x[is.na(.x)] = median(.x, na.rm = TRUE); .x}) %>% t
    )
  else
    .matrix
}

# #' @importFrom magrittr %>%
# #' @export
# magrittr::`%>%`

#' @importFrom tibble tibble
#' @export
tibble::tibble

#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble

#' get_x_y_annotation_columns
#'
#' @keywords internal
#'
#' @importFrom magrittr equals
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .horizontal The name of the column horizontally presented in the heatmap
#' @param .vertical The name of the column vertically presented in the heatmap
#' @param .abundance The name of the transcript/gene abundance column
#' @param .abundance_scaled The name of the transcript/gene scaled abundance column
#'
#' @description This function recognise what are the sample-wise columns and transcrip-wise columns
#'
#' @return A list
#'
get_x_y_annotation_columns = function(.data, .horizontal, .vertical, .abundance, .abundance_scaled){


  # Comply with CRAN NOTES
  . = NULL

  # Make col names
  .horizontal = enquo(.horizontal)
  .vertical = enquo(.vertical)
  .abundance = enquo(.abundance)
  .abundance_scaled = enquo(.abundance_scaled)

  # x-annotation df
  n_x = .data %>% distinct(!!.horizontal) %>% nrow
  n_y = .data %>% distinct(!!.vertical) %>% nrow

  # Sample wise columns
  horizontal_cols=
    .data %>%
    select(-!!.horizontal, -!!.vertical, -!!.abundance) %>%
    colnames %>%
    map(
      ~
        .x %>%
        ifelse_pipe(
          .data %>%
            distinct(!!.horizontal, !!as.symbol(.x)) %>%
            nrow %>%
            equals(n_x),
          ~ .x,
          ~ NULL
        )
    ) %>%

    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist

  # Transcript wise columns
  vertical_cols=
    .data %>%
    select(-!!.horizontal, -!!.vertical, -!!.abundance, -horizontal_cols) %>%
    colnames %>%
    map(
      ~
        .x %>%
        ifelse_pipe(
          .data %>%
            distinct(!!.vertical, !!as.symbol(.x)) %>%
            nrow %>%
            equals(n_y),
          ~ .x,
          ~ NULL
        )
    ) %>%

    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist

  # Counts wise columns, at the moment scaled counts is treated as special and not accounted for here
  counts_cols =
    .data %>%
    select(-!!.horizontal, -!!.vertical, -!!.abundance) %>%

    # Exclude horizontal
    ifelse_pipe(!is.null(horizontal_cols),  ~ .x %>% select(-horizontal_cols)) %>%

    # Exclude vertical
    ifelse_pipe(!is.null(vertical_cols),  ~ .x %>% select(-vertical_cols)) %>%

    # Exclude scaled counts if exist
    ifelse_pipe(.abundance_scaled %>% quo_is_symbol,  ~ .x %>% select(-!!.abundance_scaled) ) %>%

    # Select colnames
    colnames %>%

    # select columns
    map(
      ~
        .x %>%
        ifelse_pipe(
          .data %>%
            distinct(!!.vertical, !!.horizontal, !!as.symbol(.x)) %>%
            nrow %>%
            equals(n_x * n_y),
          ~ .x,
          ~ NULL
        )
    ) %>%

    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist

  list(  horizontal_cols = horizontal_cols,  vertical_cols = vertical_cols, counts_cols = counts_cols )
}

get_specific_annotation_columns = function(.data, .col){


  # Comply with CRAN NOTES
  . = NULL

  # Make col names
  .col = enquo(.col)

  # x-annotation df
  n_x = .data %>% distinct(!!.col) %>% nrow

  # Sample wise columns
  .data %>%
  select(-!!.col) %>%
  colnames %>%
  map(
    ~
      .x %>%
      ifelse_pipe(
        .data %>%
          distinct(!!.col, !!as.symbol(.x)) %>%
          nrow %>%
          equals(n_x),
        ~ .x,
        ~ NULL
      )
  ) %>%

  # Drop NULL
  {	(.)[lengths((.)) != 0]	} %>%
  unlist

}

#' log10_reverse_trans
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description it perform log scaling and reverse the axis. Useful to plot negative log probabilities. To not be used directly but with ggplot (e.g. scale_y_continuous(trans = "log10_reverse") )
#'
#' @importFrom scales trans_new
#' @importFrom scales log_breaks
#'
#' @return A scales object
#'
#' @examples
#'
#' library(ggplot2)
#' library(tibble)
#'
#' tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) %>%
#'  ggplot(aes(fold_change , pvalue)) +
#'  geom_point() +
#'  scale_y_continuous(trans = "log10_reverse")
#'
#' @export
log10_reverse_trans <- function() {
	trans <- function(x) -log10(x)
	inv <- function(x) 10^(-x)

	trans_new("log10_reverse", trans, inv, log_breaks(base = 10))
}

#' logit scale
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description it perform logit scaling with right axis formatting. To not be used directly but with ggplot (e.g. scale_y_continuous(trans = "log10_reverse") )
#'
#' @importFrom scales label_scientific
#' @importFrom scales extended_breaks
#' @importFrom stats qlogis plogis
#'
#' @return A scales object
#'
#' @examples
#'
#' library(ggplot2)
#' library(tibble)
#'
#' tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) %>%
#'  ggplot(aes(fold_change , pvalue)) +
#'  geom_point() +
#'  scale_y_continuous(trans = "log10_reverse")
#'
#' @export
logit_trans <- function(){


	if (find.package("functional", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing functional needed for analyses")
		install.packages("functional", repos = "https://cloud.r-project.org")
	}

	trans <- qlogis
	inv <- plogis

	trans_new("logit",
						transform = trans,
						inverse = inv,
						breaks = functional::Compose(trans, extended_breaks(), inv),
						format = label_scientific(digits = 2)
	)
}


#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {

	v = quo_name(quo_squash(v))
	gsub('^c\\(|`|\\)$', '', v) %>%
		strsplit(', ') %>%
		unlist
}

#' Drop lowly transcribed genes for TMM normalization
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats median
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param factor_of_interest The name of the column of the factor of interest
#' @param minimum_counts A positive integer. Minimum counts required for at least some samples.
#' @param minimum_proportion A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#'
#' @return A tibble filtered
add_scaled_counts_bulk.get_low_expressed <- function(.data,
																										 .sample = `sample`,
																										 .transcript = `transcript`,
																										 .abundance = `count`,
																										 factor_of_interest = NULL,
																										 minimum_counts = 10,
																										 minimum_proportion = 0.7) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	factor_of_interest = enquo(factor_of_interest)

	# Check if package is installed, otherwise install
	if (find.package("edgeR", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing edgeR needed for differential transcript abundance analyses")
		if (!requireNamespace("BiocManager", quietly = TRUE))
			install.packages("BiocManager", repos = "https://cloud.r-project.org")
		BiocManager::install("edgeR", ask = FALSE)
	}

	# Check if factor_of_interest is continuous and exists
	string_factor_of_interest =

		factor_of_interest %>%
		ifelse2_pipe(
			quo_is_symbol(.) &&
				select(.data, !!(.)) %>% lapply(class) %>%	as.character() %in% c("numeric", "integer", "double"),
			quo_is_symbol(.),
			~ {
				message("tidybulk says: The factor of interest is continuous (e.g., integer,numeric, double). The data will be filtered without grouping.")
				NULL
			},
			~ .data %>%
				distinct(!!.sample, !!factor_of_interest) %>%
				arrange(!!.sample) %>%
				pull(!!factor_of_interest),
			~ NULL
		)

	if (minimum_counts < 0)
		stop("The parameter minimum_counts must be > 0")
	if (minimum_proportion < 0 |	minimum_proportion > 1)
		stop("The parameter minimum_proportion must be between 0 and 1")

	.data %>%
		select(!!.sample,!!.transcript, !!.abundance) %>%
		spread(!!.sample, !!.abundance) %>%

		# Drop if transcript have missing value
		drop_na() %>%

		# If I don't have any transcript with all samples give meaningful error
		when(
			nrow(.) == 0 ~ stop("tidybulk says: you don't have any transcript that is in all samples. Please consider using impute_missing_abundance."),
			~ (.)
		) %>%

		# Call edgeR
		as_matrix(rownames = !!.transcript) %>%
		edgeR::filterByExpr(
			min.count = minimum_counts,
			group = string_factor_of_interest,
			min.prop = minimum_proportion
		) %>%
		not() %>%
		which %>%
		names %>%

		# Attach attributes
		reattach_internals(.data)
}

#' Calculate the norm factor with calcNormFactor from limma
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#'
#' @param .data A tibble
#' @param reference A reference matrix, not sure if used anymore
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#'
#'
#' @return A list including the filtered data frame and the normalization factors
add_scaled_counts_bulk.calcNormFactor <- function(.data,
																									reference = NULL,
																									.sample = `sample`,
																									.transcript = `transcript`,
																									.abundance = `count`,
																									method) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	error_if_log_transformed(.data,!!.abundance)

	# Get data frame for the highly transcribed transcripts
	df.filt <-
		.data %>%
		# dplyr::filter(!(!!.transcript %in% gene_to_exclude)) %>%
		droplevels() %>%
		select(!!.sample, !!.transcript, !!.abundance)

	df.filt.spread =
	  df.filt %>%
	  tidyr::spread(!!.sample,!!.abundance) %>%
	  tidyr::drop_na() %>%
	  dplyr::select(-!!.transcript)

	# If not enough genes, warning
	if(nrow(df.filt.spread)<100) warning(warning_for_scaling_with_few_genes)

	# scaled data set
	nf =
		tibble::tibble(
			# Sample factor
			sample = factor(levels(df.filt %>% pull(!!.sample))),

			# scaled data frame
			nf = edgeR::calcNormFactors(
			  df.filt.spread,
				refColumn = which(reference == factor(levels(
					df.filt %>% pull(!!.sample)
				))),
				method = method
			)
		) %>%

		setNames(c(quo_name(.sample), "nf")) %>%

		# Add the statistics about the number of genes filtered
		dplyr::left_join(
			df.filt %>%
				dplyr::group_by(!!.sample) %>%
				dplyr::summarise(tot_filt = sum(!!.abundance, na.rm = TRUE)) %>%
				dplyr::mutate(!!.sample := as.factor(as.character(!!.sample))),
			by = quo_name(.sample)
		)

	# Return
	list(
		# gene_to_exclude = gene_to_exclude,
		nf = nf
	) %>%

		# Attach attributes
		reattach_internals(.data)
}

# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

# Raise to the power
pow = function(a,b){	a^b }

outersect <- function(x, y) {
	sort(c(setdiff(x, y),
				 setdiff(y, x)))
}

do_validate = function(){

	if(!"tidybulk_do_validate" %in% names(options())) TRUE
	else getOption("tidybulk_do_validate")

}

#' @importFrom stringr str_remove
#' @importFrom stringr str_replace_all
#'
multivariable_differential_tissue_composition = function(
	deconvoluted,
	method,
	.my_formula,
	min_detected_proportion
){
	results_regression =
		deconvoluted %>%

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
				if (find.package("survival", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing betareg needed for analyses")
					install.packages("survival", repos = "https://cloud.r-project.org")
				}

				if (find.package("boot", quiet = TRUE) %>% length %>% equals(0)) {
					message("Installing boot needed for analyses")
					install.packages("boot", repos = "https://cloud.r-project.org")
				}

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

univariable_differential_tissue_composition = function(
	deconvoluted,
	method,
	.my_formula,
	min_detected_proportion
){
		deconvoluted %>%

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
							if (find.package("survival", quiet = TRUE) %>% length %>% equals(0)) {
								message("Installing betareg needed for analyses")
								install.packages("survival", repos = "https://cloud.r-project.org")
							}

							if (find.package("boot", quiet = TRUE) %>% length %>% equals(0)) {
								message("Installing boot needed for analyses")
								install.packages("boot", repos = "https://cloud.r-project.org")
							}

							(.) %>%
								mutate(.proportion_0_corrected = .proportion_0_corrected  %>% boot::logit()) %>%
								survival::coxph(.my_formula, .)	%>%
								broom::tidy() %>%
								select(-term)
						} ,
						~ {
							# Check if package is installed, otherwise install
							if (find.package("betareg", quiet = TRUE) %>% length %>% equals(0)) {
								message("Installing betareg needed for analyses")
								install.packages("betareg", repos = "https://cloud.r-project.org")
							}
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

univariable_differential_tissue_stratification = function(
	deconvoluted,
	method,
	.my_formula
){

	# Check if package is installed, otherwise install
	if (find.package("survival", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing survival needed for analyses")
		install.packages("survival", repos = "https://cloud.r-project.org")
	}

	# Check if package is installed, otherwise install
	if (find.package("survminer", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing survminer needed for analyses")
		install.packages("survminer", repos = "https://cloud.r-project.org")
	}


	if (find.package("broom", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing broom needed for analyses")
		install.packages("broom", repos = "https://cloud.r-project.org")
	}

	deconvoluted %>%

		# Test
		pivot_longer(
			names_prefix = sprintf("%s: ", method),
			cols = starts_with(method),
			names_to = ".cell_type",
			values_to = ".proportion"
		) %>%

		# Test survival
		tidyr::nest(cell_type_proportions = -.cell_type) %>%
		mutate(surv_test = map(
			cell_type_proportions,
			~ {

				data = .x %>%
					mutate(.high_cellularity = .proportion > median(.proportion))

				if(data %>%
					 distinct(.high_cellularity) %>%
					 nrow %>%
					 equals(1)
				) return(NULL)

				# See if regression if censored or not
				fit = survival::survdiff(data = data, .my_formula)

				p =
					survminer::surv_fit(data = data, .my_formula) %>%
					survminer::ggsurvplot(
						fit=.,
						data = data,
						risk.table = FALSE,
						conf.int = T,
						palette = c("#ed6f68",  "#5366A0" ),
						legend = "none",
						pval = T
					)

				fit %>%
					broom::tidy() %>%
					select(-N, -obs) %>%
					spread(.high_cellularity, exp) %>%
					setNames(c(".low_cellularity_expected", ".high_cellularity_expected")) %>%
					mutate(pvalue = 1 - pchisq(fit$chisq, length(fit$n) - 1)) %>%
					mutate(plot = list(p))

			}
		)) %>%

		unnest(surv_test, keep_empty = TRUE)
}

univariable_differential_tissue_stratification_SE = function(
	deconvoluted,
	method,
	.my_formula
){

	# Check if package is installed, otherwise install
	if (find.package("survival", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing survival needed for analyses")
		install.packages("survival", repos = "https://cloud.r-project.org")
	}

	# Check if package is installed, otherwise install
	if (find.package("survminer", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing survminer needed for analyses")
		install.packages("survminer", repos = "https://cloud.r-project.org")
	}


	if (find.package("broom", quiet = TRUE) %>% length %>% equals(0)) {
		message("Installing broom needed for analyses")
		install.packages("broom", repos = "https://cloud.r-project.org")
	}

	deconvoluted %>%

		pivot_sample() %>%

		# Test
		pivot_longer(
			names_prefix = sprintf("%s: ", method),
			cols = starts_with(method),
			names_to = ".cell_type",
			values_to = ".proportion"
		) %>%

		# Test survival
		tidyr::nest(cell_type_proportions = -.cell_type) %>%
		mutate(surv_test = map(
			cell_type_proportions,
			~ {

				data = .x %>%
					mutate(.high_cellularity = .proportion > median(.proportion))

				if(data %>%
					 distinct(.high_cellularity) %>%
					 nrow %>%
					 equals(1)
				) return(NULL)

				# See if regression if censored or not
				fit = survival::survdiff(data = data, .my_formula)

				p =
					survminer::surv_fit(data = data, .my_formula) %>%
					survminer::ggsurvplot(
						fit=.,
						data = data,
						risk.table = FALSE,
						conf.int = T,
						palette = c("#ed6f68",  "#5366A0" ),
						legend = "none",
						pval = T
					)

				fit %>%
					broom::tidy() %>%
					select(-N, -obs) %>%
					spread(.high_cellularity, exp) %>%
					setNames(c(".low_cellularity_expected", ".high_cellularity_expected")) %>%
					mutate(pvalue = 1 - pchisq(fit$chisq, length(fit$n) - 1)) %>%
					mutate(plot = list(p))

			}
		)) %>%

		unnest(surv_test, keep_empty = TRUE)
}

# Function that rotates a 2D space of a arbitrary angle
rotation = function(m, d) {
	r = d * pi / 180
	((dplyr::bind_rows(
		c(`1` = cos(r), `2` = -sin(r)),
		c(`1` = sin(r), `2` = cos(r))
	) %>% as_matrix) %*% m)
}

#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
get_special_datasets <- function(SummarizedExperiment_object) {
	if (
		"RangedSummarizedExperiment" %in% .class2(SummarizedExperiment_object) &

		rowRanges(SummarizedExperiment_object) %>%
		as.data.frame() %>%
		nrow() %>%
		gt(0)
	) {
		rowRanges(SummarizedExperiment_object) %>%
			as.data.frame() %>%

			# Take off rowData columns as there is a recursive anomaly within gene ranges
			suppressWarnings(
				select(-one_of(colnames(rowData(SummarizedExperiment_object))))
			) %>%
			tibble::as_tibble(rownames="feature") %>%
			list()
	} else {
		tibble() %>% list()
	}
}

combineByRow <- function(m, fun = NULL) {
  # Shown here
  #https://stackoverflow.com/questions/8139301/aggregate-rows-in-a-large-matrix-by-rowname

  m <- m[ order(rownames(m)), ]

  ## keep track of previous row name
  prev <- rownames(m)[1]
  i.start <- 1
  i.end <- 1

  ## cache the rownames -- profiling shows that it takes
  ## forever to look at them
  m.rownames <- rownames(m)
  stopifnot(all(!is.na(m.rownames)))


  ## go through matrix in a loop, as we need to combine some unknown
  ## set of rows
  for (i in 2:(1+nrow(m))) {

    curr <- m.rownames[i]

    ## if we found a new row name (or are at the end of the matrix),
    ## combine all rows and mark invalid rows
    if (prev != curr || is.na(curr)) {
      if (i.start < i.end) {
        m[i.start,] <- apply(m[i.start:i.end,], 2, fun)
        m.rownames[(1+i.start):i.end] <- NA
      }

      prev <- curr
      i.start <- i
    } else {
      i.end <- i
    }
  }

  m[ which(!is.na(m.rownames)),]
}

filter_genes_on_condition = function(.data, .subset_for_scaling){

  .subset_for_scaling = enquo(.subset_for_scaling)

  # Get genes from condition
  my_genes =
    rowData(.data) %>%
    as_tibble(rownames=".feature") %>%
    filter(!!.subset_for_scaling) %>%
    pull(.feature)

  .data[rownames(.data) %in% my_genes,]

}


which_NA_matrix = function(.data){
  is_na <- which(is.na(.data), arr.ind=TRUE)
  which_is_NA = fill_matrix_with_FALSE(.data )
  which_is_NA[is_na] <- TRUE
  which_is_NA
}

fill_matrix_with_FALSE = function(.data){
  .data[,] = FALSE
  .data
}

rowMedians = function(.data, na.rm){
  apply(.data, 1, median, na.rm=na.rm)
}

fill_NA_matrix_with_factor_colwise = function(.data, factor){

  rn = rownames(.data)
  cn = colnames(.data)

  .data %>%
    t %>%
    split.data.frame(factor) %>%
    map(~ t(.x)) %>%

    # Fill
    map(
      ~ {
        k <- which(is.na(.x), arr.ind=TRUE)
        .x[k] <- rowMedians(.x, na.rm=TRUE)[k[,1]]
        .x
      }
    ) %>%

    # Add NA factors if any
    when(
      is.na(factor) %>% length() %>% gt(0) ~ (.) %>% c(list(.data[,is.na(factor)])),
      ~ (.)
    ) %>%

    # Merge
    reduce(cbind) %>%

    # Reorder rows and column as it was
    .[rn, cn]

}
