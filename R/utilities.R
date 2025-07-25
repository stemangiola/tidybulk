


#' Internal error function for tidybulk
#'
#' @keywords internal
#' @noRd
#'
#' @return Stops execution with a helpful error message
my_stop = function() {
  stop("
        tidybulk says: The function does not know what your sample, transcript and counts columns are.
	You might need to specify the arguments .sample, .transcript and/or .abundance. 
	Please read the documentation of this function for more information.
      ")
}





#' Check whether there are NA counts
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
#' @import tibble
#'
#' @param .data A tibble of read counts
#' @param abundance A character name of the read count column
#' @param .abundance DEPRECATED
#'
#' @return A tbl
#'
error_if_counts_is_na = function(.data, abundance, .abundance = NULL) {
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "error_if_counts_is_na(.abundance)", "error_if_counts_is_na(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }

  # Do the check
  if (.data %>% filter(!!abundance %>% is.na) %>% nrow() %>% gt(0))
    stop("tidybulk says: You have NA values in your counts")

  # If all good return original data frame
  .data
}

#' Check whether there are NA counts
#'
#' @keywords internal
#' @noRd
#'
#' 
#' 
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
#' @noRd
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



#' Scale design matrix
#'
#' @keywords internal
#' @noRd
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
    select(`(Intercept)`, any_of(parse_formula(.formula)))
}

get_tt_columns = function(.data){
  if(
    .data %>% attr("internals") %>% is.list() &&
    "tt_columns" %in% names(.data %>% attr("internals"))
    ) #& "internals" %in% (.data %>% attr("internals") %>% names()))
    .data %>% attr("internals") %$% tt_columns
  else NULL
}

#' @importFrom rlang quo_is_symbol
#'
add_tt_columns = function(.data,
                          .sample,
                          .transcript,
                          abundance,
                          .abundance = NULL,
													.abundance_scaled = NULL,
													.abundance_adjusted = NULL){

  # Make col names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  .abundance_scaled = enquo(.abundance_scaled)
  .abundance_adjusted = enquo(.abundance_adjusted)
  tt_list <- list(
    .sample = .sample,
    .transcript = .transcript,
    .abundance = .abundance
  )
  if (.abundance_scaled %>% quo_is_symbol) {
    tt_list <- c(tt_list, list(.abundance_scaled = .abundance_scaled))
  }
  if (.abundance_adjusted %>% quo_is_symbol) {
    tt_list <- c(tt_list, list(.abundance_adjusted = .abundance_adjusted))
  }
  .data %>% attach_to_internals(tt_list, "tt_columns")

}

initialise_tt_internals = function(.data){
  if ("internals" %in% (attributes(.data) %>% names) %>% not()) {
    .data <- add_attr(.data, list(), "internals")
  }
  .data
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
#' @noRd
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
#' @noRd
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
#' @noRd
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
drop_class = function(var, name) {
  class(var) <- class(var)[!class(var)%in%name]
  var
}

#' From rlang deprecated
#'
#' @keywords internal
#' @noRd
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
#' @noRd
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
#' @noRd
#'
#' @importFrom rlang quo_is_symbol
#' @importFrom rlang quo_is_symbolic
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param abundance A character name of the read count column
#' @param .abundance DEPRECATED
#'
#' @return A list of column enquo or error
get_sample_transcript_counts = function(.data, .sample, .transcript, abundance, .abundance = NULL){
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_sample_transcript_counts(.abundance)", "get_sample_transcript_counts(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }

    if( quo_is_symbolic(.sample) ) .sample = .sample
    else if(".sample" %in% (.data %>% get_tt_columns() %>% names))
      .sample =  get_tt_columns(.data)$.sample
    else my_stop()

    if( quo_is_symbolic(.transcript) ) .transcript = .transcript
    else if(".transcript" %in% (.data %>% get_tt_columns() %>% names))
      .transcript =  get_tt_columns(.data)$.transcript
    else my_stop()

    if(  quo_is_symbolic(abundance) ) abundance = abundance
    else if(".abundance" %in% (.data %>% get_tt_columns() %>% names))
      abundance = get_tt_columns(.data)$.abundance
    else my_stop()

    list(.sample = .sample, .transcript = .transcript, .abundance = abundance)

}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param abundance A character name of the read count column
#' @param .abundance DEPRECATED
#'
#' @return A list of column enquo or error
get_sample_counts = function(.data, .sample, abundance, .abundance = NULL){
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_sample_counts(.abundance)", "get_sample_counts(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }

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
#' @noRd
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
#' @noRd
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .transcript A character name of the transcript column
#'
#' @return A list of column enquo or error
get_transcript = function(.data, .transcript){




  if( .transcript %>% quo_is_symbol() ) .transcript = .transcript
  else if(".transcript" %in% (.data %>% get_tt_columns() %>% names))
    .transcript =  get_tt_columns(.data)$.transcript
  else my_stop()

  list(.transcript = .transcript)

}


#' Get column names either from user or from attributes
#'
#' @keywords internal
#' @noRd
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
#' @noRd
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
#' @noRd
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
        tidybulk says: The function does not know what your sample, transcript and counts columns are.
	You might need to specify the arguments .sample, .transcript and/or .abundance 
	Please read the documentation of this function for more information.
      ")
  }
}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .element A character name of the sample column
#' @param .feature A character name of the transcript/gene column
#' @param abundance A character name of the read count column
#' @param .abundance DEPRECATED
#'
#' @return A list of column enquo or error
#'
get_elements_features_abundance = function(.data, .element, .feature, abundance, .abundance = NULL, of_samples = TRUE){
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_elements_features_abundance(.abundance)", "get_elements_features_abundance(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
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
#' @noRd
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
        tidybulk says: The function does not know what your sample, transcript and counts columns are.
	You might need to specify the arguments .sample, .transcript and/or .abundance. 
	Please read the documentation of this function for more information.
      ")
  }
}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param abundance A character name of the abundance column
#' @param .abundance DEPRECATED
#'
#' @return A list of column enquo or error
get_abundance_norm_if_exists = function(.data, abundance, .abundance = NULL){
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_abundance_norm_if_exists(.abundance)", "get_abundance_norm_if_exists(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }
  if (.abundance %>% quo_is_symbol()) {
    return(list(
      .abundance = .abundance
    ))
  } else {
    if (!is.null(get_tt_columns(.data))) {
      if (".abundance_scaled" %in% names(get_tt_columns(.data)) &&
          quo_name(get_tt_columns(.data)$.abundance_scaled) %in% colnames(.data)) {
        return(list(
          .abundance = get_tt_columns(.data)$.abundance_scaled
        ))
      } else {
        return(list(
          .abundance = get_tt_columns(.data)$.abundance
        ))
      }
    } else {
      stop("
        tidybulk says: The function does not know what your sample, transcript and counts columns are.
	You might need to specify the arguments .sample, .transcript and/or .abundance. 
	Please read the documentation of this function for more information.
      ")
    }
  }
}

#' Sub function of remove_redundancy_elements_though_reduced_dimensions
#'
#' @keywords internal
#' @noRd
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
#' @noRd
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



#' get_x_y_annotation_columns
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .horizontal The name of the column horizontally presented in the heatmap
#' @param .vertical The name of the column vertically presented in the heatmap
#' @param abundance A character name of the transcript/gene abundance column
#' @param .abundance DEPRECATED
#' @param .abundance_scaled The name of the transcript/gene scaled abundance column
#'
#' @description This function recognise what are the sample-wise columns and transcrip-wise columns
#'
#' @return A list
#'
get_x_y_annotation_columns = function(.data, .horizontal, .vertical, abundance, .abundance = NULL, .abundance_scaled = NULL){
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_x_y_annotation_columns(.abundance)", "get_x_y_annotation_columns(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }


  # Comply with CRAN NOTES
  . = NULL

  # Make col names
  .horizontal = enquo(.horizontal)
  .vertical = enquo(.vertical)
  .abundance = enquo(.abundance)
  .abundance_scaled = enquo(.abundance_scaled)

  # x-annotation df
  n_x = .data %>% select(!!.horizontal) |> distinct() |> nrow()
  n_y = .data %>% select(!!.vertical) |> distinct() |> nrow()

  # Sample wise columns
  horizontal_cols=
    .data %>%
    select(-!!.horizontal, -!!.vertical, -!!.abundance) %>%
    colnames %>%
    map(
      ~
        .x %>%
        when(
          .data %>%
            select(!!.horizontal, !!as.symbol(.x)) %>%
            distinct() |>
            nrow() %>%
            equals(n_x) ~ .x,
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
      ~{
          if(.data %>%
            select(!!.vertical, !!as.symbol(.x)) |>
            distinct() |>
            nrow() %>%
            equals(n_y)) {
          .x
        } else {
          NULL
        }
      }
       
    ) %>%

    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist

  # Counts wise columns, at the moment scaled counts is treated as special and not accounted for here
  counts_cols =
    .data %>%
    select(-!!.horizontal, -!!.vertical, -!!.abundance)

  if (!is.null(horizontal_cols)) {
    counts_cols <- counts_cols %>% select(-horizontal_cols)
  }
  if (!is.null(vertical_cols)) {
    counts_cols <- counts_cols %>% select(-vertical_cols)
  }
  if (.abundance_scaled %>% quo_is_symbol) {
    counts_cols <- counts_cols %>% select(-!!.abundance_scaled)
  }
  counts_cols <- counts_cols %>% colnames
  counts_cols <- counts_cols %>% map(
    ~ {
      n_unique <- .data %>% select(!!.vertical, !!.horizontal, !!as.symbol(.x)) %>% distinct() %>% nrow()
      if (n_unique == n_x * n_y) {
        .x
      } else {
        NULL
      }
    }
  )
  counts_cols <- counts_cols[lengths(counts_cols) != 0]
  counts_cols <- unlist(counts_cols)

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
    ~ {
      n_unique <- .data %>% distinct(!!.col, !!as.symbol(.x)) %>% nrow()
      if (n_unique == n_x) {
        .x
      } else {
        NULL
      }
    }
  ) %>%

  # Drop NULL
  { (.)[lengths((.)) != 0] } %>%
  unlist

}




#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#' @noRd
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

# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }


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

#'
#' @keywords internal
#' @noRd
#'
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

univariable_differential_tissue_stratification = function(
	deconvoluted,
	method,
	.my_formula
){

	# Check if package is installed, otherwise install
  check_and_install_packages(c("survival", "survminer"))



  check_and_install_packages("broom")
  

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
					 nrow() %>%
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
						conf.int = TRUE,
						palette = c("#ed6f68",  "#5366A0" ),
						legend = "none",
						pval = TRUE
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
  check_and_install_packages(c("survival", "survminer", "broom"))
  

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
					 nrow() %>%
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
						conf.int = TRUE,
						palette = c("#ed6f68",  "#5366A0" ),
						legend = "none",
						pval = TRUE
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
	((bind_rows(
		c(`1` = cos(r), `2` = -sin(r)),
		c(`1` = sin(r), `2` = cos(r))
	) %>% as_matrix) %*% m)
}

combineByRow <- function(m, fun = NULL) {
  # Shown here
  #https://stackoverflow.com/questions/8139301/aggregate-rows-in-a-large-matrix-by-rowname

  m <- m[ order(rownames(m)), ,drop=FALSE]

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
        m[i.start,] <- apply(m[i.start:i.end,,drop=FALSE], 2, fun)
        m.rownames[(1+i.start):i.end] <- NA
      }

      prev <- curr
      i.start <- i
    } else {
      i.end <- i
    }
  }

  m[ which(!is.na(m.rownames)),,drop=FALSE]
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

select_non_standard_column_class = function(.x){
  !is.numeric(.x) & !is.character(.x) & !is.factor(.x) & !is.logical(.x)
}

get_special_column_name_symbol = function(name){
  list(name = name, symbol = as.symbol(name))
}

feature__ =  get_special_column_name_symbol(".feature")
sample__ = get_special_column_name_symbol(".sample")

check_and_install_packages <- function(packages) {
  # Separate GitHub packages from CRAN/Bioconductor packages
  github_packages <- packages[grepl("/", packages)]
  regular_packages <- packages[!grepl("/", packages)]
  
  # Check if regular packages are installed
  missing_regular_packages <- regular_packages[!sapply(regular_packages, requireNamespace, quietly = TRUE)]
  
  # Check if GitHub packages are installed
  missing_github_packages <- github_packages[!sapply(gsub(".*/", "", github_packages), requireNamespace, quietly = TRUE)]
  
  # Combine all missing packages
  missing_packages <- c(missing_regular_packages, missing_github_packages)
  
  # If any packages are missing, print installation instructions
  if (length(missing_packages) > 0) {
    stop(
      "tidybulk says: The following packages are required:\n",
      paste("  -", missing_packages, collapse = "\n"), "\n",
      "Please install them by running:\n",
      "  if (!requireNamespace('BiocManager', quietly = TRUE))\n",
      "    install.packages('BiocManager', repos = 'https://cloud.r-project.org')\n",
      paste0(
        "  BiocManager::install(c(", 
        paste0("'", missing_packages, "'", collapse = ", "), 
        "), ask = FALSE)"
      )
    )
  }
}

#' Drop Environment from a Quosure
#'
#' Takes a quosure and resets its environment to `emptyenv()` without altering
#' its expression.
#'
#' @param q A quosure object to have its environment stripped.
#' @return A quosure with the same expression but environment set to `emptyenv()`.
#'
#' @importFrom rlang is_quosure
#' @importFrom rlang quo_set_env
#'
#' @examples
#' library(rlang)
#'
#' q <- quo(x + y)
#' environment(q)
#'
#' q_stripped <- drop_enquo_env(q)
#' identical(quo_get_env(q_stripped), emptyenv()) # TRUE
#'
#' @noRd
drop_enquo_env <- function(q) {
  if (!rlang::is_quosure(q)) {
    stop("`q` must be a quosure.")
  }
  rlang::quo_set_env(q, emptyenv())
}

#' Perform linear equation system analysis through linear least squares regression
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats lsfit
#'
#' @param mix A data frame containing mixture expression profiles
#' @param reference A data frame containing reference expression profiles for cell types (default = X_cibersort)
#' @param intercept Logical indicating whether to include intercept in the regression model (default = TRUE)
#'
#' @details
#' Performs cell type deconvolution using linear least squares regression. The function:
#' 1. Identifies common markers between mixture and reference
#' 2. Normalizes expression data
#' 3. Fits linear model with or without intercept
#' 4. Constrains results to non-negative values and normalizes to sum to 1
#'
#' @return A data frame containing estimated cell type proportions
#'
#'
run_llsr = function(mix, reference = X_cibersort,  intercept= TRUE) {
  # Get common markers
  markers = intersect(rownames(mix), rownames(reference))
  
  X <- (reference[markers, , drop = FALSE])
  Y <- (mix[markers, , drop = FALSE])
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  X <- (X - mean(X)) / sd(X)
  Y <- apply(Y, 2, function(mc) (mc - mean(mc)) / sd(mc)  )
  # Y <- (Y - mean(y)) / sd(Y)
  
  if(intercept)
    results <- t(data.frame(lsfit(X, Y)$coefficients)[-1, , drop = FALSE])
  else
    results <- t(data.frame(lsfit(X, Y, intercept=FALSE)$coefficients))
  results[results < 0] <- 0
  results <- results / apply(results, 1, sum)
  rownames(results) = colnames(Y)
  
  results
}

#' Perform linear equation system analysis through llsr
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats lsfit
#'
#' @param mix A data frame
#' @param reference A data frame
#'
#' @return A data frame
#'
#'
run_epic = function(mix, reference = NULL) {
  
  # Check if 'EPIC' package is installed, otherwise stop with instructions
  check_and_install_packages("EPIC")
  
  if("EPIC" %in% .packages() %>% not) stop("tidybulk says: Please install and then load the package EPIC manually (i.e. library(EPIC)). This is because EPIC is not in Bioconductor or CRAN so it is not possible to seamlessly make EPIC part of the dependencies.")
  
  # Get common markers
  if( reference  |> is("data.frame") | reference  |> is("matrix")){
    markers = intersect(rownames(mix), rownames(reference))
    
    X <- (reference[markers, , drop = FALSE])
    Y <- (mix[markers, , drop = FALSE])
    
    if(!is.null(reference))
      reference = list(
        refProfiles = X,
        sigGenes = rownames(X)
      )
  } else { Y <- mix }
  
  # Check if it is not matrix or data.frame, for example DelayedMatrix
  if(!is(Y, "matrix") & !is(Y, "data.frame"))
    Y = as.matrix(Y)
  
  results <- EPIC(Y, reference = reference)$cellFractions %>% data.frame()
  #results[results < 0] <- 0
  #results <- results / apply(results, 1, sum)
  rownames(results) = colnames(Y)
  
  results
}



counts_scaled_exist_SE = function(.data){
  
  ("tt_columns" %in% (.data %>%
                        attr("internal") %>% names())) &&
    (
      .data %>%
        attr("internal") %$%
        tt_columns %>%
        names() %>%
        grep("scaled", .) %>%
        length() %>%
        equals(1)
    )
}

get_assay_scaled_if_exists_SE = function(.data){
  if(counts_scaled_exist_SE(.data))
    .data %>%
    attr("internal") %$%
    tt_columns %$%
    .abundance_scaled %>%
    quo_name()
  else
    .data %>%
    assays() %>%
    names() %>%
    head(1)
}

filter_if_abundant_were_identified = function(.data){
  # Filter abundant if performed
  if (".abundant" %in% (rowData(.data) %>% colnames())) {
    .data[rowData(.data)[,".abundant"],]
  } else {
    warning("tidybulk says: highly abundant transcripts were not identified (i.e. identify_abundant()) or filtered (i.e., keep_abundant), therefore this operation will be performed on unfiltered data. In rare occasions this could be wanted. In standard whole-transcriptome workflows is generally unwanted.")
    .data
  }
}



#' Resolve complete confounders of non-interest in a data frame
#'
#' This internal function identifies and resolves complete confounders among specified columns (factors of non-interest) in a data frame.
#' For each specified column, a new column with the suffix "___altered" is created to preserve the original values.
#' The function then examines all pairs of these altered columns to detect and resolve complete confounding relationships.
#' Additionally, it checks for columns with only one unique value, which cannot be estimated in a linear model, and issues a warning if any are found.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stringr str_remove
#' @importFrom stringr str_subset
#' @importFrom dplyr across
#' @importFrom dplyr n_distinct
#' @importFrom tibble as_tibble
#'
#' @param df A data frame containing the data to be checked for confounders.
#' @param ... Unquoted column names to consider as confounders (factors of non-interest).
#'
#' @return A data frame with new columns (with "___altered" suffix) where confounders have been resolved.
#'
.resolve_complete_confounders_of_non_interest_df <- function(df, ...) {
  
  # Step 1: Create new columns with ___altered suffix for each specified confounder
  df <- df |>
    mutate(across(c(...), ~ .x, .names = "{col}___altered"))
  
  # Step 2: Identify all pairs of altered columns to check for complete confounding
  combination_of_factors_of_NON_interest <-
    df |>
    select(ends_with("___altered")) |>
    colnames() |>
    combn(2) |>
    t() |>
    as_tibble() |>
    set_names(c("factor_1", "factor_2"))
  
  # Step 3: Inform the user about the new columns created
  message(
    "tidybulk says: New columns created with resolved confounders: ",
    paste0(colnames(df) |> str_subset("___altered"), collapse = ", ")
  )
  
  # Step 4: For each pair of altered columns, resolve complete confounding if present
  for (i in seq_len(nrow(combination_of_factors_of_NON_interest))) {
    df <-
      df |>
      resolve_complete_confounders_of_non_interest_pair_df(
        !!as.symbol(combination_of_factors_of_NON_interest[i, ]$factor_1),
        !!as.symbol(combination_of_factors_of_NON_interest[i, ]$factor_2)
      )
  }
  
  # Step 5: Check for columns with only one unique value (cannot be estimated in a linear model)
  single_value_cols <- df |>
    select(ends_with("___altered")) |>
    summarise(across(everything(), ~ n_distinct(.x))) |>
    pivot_longer(everything()) |>
    filter(value == 1) |>
    pull(name)
  
  if (length(single_value_cols) > 0) {
    warning(
      "tidybulk says: The following columns have only one unique value and cannot be estimated by a linear model: ",
      paste(single_value_cols, collapse = ", ")
    )
  }
  
  # Step 6: Return the modified data frame with resolved confounders
  df
}

#' Resolve Complete Confounders of Non-Interest
#'
#' This function processes a SummarizedExperiment object to handle confounders
#' that are not of interest in the analysis. It deals with two factors, adjusting
#' the data by nesting and summarizing over these factors.
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_name
#' @importFrom rlang !!
#' @importFrom tidyr unnest
#' @importFrom tidyr nest
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom purrr map_int
#' @importFrom dplyr if_else
#'
#' @param se A SummarizedExperiment object that contains the data to be processed.
#' @param .factor_1 A symbol or quosure representing the first factor variable in `se`.
#' @param .factor_2 A symbol or quosure representing the second factor variable in `se`.
#' @return A modified SummarizedExperiment object with confounders resolved.
#' @examples
#' # Not run:
#' # se is a SummarizedExperiment object
#' resolve_complete_confounders_of_non_interest(se, .factor_1 = factor1, .factor_2 = factor2)
#' @noRd
resolve_complete_confounders_of_non_interest_pair_df <- function(df, .factor_1, .factor_2){
  
  # Fix NOTEs
  . = NULL
  
  .factor_1 = enquo(.factor_1)
  .factor_2 = enquo(.factor_2)
  
  cd =
    df |>
    as_tibble() |>
    rowid_to_column() |>
    distinct(rowid, !!.factor_1, !!.factor_2) |>
    
    nest(se_data = -c(!!.factor_1, !!.factor_2)) |>
    nest(data = -!!.factor_1) |>
    mutate(n1 = map_int(data, ~ .x |> distinct(!!.factor_2) |> nrow())) |>
    unnest(data) |>
    
    # Nest data excluding .factor_2 and count distinct .factor_1 values
    nest(data = -!!.factor_2) |>
    mutate(n2 = map_int(data, ~ .x |> distinct(!!.factor_1) |> nrow())) |>
    unnest(data)
  
  # Choose a dummy value for .factor_2 based on sorting by n1 + n2
  dummy_factor_2 <- cd |> arrange(desc(n1 + n2)) |> slice(1) |> pull(!!.factor_2)
  
  # Messages if I have confounders
  if(cd |> filter(n1 + n2 < 3) |> nrow() > 0){
    
    message(sprintf("tidybulk says: IMPORTANT! the columns %s and %s, have been corrected for complete confounders and now are NOT interpretable. \n      They cannot be used in hypothesis testing. However they can be used in the model to capture the unwanted variability in the data.", quo_name(.factor_1), quo_name(.factor_2)))
    
    message(sprintf(
      "tidybulk says: The value(s) %s in column %s from sample(s) %s, has been changed to %s.",
      cd |> filter(n1 + n2 < 3) |> pull(!!.factor_2),
      quo_name(.factor_2),
      cd |> filter(n1 + n2 < 3) |> pull(se_data) |> map(colnames) |> unlist(),
      dummy_factor_2
    ))
    
    # Replace .factor_2 with dummy_factor_2 where n1 + n2 is less than 3 and unnest
    cd = cd |>
      mutate(!!.factor_2 := if_else(n1 + n2 < 3, dummy_factor_2, !!.factor_2))
  }
  
  df[,c(quo_name(.factor_1), quo_name(.factor_2))] =
    cd |>
    unnest(se_data) |>
    arrange(rowid) |>
    select(!!.factor_1, !!.factor_2)
  
  df
}

# Global variable declarations to avoid R CMD check warnings
globalVariables(c("transcript", "read count", "n", "m", ".", ".feature", ".abundance_scaled", 
                  "tt_columns", "seurat_clusters", "cluster", "tagwise.dispersion", "Component", 
                  "Component value", "sdev", "name", "value", "x", "Y", "x"))
