

#' Check whether there are NA counts
#'
#' @keywords internal
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
check_if_wrong_input <- function(.data, list_input, expected_type) {
	# Do the check
	if (list_input %>%
			map( ~ .x %>% class() %>% `[` (1)) %>%
			unlist %>%
			equals(expected_type) %>%
			not())
		stop("tidybulk says: You have passed the wrong argument to the function. Please check again.")

	# If all good return original data frame
	.data
}

#' Check whether there are duplicated genes/transcripts
#'
#' @keywords internal
#'
#' 
#' 
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
check_if_duplicated_genes <- function(.data,
																			.sample = `sample`,
																			.transcript = `transcript`,
																			.abundance = `read count`) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	paste(pull(.data, !!.sample), pull(.data, !!.transcript) ) %>% duplicated() %>% which() %>% length() %>% equals(0)

}

#' Check whether there are NA counts
#'
#' @keywords internal
#'
#' 
#' 
#' @import tibble
#'
#' @param .data A tibble of read counts
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
#'
check_if_counts_is_na = function(.data, .abundance) {
	.abundance = enquo(.abundance)

	.data %>% filter(!!.abundance %>% is.na) %>% nrow() %>% equals(0)

}

check_if_transcript_is_na = function(.data, .transcript) {
	.transcript = enquo(.transcript)

	.data %>% filter(!!.transcript %>% is.na) %>% nrow() %>% equals(0)

}

check_if_column_missing = function(.data, .sample, .transcript, .abundance) {
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Check that the intersection is length 3
	.data %>% colnames %>%
		intersect(c(
			quo_name(.sample),
			quo_name(.transcript),
			quo_name(.abundance)
		)) %>%
		length %>%
		equals(3)
}

#' @importFrom dplyr pull
column_type_checking = function(.data, .sample, .transcript, .abundance) {
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>% pull(!!.sample) %>% class %in% c("character", "factor") &
		.data %>% pull(!!.transcript) %>% class %in% c("character", "factor") &
		.data %>% pull(!!.abundance) %>% class %in% c("integer", "numeric", "double")

}

check_if_attribute_present = function(.data) {
	"internals" %in% (.data %>% attributes %>% names) &&
	"tt_columns" %in% (.data %>% attr("internals")  %>% names)
}

eliminate_sparse_transcripts <- function(.data, .transcript){
	# Parse column names
	.transcript = enquo(.transcript)
  
	warning("tidybulk says: Some transcripts have been omitted from the analysis ",
	        "because not present in every sample.")

	.data %>%
		add_count(!!.transcript, name = "my_n") %>%
		filter(my_n == max(my_n)) %>%
		select(-my_n)
}

check_if_data_rectangular = function(.data, .sample, .transcript, .abundance){
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	is_rectangular =
		.data %>%
		distinct(!!.sample, !!.transcript, !!.abundance) %>%
		count(!!.sample) %>%
		count(n, name = "nn") %>%
		nrow() %>%
		equals(1)

	is_rectangular

}

warning_if_data_is_not_rectangular <- function(.data, .sample, .transcript, .abundance) {
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	if(!check_if_data_rectangular(.data, !!.sample, !!.transcript, !!.abundance))
		warning("tidybulk says: the data does not have the same number of transcript ",
		"per sample. The data set is not rectangular.")

}

error_if_data_is_not_rectangular = function(.data, .sample, .transcript, .abundance){

	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	if(!check_if_data_rectangular(.data, !!.sample, !!.transcript, !!.abundance))
		stop("tidybulk says: the data must have the same number of transcript per sample. ",
		"Check again that you have not filtered single observations accidentally. ",
		"If you have missing data you can use fill_missing_abundance() or impute_missing_abundance()")
}

tidybulk_to_tbl = function(.data) {
	.data %>%	drop_class(c("tidybulk", "tt"))
}


validation_default <- function(.data,
															.sample,
															.transcript,
															.abundance,
															type = "hard",
															skip_dupli_check = FALSE) {
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Type check
	is_missing = check_if_column_missing(.data,!!.sample,!!.transcript,!!.abundance)
	if (type == "hard" &
			!is_missing)
		stop(
			"tidybulk says: One or more columns that should include sample identifier, ",
			"transcript identified or transcript abundance are missing from your data frame."
		)
	if (type == "soft" & !is_missing) {
		warning(
			"tidybulk says: One or more columns that should include sample identifier, ",
			"transcript identified or transcript abundance are missing from your data frame. ",
			"The tidybulk object has been converted to a `tbl`"
		)
		return(.data %>% tidybulk_to_tbl)
	}

	# Type check
	is_type = column_type_checking(.data,!!.sample,!!.transcript,!!.abundance)
	if (type == "hard" &
			!is_type)
		stop(
			"tidybulk says: The column provided as .sample .transcript or .abundance do not ",
			"comply with the required types (<FACTOR/CHARACTER>, <FACTOR/CHARACTER>, <NUMERIC>)."
		)
	if (type == "soft" & !is_type) {
		warning(
			"tidybulk says: The column provided as .sample .transcript or .abundance do not ",
			"comply with the required types. The tidybulk object has been converted to a `tbl`"
		)
		return(.data %>% tidybulk_to_tbl)
	}

	# Check if duplicated genes
	if (!skip_dupli_check) {
		is_unique = check_if_duplicated_genes(.data,!!.sample,!!.transcript,!!.abundance)
		if (type == "hard" & !is_unique){

			dup = paste(pull(.data, !!.sample), pull(.data, !!.transcript) )

			stop(
				"tidybulk says: Your dataset include duplicated sample/gene pairs. ",
				dup[dup %>% duplicated()] %>% head(30) %>% paste(collapse=", "),
				"Please, remove redundancies before proceeding (e.g., aggregate_duplicates())."
			)
		}
		if (type == "soft" & !is_unique) {
			warning(
				"tidybulk says: Your dataset include duplicated sample/gene pairs. ",
				dup[dup %>% duplicated()] %>% paste(collapse=", "),
				" Please, remove redundancies before proceeding (e.g., aggregate_duplicates()). The tidybulk object has been converted to a `tbl`"
			)
			return(.data %>% tidybulk_to_tbl)
		}
	}

	# Check if NA in counts
	is_count_good = check_if_counts_is_na(.data,!!.abundance)
	if (type == "hard" &
			!is_count_good)
		stop("tidybulk says: You have NA values in your counts. Please check your data frame.")
	if (type == "soft" & !is_count_good) {
		warning(
			"tidybulk says: You have NA values in your counts. The tidybulk object has been converted to a `tbl`"
		)
		return(.data %>% tidybulk_to_tbl)
	}

}

validation <- function(.data,
											 .sample = NULL,
											 .transcript = NULL,
											 .abundance = NULL,
											 type = "hard",
											 skip_dupli_check = FALSE) {
	UseMethod("validation", .data)
}

validation.default <- validation_default

validation.tidybulk <- function(.data,
														 .sample = NULL,
														 .transcript = NULL,
														 .abundance = NULL,
														 type = "hard",
														 skip_dupli_check = FALSE) {
	# Check if attribute is present
	is_attr = check_if_attribute_present(.data)
	if (type == "hard" &
			!is_attr)
		stop(
			"tidybulk says: The object provided has tidybulk class but no attribute ",
			"containing the column names (attr(., \"internals\")). You must have used ",
			"an external function that eliminated the attributes. Insert a valid ",
			"tidybulk object or provide `.sample`, `.transcript`, `.abundance` column names as arguments "
		)
	if (type == "soft" & !is_attr) {
		warning(
			"tidybulk says: The object provided has tidybulk class but no attribute ",
			"containing the column names (attr(., \"internals\")). You must have used an ",
			"external function that eliminated the attributes. The tidybulk object has been converted to a `tbl`"
		)
		return(.data %>% tidybulk_to_tbl)
	}

	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
	.sample = col_names$.sample
	.transcript = col_names$.transcript
	.abundance = col_names$.abundance

	validation_default(
		.data,
		!!.sample,
		!!.transcript,
		!!.abundance,
		type = type,
		skip_dupli_check = skip_dupli_check
	)

}

validate_signature = function(.data, reference, .transcript){

	.transcript = enquo(.transcript)

	overlapping_genes = .data %>%     pull(!!.transcript) %in% rownames(reference) %>%  which

	if(length(overlapping_genes) == 0  )
	  stop(sprintf(
	    "\ntidybulk says: You have NO genes in common between the query data and the reference data. Please check again your input dataframes\nthe genes in the reference look like this %s", 
	    paste(rownames(reference)[1:10], collapse = ", ")
	  ))

	if ( length(overlapping_genes) %>%	st(50) )
	  warning(sprintf(
	    "\ntidybulk says: You have less than 50 genes in common between the query data and the reference data. Please check again your input dataframes\nthe genes in the reference look like this %s", 
	    paste(rownames(reference)[1:10], collapse = ", ")
	  ))

	# Check if rownames exist
	if (reference %>% sapply(class) %in% c("numeric", "double", "integer") %>% not() %>% any)
	  stop("tidybulk says: your reference has non-numeric/integer columns.")


}

validate_signature_SE <- function(assay, reference) {
  overlapping_genes = (rownames(assay)  %in% rownames(reference)) %>%  which

  if(length(overlapping_genes) == 0  )
    stop(sprintf(
      "\ntidybulk says: You have NO genes in common between the query data and ",
      "the reference data. Please check again your input dataframes\nthe genes in the reference look like this %s", paste(rownames(reference)[1:10], collapse = ", ")
    ))

	if ( length(overlapping_genes) %>%	st(50) )
	  warning(sprintf(
			"\ntidybulk says: You have less than 50 genes in common between the query data and the reference data. Please check again your input dataframes\nthe genes in the reference look like this %s", 
			paste(rownames(reference)[1:10], collapse = ", ")
		))

	# Check if rownames exist
	if (reference %>% sapply(class) %in% c("numeric", "double", "integer") %>% not() %>% any)
		stop("tidybulk says: your reference has non-numeric/integer columns.")

}
