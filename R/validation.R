

#' Check whether there are NA counts
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
check_if_wrong_input = function(.data, list_input, expected_type) {
	# Do the check
	if (list_input %>%
			map( ~ .x %>% class() %>% `[` (1)) %>%
			unlist %>%
			equals(expected_type) %>%
			`!`)
		stop("You have passed the wrong argument to the function. Please check again.")

	# If all good return original data frame
	.data
}

#' Check whether there are duplicated genes/transcripts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
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

	duplicates <-
		distinct(.data,!!.sample, !!.transcript, !!.abundance) %>%
		count(!!.sample, !!.transcript) %>%
		filter(n > 1) %>%
		arrange(n %>% desc())

	is_unique = duplicates %>% nrow() %>% equals(0)

	if (!is_unique) {
		writeLines("Those are the duplicated genes")
		duplicates %>% print()
	}

	is_unique
}

#' Check whether there are NA counts
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
check_if_counts_is_na = function(.data, .abundance) {
	.abundance = enquo(.abundance)

	.data %>% filter(!!.abundance %>% is.na) %>% nrow %>% equals(0)

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
	"parameters" %in% (.data %>% attributes %>% names)
}

eliminate_sparse_transcripts = function(.data, .transcript){
	# Parse column names
	.transcript = enquo(.transcript)

	warning("Some transcripts have been omitted from the analysis because not present in every sample.")

	.data %>%
		add_count(!!.transcript, name = "my_n") %>%
		filter(my_n == max(my_n)) %>%
		select(-my_n)
}

check_if_data_rectangular = function(.data, .sample, .transcript, .abundance, type = "hard"){

	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	is_rectangular =
		.data %>%
		distinct(!!.sample, !!.transcript, !!.abundance) %>%
		count(!!.sample) %>%
		count(n) %>%
		nrow %>%
		equals(1)

	is_rectangular

	# if(!is_rectangular & type == "hard") stop("ttBulk says: the data must have the same number of transcript per sample.")
	#
	# if(!is_rectangular & type == "soft") warning("ttBulk says: the data should have the same number of transcript per sample.")


	# # Eliminate sparse transcripts
	# .data %>% eliminate_sparse_transcripts(!!.transcript)


}

ttBulk_to_tbl = function(.data) {
	.data %>%	drop_class(c("ttBulk", "tt"))
}

validation_default = function(.data,
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
			"ttBulk says: One or more columns .sample .transcript or .abundance are missing from your data frame."
		)
	if (type == "soft" & !is_missing) {
		warning(
			"ttBulk says: One or more columns .sample .transcript or .abundance are missing from your data frame. The ttBulk object has been converted to a `tbl`"
		)
		return(.data %>% ttBulk_to_tbl)
	}

	# Type check
	is_type = column_type_checking(.data,!!.sample,!!.transcript,!!.abundance)
	if (type == "hard" &
			!is_type)
		stop(
			"ttBulk says: The column provided as .sample .transcript or .abundance do not comply with the required types."
		)
	if (type == "soft" & !is_type) {
		warning(
			"ttBulk says: The column provided as .sample .transcript or .abundance do not comply with the required types. The ttBulk object has been converted to a `tbl`"
		)
		return(.data %>% ttBulk_to_tbl)
	}

	# Check if duplicated genes
	if (!skip_dupli_check) {
		is_unique = check_if_duplicated_genes(.data,!!.sample,!!.transcript,!!.abundance)
		if (type == "hard" &
				!is_unique)
			stop(
				"ttBulk says: Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding (e.g., aggregate_duplicates())."
			)
		if (type == "soft" & !is_unique) {
			warning(
				"ttBulk says: Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding (e.g., aggregate_duplicates()). The ttBulk object has been converted to a `tbl`"
			)
			return(.data %>% ttBulk_to_tbl)
		}
	}

	# Check if NA in counts
	is_count_good = check_if_counts_is_na(.data,!!.abundance)
	if (type == "hard" &
			!is_count_good)
		stop("ttBulk says: You have NA values in your counts. Please check your data frame.")
	if (type == "soft" & !is_count_good) {
		warning(
			"ttBulk says: You have NA values in your counts. The ttBulk object has been converted to a `tbl`"
		)
		return(.data %>% ttBulk_to_tbl)
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

validation.default = validation_default

validation.ttBulk = function(.data,
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
			"ttBulk says: The object provided has ttBulk class but no attribute containing the column names. Insert a valid ttBulk object or provide `.sample`, `.transcript`, `.abundance` column names as arguments "
		)
	if (type == "soft" & !is_attr) {
		warning(
			"ttBulk says: The object provided has ttBulk class but no attribute containing the column names. The ttBulk object has been converted to a `tbl`"
		)
		return(.data %>% ttBulk_to_tbl)
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