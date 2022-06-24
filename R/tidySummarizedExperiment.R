change_reserved_column_names = function(.data){

	.data %>%

		setNames(
			colnames(.) %>%
				str_replace("^feature$", "feature.x") %>%
				str_replace("^sample$", "sample.x") %>%
				str_replace("^coordinate$", "coordinate.x")
		)

}

#' @importFrom dplyr select
#' @importFrom tidyselect one_of
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom tibble rowid_to_column
#'
#' @keywords internal
#' @noRd
get_special_datasets <- function(SummarizedExperiment_object) {

	SummarizedExperiment_object %>%
		rowRanges() %>%
		when(
			# If no ranges
			as.data.frame(.) %>%
				nrow() %>%
				equals(0) ~ tibble(),

			# If it is a range list (multiple rows per feature)
			class(.) %>% equals("CompressedGRangesList") ~
				tibble::as_tibble(.) %>%
				eliminate_GRanges_metadata_columns_also_present_in_Rowdata(SummarizedExperiment_object) %>%
				nest(coordinate = -group_name) %>%
				rename(feature = group_name),

			# If standard GRanges (one feature per line)
			~ {
				transcript_column =
					rowRanges(SummarizedExperiment_object) %>%
					as.data.frame() %>%
					lapply(function(x) rownames(SummarizedExperiment_object)[1] %in% x) %>%
					unlist() %>%
					which() %>%
					names()


				# Just rename
				(.) %>%

					# If transcript_column exists all good
					when(
						!is.null(transcript_column) ~  tibble::as_tibble(.) %>%
							eliminate_GRanges_metadata_columns_also_present_in_Rowdata(SummarizedExperiment_object) %>%
							rename(feature := !!transcript_column) ,

						# If transcript_column is NULL add numeric column
						~ tibble::as_tibble(.) %>%
							eliminate_GRanges_metadata_columns_also_present_in_Rowdata(SummarizedExperiment_object) %>%
							rowid_to_column(var = "feature") %>%
							mutate(feature = as.character(feature))
					) %>%

					# Always nest
					nest(coordinate = -feature)

			}
		) %>%
		list()

}

change_reserved_column_names = function(.data){

	.data %>%

		setNames(
			colnames(.) %>%
				str_replace("^feature$", "feature.x") %>%
				str_replace("^sample$", "sample.x") %>%
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
#' @keywords internal
#' @noRd
get_count_datasets <- function(SummarizedExperiment_object, feature_column_name, sample_column_name) {


	map2(
		assays(SummarizedExperiment_object) %>% as.list(),
		names(assays(SummarizedExperiment_object)),
		~ .x %>%
			tibble::as_tibble(rownames = feature_column_name, .name_repair = "minimal") %>%

			# If the matrix does not have sample names, fix column names
			when(colnames(.x) %>% is.null() ~ setNames(., c(
			  feature_column_name,  seq_len(ncol(.x))
			)),
			~ (.)
			) %>%

			gather(!!as.symbol(sample_column_name), count,-!!as.symbol(feature_column_name)) %>%
			rename(!!.y := count)
	) %>%
		reduce(left_join, by = c(feature_column_name, sample_column_name))
}
