#' @name unnest
#' @rdname unnest
#' @inherit tidyr::unnest
#' @return A tidySummarizedExperiment objector a tibble depending on input
#' 
#' @examples
#' tidybulk::se_mini |> tidybulk() |> nest( data = -.feature) |> unnest(data)
#'
#' @importFrom tidyr unnest
#' @importFrom rlang enquo
#' @export
unnest.nested_tidybulk <- function (data, cols, ..., keep_empty=FALSE, ptype=NULL, 
                                    names_sep=NULL, names_repair="check_unique", 
                                    .drop, .id, .sep, .preserve) {
	cols <- enquo(cols)
	data %>%
		drop_class(c("nested_tidybulk", "tt")) %>%
		tidyr::unnest(!!cols, ..., keep_empty = keep_empty, ptype = ptype,
								names_sep = names_sep, names_repair = names_repair) %>%

		# Attach attributes
		reattach_internals(data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")
}

#' @name nest
#' @rdname nest
#' @inherit tidyr::nest
#' @return A tt object
#' 
#' @examples
#' tidybulk::se_mini %>% tidybulk() %>% nest( data = -.feature)
#' 
#' @importFrom rlang enquos
#' @importFrom tidyr nest
#' @export
nest.tidybulk <- function (.data, ..., .names_sep = NULL) {
	cols <- enquos(...)
	col_name_data <- names(cols)
	.data %>%

		# This is needed otherwise nest goes into loop and fails
		drop_class(c("tidybulk", "tt")) %>%
		tidyr::nest(...) %>%

		# Add classes afterwards
		mutate(!!as.symbol(col_name_data) := map(
			!!as.symbol(col_name_data),
			~ .x %>%
				add_class("tt") %>%
				add_class("tidybulk")
		)) %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("nested_tidybulk")
}
