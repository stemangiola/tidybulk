#' unnest
#'
#' @importFrom tidyr unnest
#'
#' @param .data A tbl. (See tidyr)
#' @param ... Name-variable pairs of the form new_col = c(col1, col2, col3) (See tidyr)
#'
#' @return A tt object
#'
#' @examples
#'
#' unnest(tidybulk(tidybulk::counts_mini, sample, transcript, count), data = -transcript)
#'
#' @rdname tidyr-methods
#'
#' @export
unnest <- function (.data, cols, ..., keep_empty = FALSE, ptype = NULL, 
										names_sep = NULL, names_repair = "check_unique")  {
	UseMethod("unnest")
}

#' @export
#' @rdname tidyr-methods
unnest.default <-  function (.data, cols, ..., keep_empty = FALSE, ptype = NULL, 
														 names_sep = NULL, names_repair = "check_unique")
{
	tidyr::unnest(.data, cols, ..., keep_empty = keep_empty, ptype = ptype, 
								names_sep = names_sep, names_repair = names_repair)
}

#' @export
#' @rdname tidyr-methods
unnest.nested_tidybulk <- function (.data, cols, ..., keep_empty = FALSE, ptype = NULL, 
														 names_sep = NULL, names_repair = "check_unique")
{

	cols <- enquo(cols)
	
	
	.data %>%
		drop_class(c("nested_tidybulk", "tt")) %>%
		tidyr::unnest(!!cols, ..., keep_empty = keep_empty, ptype = ptype, 
								names_sep = names_sep, names_repair = names_repair) %>%

		# Attach attributes
		reattach_internals(.data) %>%
		
		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}

#' nest
#'
#' @importFrom tidyr nest
#'
#' @param .data A tbl. (See tidyr)
#' @param ... Name-variable pairs of the form new_col = c(col1, col2, col3) (See tidyr)
#'
#' @return A tt object
#'
#' @examples
#'
#' nest(tidybulk(tidybulk::counts_mini, sample, transcript, count), data = -transcript)
#'
#' @rdname tidyr-methods
#'
#' @export
nest <- function (.data, ...)  {
	UseMethod("nest")
}

#' @export
#' @rdname tidyr-methods
nest.default <-  function (.data, ...)
{
	tidyr::nest(.data, cols, ..., keep_empty = keep_empty, ptype = ptype, 
								names_sep = names_sep, names_repair = names_repair)
}

#' @export
#' @rdname tidyr-methods
nest.tidybulk <- function (.data, ...)
{
	cols <- enquos(...)
	col_name_data  = names(cols)
	
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