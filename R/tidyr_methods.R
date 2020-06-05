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
	tidyr::nest(.data, ...)
}

#' @export
#' @rdname tidyr-methods
nest.tidybulk <- function (.data, ...)
{
	warning("nest is not fully supported yet by tidybulk. The nested data frame has been reverted to tbl_df")

	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		drop_internals() %>%
		tidyr::nest(...)

	#   %>%
	#
	# 	# Attach attributes
	# 	reattach_internals(.data) %>%
	#
	# 	# Add class
	# 	add_class("tt") %>%
	# 	add_class("tidybulk")

}
