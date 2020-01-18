
#' Mutate datasets
#'
#' @param .data A tbl. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tt object
#' @export
mutate <- function(.data, ...) {
	UseMethod("mutate")
}

#' @export
mutate.default <-  function(.data, ...)
{
	dplyr::mutate(.data, ...)
}

#' @export
mutate.ttBulk <- function(.data, ...)
{
	.data %>%
		drop_class(c("ttBulk", "tt")) %>%
		dplyr::mutate(...) %>%

		# Attach attributes
		add_attr(.data %>% attr("parameters"), "parameters") %>%

		# Add class
		add_class("tt") %>%
		add_class("ttBulk")


}


#' Left join datasets
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tt object
#'
#' @export
left_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
	UseMethod("left_join")
}

#' @export
left_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
																...)
{
	dplyr::left_join(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...)
}

#' @export
left_join.ttBulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
															...)
{
	x %>%
		drop_class(c("ttBulk", "tt")) %>%
		dplyr::left_join(y, by = by, copy = copy, suffix = suffix, ...) %>%

		# Attach attributes
		add_attr(x %>% attr("parameters"), "parameters") %>%

		# Add class
		add_class("tt") %>%
		add_class("ttBulk")

}

#' Inner join datasets
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tt object
#'
#' @export
inner_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
	UseMethod("inner_join")
}

#' @export
inner_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
																 ...)
{
	dplyr::inner_join(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...)
}

#' @export
inner_join.ttBulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
															 ...)
{
	x %>%
		drop_class(c("ttBulk", "tt")) %>%
		dplyr::inner_join(y, by = by, copy = copy, suffix = suffix, ...) %>%

		# Attach attributes
		add_attr(x %>% attr("parameters"), "parameters") %>%

		# Add class
		add_class("tt") %>%
		add_class("ttBulk")

}


