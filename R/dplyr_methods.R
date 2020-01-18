#' Join datasets
#'
#' @importFrom purrr map
#' @importFrom stats setNames
#'
#' @param ... Data frames to combine
#'
#' @return A tt object
#'
#' @examples
#'
#' bind_rows(ttBulk::counts_mini, ttBulk::counts_mini)
#'
#'
#' @export
bind_rows <- function(...) {
	UseMethod("bind_rows")
}

#' @export
bind_rows.default <-  function(...)
{
	dplyr::bind_rows(...)
}

#' @export
bind_rows.ttBulk <- function(...)
{

	tts = dplyr:::flatten_bindable(rlang::dots_values(...))

	par1 = tts[[1]] %>% attr("parameters") %>% unlist
	par2 = tts[[2]] %>% attr("parameters") %>% unlist

	# Parameters of the two objects must match
	error_if_parameters_not_match(par1, par2)

	par =
		unique(c(par1 %>% names, par2 %>% names)) %>%
		map(~ switch(par1[[.x]] %>% is.null %>% sum(1), par1[[.x]], par2[[.x]])) %>%
		setNames(par1 %>% names)


	dplyr::bind_rows(...) %>%

		# Attach attributes
		add_attr(par, "parameters")

}


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

#' distinct
#' @param .data A tbl. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#' @param .keep_all If TRUE, keep all variables in .data. If a combination of ... is not distinct, this keeps the first row of values. (See dplyr)
#'
#' @return A tt object
#'
#' @examples
#'
#' distinct(ttBulk::counts_mini)
#'
#'
#' @export
distinct <- function (.data, ..., .keep_all = FALSE)  {
	UseMethod("distinct")
}

#' @export
distinct.default <-  function (.data, ..., .keep_all = FALSE)
{
	dplyr::distinct(.data, ..., .keep_all = FALSE)
}

#' @export
distinct.ttBulk <- function (.data, ..., .keep_all = FALSE)
{
	.data %>%
		drop_class(c("ttBulk", "tt")) %>%
		dplyr::distinct(..., .keep_all = .keep_all) %>%

		# Attach attributes
		add_attr(.data %>% attr("parameters"), "parameters") %>%

		# Add class
		add_class("tt") %>%
		add_class("ttBulk")

}
