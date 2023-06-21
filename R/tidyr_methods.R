#' unnest
#'
#' @importFrom tidyr unnest
#'
#' @param data A tbl. (See tidyr)
#' @param cols <[`tidy-select`][tidyr_tidy_select]> Columns to unnest.
#'   If you `unnest()` multiple columns, parallel entries must be of
#'   compatibble sizes, i.e. they're either equal or length 1 (following the
#'   standard tidyverse recycling rules).
#' @param ... <[`tidy-select`][tidyr_tidy_select]> Columns to nest, specified
#'   using name-variable pairs of the form `new_col=c(col1, col2, col3)`.
#'   The right hand side can be any valid tidy select expression.
#'
#'   \Sexpr[results=rd, stage=render]{lifecycle::badge("deprecated")}:
#'   previously you could write `df %>% nest(x, y, z)` and `df %>%
#'   unnest(x, y, z)`. Convert to `df %>% nest(data=c(x, y, z))`.
#'   and `df %>% unnest(c(x, y, z))`.
#'
#'   If you previously created new variable in `unnest()` you'll now need to
#'   do it explicitly with `mutate()`. Convert `df %>% unnest(y=fun(x, y, z))`
#'   to `df %>% mutate(y=fun(x, y, z)) %>% unnest(y)`.
#' @param names_sep If `NULL`, the default, the names will be left
#'   as is. In `nest()`, inner names will come from the former outer names;
#'   in `unnest()`, the new outer names will come from the inner names.
#'
#'   If a string, the inner and outer names will be used together. In `nest()`,
#'   the names of the new outer columns will be formed by pasting together the
#'   outer and the inner column names, separated by `names_sep`. In `unnest()`,
#'   the new inner names will have the outer names (+ `names_sep`) automatically
#'   stripped. This makes `names_sep` roughly symmetric between nesting and unnesting.
#' @param keep_empty See tidyr::unnest
#' @param names_repair See tidyr::unnest
#' @param ptype See tidyr::unnest
#' @param .drop See tidyr::unnest
#' @param .id tidyr::unnest
#' @param .sep tidyr::unnest
#' @param .preserve See tidyr::unnest
#'
#'
#' @return A tidySummarizedExperiment objector a tibble depending on input
#'
#' @examples
#'
#'
#' tidybulk::se_mini |> tidybulk() |> nest( data = -.feature) |> unnest(data)
#'
#' @rdname nest-methods
#' @name unnest
#'
#' @export
NULL

#' @export
unnest.nested_tidybulk <- function (data, cols, ..., keep_empty=FALSE, ptype=NULL, names_sep=NULL, names_repair="check_unique", .drop, .id, .sep, .preserve)
{

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
#' tidybulk::se_mini %>% tidybulk() %>% nest( data = -.feature)
#'
#' @rdname nest-methods
#' @name nest
#'
#' @export
NULL

#' @importFrom rlang enquos
#' 
#' @export
#' 
#' 
#' 
nest.tidybulk <- function (.data, ..., .names_sep = NULL)
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
