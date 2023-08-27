#' @export
dplyr::select

#' @name arrange
#' @rdname arrange
#' @inherit dplyr::arrange
#' @family single table verbs
#' @importFrom dplyr arrange
#' @export
arrange.tidybulk <- function(.data, ..., .by_group = FALSE) {
	.data |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::arrange( ..., .by_group = .by_group) |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name bind_rows
#' @rdname bind_rows
#' @inherit ttservice::bind_rows
#'
#' @examples
#' data(se_mini)
#' se_mini_tidybulk <- se_mini |> tidybulk()
#' bind_rows(se_mini_tidybulk, se_mini_tidybulk)
#'
#' tt_bind <- se_mini_tidybulk |> select(time, condition)
#' se_mini_tidybulk |> bind_cols(tt_bind)
#'
#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom ttservice bind_rows
#' @export
bind_rows.tidybulk <- function(..., .id = NULL) {
	tts <- flatten_if(dots_values(...), is_spliced) # Original that fails Bioconductor dplyr:::flatten_bindable(rlang::dots_values(...))

	par1 <- tts[[1]] |> get_tt_columns() |> unlist()
	par2 <- tts[[2]] |> get_tt_columns() |> unlist()

	# # tt_columns of the two objects must match
	# error_if_parameters_not_match(par1, par2)

	ttservice:::bind_rows.data.frame(..., .id = .id) |>
		# Attach attributes
		reattach_internals(tts[[1]])
}

#' @rdname bind_rows
#' @aliases bind_cols
#' 
#' @importFrom ttservice bind_cols
#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @export
bind_cols.tidybulk <- function(..., .id = NULL) {
	tts <- tts <- flatten_if(dots_values(...), is_spliced) # Original that fails Bioconductor dplyr:::flatten_bindable(rlang::dots_values(...))

	dplyr::bind_cols(..., .id = .id) |>
		# Attach attributes
		reattach_internals(tts[[1]])
}

#' @name distinct
#' @rdname distinct
#' @inherit dplyr::distinct
#' 
#' @examples
#' tidybulk::se_mini |> tidybulk() |> distinct()
#'
#' @importFrom dplyr distinct
#' @export
distinct.tidybulk <- function (.data, ..., .keep_all = FALSE) {
	.data |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::distinct(..., .keep_all = .keep_all) |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name filter
#' @rdname filter
#' @inherit dplyr::filter
#' 
#' @examples
#' data(se)
#' se |> tidybulk() |> filter(dex=="untrt")
#' 
#' # Learn more in ?dplyr_tidy_eval
#' 
#' @importFrom dplyr filter
#' @export
filter.tidybulk <- function (.data, ..., .preserve = FALSE) {
	.data |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::filter( ..., .preserve = .preserve) |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name group_by
#' @rdname group_by
#' @inherit dplyr::group_by
#' @seealso \code{}
#' @importFrom dplyr group_by group_by_drop_default
#' @export
group_by.tidybulk <- function (.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)) {
	.data |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::group_by( ..., .drop = .drop) |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name ungroup
#' @rdname ungroup
#' @inherit dplyr::ungroup
#' @importFrom dplyr ungroup
#' @export
ungroup.tidybulk <- function (x, ...) {
	x |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::ungroup( ...) |>

		# Attach attributes
		reattach_internals(x) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}


#' @name summarise
#' @rdname summarise
#' @aliases summarize
#' @inherit dplyr::summarise
#' @family single table verbs
#' 
#' @importFrom dplyr summarise
#' @importFrom purrr map
#' @export
summarise.tidybulk <- function (.data, ...) {
	.data |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::summarise( ...) |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name summarise
#' @rdname summarise
#' @importFrom dplyr summarize
#' @export
summarize.tidybulk <- summarise.tidybulk

#' @name mutate
#' @rdname mutate
#' @inherit dplyr::mutate
#' @family single table verbs
#' @importFrom dplyr mutate
#' @export
mutate.tidybulk <- function(.data, ...) {
	.data |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::mutate(...) |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @inherit dplyr::mutate
#' @export
mutate.nested_tidybulk <- function(.data, ...) {
	.data |>
		drop_class(c("nested_tidybulk", "tt")) |>
		dplyr::mutate(...) |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("nested_tidybulk")
}

#' @name rename
#' @rdname rename
#' @inherit dplyr::rename
#' @family single table verbs
#' @importFrom dplyr rename
#' @export
rename.tidybulk <- function(.data, ...) {
	.data |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::rename(...) |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name rowwise
#' @rdname rowwise
#' @inherit dplyr::rowwise
#' @examples
#' 
#' df <- expand.grid(x = 1:3, y = 3:1)
#' df_done <- df |> rowwise() 
#'
#' @importFrom dplyr rowwise
#' @export
rowwise.tidybulk <- function(data, ...) {
	data |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::rowwise() |>

		# Attach attributes
		reattach_internals(.data) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name left_join
#' @rdname left_join
#' @inherit dplyr::left_join
#' @examples
#' annotation <- tidybulk::se_mini |> tidybulk() |> as_tibble() |> distinct(.sample) |> mutate(source = "AU")
#' tidybulk::se_mini |> tidybulk() |> as_tibble() |> left_join(annotation)
#'
#' @importFrom dplyr left_join
#' @export
left_join.tidybulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
															...) {
	x |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::left_join(y, by = by, copy = copy, suffix = suffix, ...) |>

		# Attach attributes
		reattach_internals(x) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name inner_join
#' @rdname inner_join
#' @inherit dplyr::inner_join
#'
#' @examples
#' annotation <- tidybulk::se_mini |> tidybulk() |> as_tibble() |> distinct(.sample) |> mutate(source = "AU")
#' tidybulk::se_mini |> tidybulk() |> as_tibble() |> inner_join(annotation)
#'
#' @importFrom dplyr inner_join
#' @export
inner_join.tidybulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...) {
	x |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::inner_join(y, by = by, copy = copy, suffix = suffix, ...) |>

		# Attach attributes
		reattach_internals(x) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @name right_join
#' @rdname right_join
#' @inherit dplyr::right_join
#' 
#' @examples
#' annotation <- tidybulk::se_mini |> tidybulk() |> as_tibble() |> distinct(.sample) |> mutate(source = "AU")
#' tidybulk::se_mini |> tidybulk() |> as_tibble() |> right_join(annotation)
#'
#' @importFrom dplyr right_join
#' @export
right_join.tidybulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
															 ...) {
	x |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::right_join(y, by = by, copy = copy, suffix = suffix, ...) |>

		# Attach attributes
		reattach_internals(x) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}


#' @name full_join
#' @rdname full_join
#' @inherit dplyr::full_join
#' 
#' @examples
#' annotation <- tidybulk::se_mini |> tidybulk() |> as_tibble() |> distinct(.sample) |> mutate(source = "AU")
#' tidybulk::se_mini |> tidybulk() |> as_tibble() |> full_join(annotation)
#'
#' @importFrom dplyr full_join
#' @export
full_join.tidybulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
															 ...) {
	x |>
		drop_class(c("tidybulk", "tt")) |>
		dplyr::full_join(y, by = by, copy = copy, suffix = suffix, ...) |>

		# Attach attributes
		reattach_internals(x) |>

		# Add class
		add_class("tt") |>
		add_class("tidybulk")
}

#' @importFrom dplyr do
#' @export
dplyr::do
