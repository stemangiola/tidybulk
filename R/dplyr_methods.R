
#' Arrange rows by column values
#'
#'
#' @description
#' `arrange()` order the rows of a data frame rows by the values of selected
#' columns.
#'
#' Unlike other dplyr verbs, `arrange()` largely ignores grouping; you
#' need to explicit mention grouping variables (or use  `by_group = TRUE`)
#' in order to group by them, and functions of variables are evaluated
#' once per data frame, not once per group.
#'
#' @details
#' ## Locales
#' The sort order for character vectors will depend on the collating sequence
#' of the locale in use: see [locales()].
#'
#' ## Missing values
#' Unlike base sorting with `sort()`, `NA` are:
#' * always sorted to the end for local data, even when wrapped with `desc()`.
#' * treated differently for remote data, depending on the backend.
#'
#' @return
#' An object of the same type as `.data`.
#'
#' * All rows appear in the output, but (usually) in a different place.
#' * Columns are not modified.
#' * Groups are not modified.
#' * Data frame attributes are preserved.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @export
#' @param .data A data frame, data frame extension (e.g. a tibble), or a
#'   lazy data frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for
#'   more details.
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Variables, or functions or
#'   variables. Use [desc()] to sort a variable in descending order.
#' @family single table verbs
#' @examples
#' `%>%` = magrittr::`%>%`
#' arrange(mtcars, cyl, disp)
arrange <- function(.data, ..., .by_group = FALSE) {
	UseMethod("arrange")
}

#' @param .by_group If `TRUE`, will sort first by grouping variable. Applies to
#'   grouped data frames only.
#' @rdname dplyr-methods
#' @export
#'
############# START ADDED tidybulk ###################################

arrange.default <- function(.data, ..., .by_group = FALSE) {

	dplyr::arrange(.data, ..., .by_group = .by_group)

}

#' @export
arrange.tidybulk <- function(.data, ..., .by_group = FALSE) {

	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::arrange( ..., .by_group = .by_group) %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}

############# END ADDED tidybulk #####################################

#' Efficiently bind multiple data frames by row and column
#'
#' This is an efficient implementation of the common pattern of
#' `do.call(rbind, dfs)` or `do.call(cbind, dfs)` for binding many
#' data frames into one.
#'
#' The output of `bind_rows()` will contain a column if that column
#' appears in any of the inputs.
#'
#' @param ... Data frames to combine.
#'
#'   Each argument can either be a data frame, a list that could be a data
#'   frame, or a list of data frames.
#'
#'   When row-binding, columns are matched by name, and any missing
#'   columns will be filled with NA.
#'
#'   When column-binding, rows are matched by position, so all data
#'   frames must have the same number of rows. To match by value, not
#'   position, see [mutate-joins].
#' @param .id Data frame identifier.
#'
#'   When `.id` is supplied, a new column of identifiers is
#'   created to link each row to its original data frame. The labels
#'   are taken from the named arguments to `bind_rows()`. When a
#'   list of data frames is supplied, the labels are taken from the
#'   names of the list. If no names are found a numeric sequence is
#'   used instead.
#' @return `bind_rows()` and `bind_cols()` return the same type as
#'   the first input, either a data frame, `tbl_df`, or `grouped_df`.
#' @examples
#' `%>%` = magrittr::`%>%`
#' one <- mtcars[1:4, ]
#' two <- mtcars[11:14, ]
#'
#' # You can supply data frames as arguments:
#' bind_rows(one, two)
#'
#' @name bind
NULL


#' @rdname dplyr-methods
#' @export
#'
############# START ADDED tidybulk #####################################

bind_rows <- function(..., .id = NULL) {
	UseMethod("bind_rows")
}

#' @export
bind_rows.default <-  function(..., .id = NULL)
{
	dplyr::bind_rows(..., .id = .id)
}

#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' 
#' @export
#' 
bind_rows.tidybulk <- function(..., .id = NULL)
{

	tts = flatten_if(dots_values(...), is_spliced) # Original that fails Bioconductor dplyr:::flatten_bindable(rlang::dots_values(...))

	par1 = tts[[1]] %>% get_tt_columns() %>% unlist
	par2 = tts[[2]] %>% get_tt_columns() %>% unlist

	# # tt_columns of the two objects must match
	# error_if_parameters_not_match(par1, par2)

	dplyr::bind_rows(..., .id = .id) %>%

		# Attach attributes
		reattach_internals(tts[[1]])

}

############# END ADDED tidybulk #####################################

#' @export
#' @rdname dplyr-methods
############# START ADDED tidybulk #####################################

bind_cols <- function(..., .id = NULL) {
	UseMethod("bind_cols")
}

#' @export
bind_cols.default <-  function(..., .id = NULL)
{
	dplyr::bind_cols(..., .id = .id)
}

#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' 
#' @export
#' 
bind_cols.tidybulk <- function(..., .id = NULL)
{

	tts = 	tts = flatten_if(dots_values(...), is_spliced) # Original that fails Bioconductor dplyr:::flatten_bindable(rlang::dots_values(...))

	dplyr::bind_cols(..., .id = .id) %>%

		# Attach attributes
		reattach_internals(tts[[1]])

}

############# END ADDED tidybulk #####################################
############# START ADDED tidybulk #####################################

# #' @importFrom dplyr arrange_all
# #' @export
# dplyr::arrange_all
#
# #' @importFrom dplyr arrange_at
# #' @export
# dplyr::arrange_at
#
# #' @importFrom dplyr arrange_if
# #' @export
# dplyr::arrange_if

############# END ADDED tidybulk #####################################
############# START ADDED tidybulk #####################################

#' distinct
#' @param .data A tbl. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#' @param .keep_all If TRUE, keep all variables in .data. If a combination of ... is not distinct, this keeps the first row of values. (See dplyr)
#'
#' @return A tt object
#'
#' @examples
#'
#' distinct(tidybulk::counts_mini)
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
distinct.tidybulk <- function (.data, ..., .keep_all = FALSE)
{
	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::distinct(..., .keep_all = .keep_all) %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}
############# END ADDED tidybulk #####################################

############# START ADDED tidybulk #####################################

# #' @importFrom dplyr distinct_all
# #' @export
# dplyr::distinct_all
#
# #' @importFrom dplyr distinct_at
# #' @export
# dplyr::distinct_at
#
# #' @importFrom dplyr distinct_if
# #' @export
# dplyr::distinct_if

############# END ADDED tidybulk #####################################

#' Subset rows using column values
#'
#' `filter()` retains the rows where the conditions you provide a `TRUE`. Note
#' that, unlike base subsetting with `[`, rows where the condition evaluates
#' to `NA` are dropped.
#'
#' dplyr is not yet smart enough to optimise filtering optimisation
#' on grouped datasets that don't need grouped calculations. For this reason,
#' filtering is often considerably faster on [ungroup()]ed data.
#'
#' @section Useful filter functions:
#'
#' * [`==`], [`>`], [`>=`] etc
#' * [`&`], [`|`], [`!`], [xor()]
#' * [is.na()]
#' * [between()], [near()]
#'
#' @section Grouped tibbles:
#'
#' Because filtering expressions are computed within groups, they may
#' yield different results on grouped tibbles. This will be the case
#' as soon as an aggregating, lagging, or ranking function is
#' involved. Compare this ungrouped filtering:
#'
#'
#' The former keeps rows with `mass` greater than the global average
#' whereas the latter keeps rows with `mass` greater than the gender
#'
#' average.
#' @family single table verbs
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Logical predicates defined in
#'   terms of the variables in `.data`.
#'   Multiple conditions are combined with `&`. Only rows where the
#'   condition evaluates to `TRUE` are kept.
#' @param .preserve when `FALSE` (the default), the grouping structure
#'   is recalculated based on the resulting data, otherwise it is kept as is.
#' @return
#' An object of the same type as `.data`.
#'
#' * Rows are a subset of the input, but appear in the same order.
#' * Columns are not modified.
#' * The number of groups may be reduced (if `.preserve` is not `TRUE`).
#' * Data frame attributes are preserved.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @seealso [filter_all()], [filter_if()] and [filter_at()].
#' @export
#' @examples
#'
#' # Learn more in ?dplyr_tidy_eval
############# START ADDED tidybulk #####################################
#' @export
filter <- function (.data, ..., .preserve = FALSE)  {
	UseMethod("filter")
}

#' @export
filter.default <-  function (.data, ..., .preserve = FALSE)
{
	dplyr::filter(.data, ..., .preserve = .preserve)
}

#' @export
filter.tidybulk <- function (.data, ..., .preserve = FALSE)
{
	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::filter( ..., .preserve = .preserve) %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}
############# END ADDED tidybulk #####################################


############# START ADDED tidybulk #####################################

# #' @importFrom dplyr filter_all
# #' @export
# dplyr::filter_all
#
# #' @importFrom dplyr filter_at
# #' @export
# dplyr::filter_at
#
# #' @importFrom dplyr filter_if
# #' @export
# dplyr::filter_if

############# END ADDED tidybulk #####################################

#' Group by one or more variables
#'
#' @description
#' Most data operations are done on groups defined by variables.
#' `group_by()` takes an existing tbl and converts it into a grouped tbl
#' where operations are performed "by group". `ungroup()` removes grouping.
#'
#' @family grouping functions
#' @inheritParams arrange
#' @param ... In `group_by()`, variables or computations to group by.
#'   In `ungroup()`, variables to remove from the grouping.
#' @param .add When `FALSE`, the default, `group_by()` will
#'   override existing groups. To add to the existing groups, use
#'   `.add = TRUE`.
#'
#'   This argument was previously called `add`, but that prevented
#'   creating a new grouping variable called `add`, and conflicts with
#'   our naming conventions.
#' @param .drop When `.drop = TRUE`, empty groups are dropped. See [group_by_drop_default()] for
#'   what the default value is for this argument.
#' @return A [grouped data frame][grouped_df()], unless the combination of `...` and `add`
#'   yields a non empty set of grouping columns, a regular (ungrouped) data frame
#'   otherwise.
#' @section Methods:
#' These function are **generic**s, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' Methods available in currently loaded packages:
#'
#' @export
#' @examples
#' `%>%` = magrittr::`%>%`
#' by_cyl <- mtcars %>% group_by(cyl)
#'

############# START ADDED tidybulk #####################################
#' @export
group_by <- function (.data, ..., .add = FALSE, .drop = group_by_drop_default(.data))  {
	UseMethod("group_by")
}

#' @export
group_by.default <-  function (.data, ..., .add = FALSE, .drop = group_by_drop_default(.data))
{
	dplyr::group_by(.data, ...,  .drop = .drop)
}

#' @export
group_by.tidybulk <- function (.data, ..., .add = FALSE, .drop = group_by_drop_default(.data))
{
	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::group_by( ..., .drop = .drop) %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}
############# END ADDED tidybulk #####################################


#' @rdname dplyr-methods
#' @export
#' @param x A [tbl()]
ungroup <- function(x, ...) {
	UseMethod("ungroup")
}
############# START ADDED tidybulk #####################################

#' @export
ungroup.default <-  function (x, ...)
{
	dplyr::ungroup(x, ...)
}

#' @export
ungroup.tidybulk <- function (x, ...)
{
	x %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::ungroup( ...) %>%

		# Attach attributes
		reattach_internals(x) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}
############# END ADDED tidybulk #####################################

############# START ADDED tidybulk #####################################

# #' @importFrom dplyr group_by_all
# #' @export
# dplyr::group_by_all
#
# #' @importFrom dplyr group_by_at
# #' @export
# dplyr::group_by_at
#
# #' @importFrom dplyr group_by_if
# #' @export
# dplyr::group_by_if

############# END ADDED tidybulk #####################################

#' Summarise each group to fewer rows
#'
#' @description
#' `summarise()` creates a new data frame. It will have one (or more) rows for
#' each combination of grouping variables; if there are no grouping variables,
#' the output will have a single row summarising all observations in the input.
#' It will contain one column for each grouping variable and one column
#' for each of the summary statistics that you have specified.
#'
#' `summarise()` and `summarize()` are synonyms.
#'
#' @section Useful functions:
#'
#' * Center: [mean()], [median()]
#' * Spread: [sd()], [IQR()], [mad()]
#' * Range: [min()], [max()], [quantile()]
#' * Position: [first()], [last()], [nth()],
#' * Count: [n()], [n_distinct()]
#' * Logical: [any()], [all()]
#'
#' @section Backend variations:
#'
#' The data frame backend supports creating a variable and using it in the
#' same summary. This means that previously created summary variables can be
#' further transformed or combined within the summary, as in [mutate()].
#' However, it also means that summary variables with the same names as previous
#' variables overwrite them, making those variables unavailable to later summary
#' variables.
#'
#' This behaviour may not be supported in other backends. To avoid unexpected
#' results, consider using new names for your summary variables, especially when
#' creating multiple summaries.
#'
#' @export
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Name-value pairs of summary
#'   functions. The name will be the name of the variable in the result.
#'
#'   The value can be:
#'
#'   * A vector of length 1, e.g. `min(x)`, `n()`, or `sum(is.na(y))`.
#'   * A vector of length `n`, e.g. `quantile()`.
#'   * A data frame, to add multiple columns from a single expression.
#' @family single table verbs
#' @return
#' An object _usually_ of the same type as `.data`.
#'
#' * The rows come from the underlying `group_keys()`.
#' * The columns are a combination of the grouping keys and the summary
#'   expressions that you provide.
#' * If `x` is grouped by more than one variable, the output will be another
#'   [grouped_df] with the right-most group removed.
#' * If `x` is grouped by one variable, or is not grouped, the output will
#'   be a [tibble].
#' * Data frame attributes are **not** preserved, because `summarise()`
#'   fundamentally creates a new data frame.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @examples
#' `%>%` = magrittr::`%>%`
#' # A summary applied to ungrouped tbl returns a single row
#' mtcars %>%
#'   summarise(mean = mean(disp))
#'
############# START ADDED tidybulk #####################################
#' @export
summarise <- function (.data, ...)  {
	UseMethod("summarise")
}

#' @export
summarise.default <-  function (.data, ...)
{
	dplyr::summarise(.data, ...)
}

#' @export
summarise.tidybulk <- function (.data, ...)
{
	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::summarise( ...) %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}
############# END ADDED tidybulk #####################################

############# START ADDED tidybulk #####################################

# #' @importFrom dplyr summarize_all
# #' @export
# dplyr::summarize_all
#
# #' @importFrom dplyr summarize_at
# #' @export
# dplyr::summarize_at
#
# #' @importFrom dplyr summarize_if
# #' @export
# dplyr::summarize_if

############# END ADDED tidybulk #####################################

# #' @rdname dplyr-methods
# #' @export
# summarize_all <- summarise_all
# #' @rdname dplyr-methods
# #' @export
# summarize_if <- summarise_if
# #' @rdname dplyr-methods
# #' @export
# summarize_at <- summarise_at

#' Create, modify, and delete columns
#'
#' `mutate()` adds new variables and preserves existing ones;
#' `transmute()` adds new variables and drops existing ones.
#' New variables overwrite existing variables of the same name.
#' Variables can be removed by setting their value to `NULL`.
#'
#' @section Useful mutate functions:
#'
#' * [`+`], [`-`], [log()], etc., for their usual mathematical meanings
#'
#' * [lead()], [lag()]
#'
#' * [dense_rank()], [min_rank()], [percent_rank()], [row_number()],
#'   [cume_dist()], [ntile()]
#'
#' * [cumsum()], [cummean()], [cummin()], [cummax()], [cumany()], [cumall()]
#'
#' * [na_if()], [coalesce()]
#'
#' * [if_else()], [recode()], [case_when()]
#'
#' @section Grouped tibbles:
#'
#' Because mutating expressions are computed within groups, they may
#' yield different results on grouped tibbles. This will be the case
#' as soon as an aggregating, lagging, or ranking function is
#' involved. Compare this ungrouped mutate:
#'

#' With the grouped equivalent:
#'
#' The former normalises `mass` by the global average whereas the
#' latter normalises by the averages within gender levels.
#'
#' @export
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Name-value pairs.
#'   The name gives the name of the column in the output.
#'
#'   The value can be:
#'
#'   * A vector of length 1, which will be recycled to the correct length.
#'   * A vector the same length as the current group (or the whole data frame
#'     if ungrouped).
#'   * `NULL`, to remove the column.
#'   * A data frame or tibble, to create multiple columns in the output.
#'
#' @family single table verbs
#' @return
#' An object of the same type as `.data`.
#'
#' For `mutate()`:
#'
#' * Rows are not affected.
#' * Existing columns will be preserved unless explicitly modified.
#' * New columns will be added to the right of existing columns.
#' * Columns given value `NULL` will be removed
#' * Groups will be recomputed if a grouping variable is mutated.
#' * Data frame attributes are preserved.
#'
#' For `transmute()`:
#'
#' * Rows are not affected.
#' * Apart from grouping variables, existing columns will be remove unless
#'   explicitly kept.
#' * Column order matches order of expressions.
#' * Groups will be recomputed if a grouping variable is mutated.
#' * Data frame attributes are preserved.
#' @section Methods:
#' These function are **generic**s, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' Methods available in currently loaded packages:
#'
#' @examples
#' `%>%` = magrittr::`%>%`
#' # Newly created variables are available immediately
#' mtcars %>% as_tibble() %>% mutate(
#'   cyl2 = cyl * 2,
#'   cyl4 = cyl2 * 2
#' )
#'
############# START ADDED tidybulk #####################################
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
mutate.tidybulk <- function(.data, ...)
{
	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::mutate(...) %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")


}
############# END ADDED tidybulk #####################################

############# START ADDED tidybulk #####################################

# #' @importFrom dplyr mutate_all
# #' @export
# dplyr::mutate_all
#
# #' @importFrom dplyr mutate_at
# #' @export
# dplyr::mutate_at
#
# #' @importFrom dplyr mutate_if
# #' @export
# dplyr::mutate_if

############# END ADDED tidybulk #####################################

#' Rename columns
#'
#' Rename individual variables using `new_name = old_name` syntax.
#'
#' @section Scoped selection and renaming:
#'
#' Use the three scoped variants ([rename_all()], [rename_if()], [rename_at()])
#' to renaming a set of variables with a function.
#'
#' @inheritParams arrange
#' @param ... <[`tidy-select`][dplyr_tidy_select]> Use `new_name = old_name`
#'   to rename selected variables.
#' @return
#' An object of the same type as `.data`.
#' * Rows are not affected.
#' * Column names are changed; column order is preserved
#' * Data frame attributes are preserved.
#' * Groups are updated to reflect new names.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @family single table verbs
#' @export
#' @examples
#' `%>%` = magrittr::`%>%`
#' iris <- as_tibble(iris) # so it prints a little nicer
#' rename(iris, petal_length = Petal.Length)
############# START ADDED tidybulk #####################################
#' @export
rename <- function(.data, ...) {
	UseMethod("rename")
}

#' @export
rename.default <-  function(.data, ...)
{
	dplyr::rename(.data, ...)
}

#' @export
rename.tidybulk <- function(.data, ...)
{
	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::rename(...) %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")


}
############# END ADDED tidybulk #####################################

#' Group input by rows
#'
#'
#' See [this repository](https://github.com/jennybc/row-oriented-workflows)
#' for alternative ways to perform row-wise operations.
#'
#' `rowwise()` is used for the results of [do()] when you
#' create list-variables. It is also useful to support arbitrary
#' complex operations that need to be applied to each row.
#'
#' Currently, rowwise grouping only works with data frames. Its
#' main impact is to allow you to work with list-variables in
#' [summarise()] and [mutate()] without having to
#' use \code{[[1]]}. This makes `summarise()` on a rowwise tbl
#' effectively equivalent to [plyr::ldply()].
#'
#' @param .data Input data frame.
#'
#' @return A `tbl`
#'
#'   A `tbl`
#'
#' @export
#' @examples
#' `%>%` = magrittr::`%>%`
#' df <- expand.grid(x = 1:3, y = 3:1)
#' df_done <- df %>% rowwise() %>% do(i = seq(.$x, .$y))
############# START ADDED tidybulk #####################################
#' @export
rowwise <- function(.data) {
	UseMethod("rowwise")
}

#' @export
rowwise.default <-  function(.data)
{
	dplyr::rowwise(.data)
}

#' @export
rowwise.tidybulk <- function(.data)
{
	.data %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::rowwise() %>%

		# Attach attributes
		reattach_internals(.data) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")


}
############# END ADDED tidybulk #####################################



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
#'
#' @examples
#'`%>%` = magrittr::`%>%`
#' annotation = tidybulk::counts %>% distinct(sample) %>% mutate(source = "AU")
#' tidybulk::counts %>% left_join(annotation)
#'
left_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
	UseMethod("left_join")
}

#' @export
left_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
																...)
{
	dplyr::left_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}

#' @export
left_join.tidybulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
															...)
{
	x %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::left_join(y, by = by, copy = copy, suffix = suffix, ...) %>%

		# Attach attributes
		reattach_internals(x) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

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
#' @examples
#'`%>%` = magrittr::`%>%`
#' annotation = tidybulk::counts %>% distinct(sample) %>% mutate(source = "AU")
#' tidybulk::counts %>% inner_join(annotation)
#'
#' @export
inner_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
	UseMethod("inner_join")
}

#' @export
inner_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),				 ...)
{
	dplyr::inner_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}

#' @export
inner_join.tidybulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)
{
	x %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::inner_join(y, by = by, copy = copy, suffix = suffix, ...) %>%

		# Attach attributes
		reattach_internals(x) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}

#' Right join datasets
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
#' @examples
#'`%>%` = magrittr::`%>%`
#' annotation = tidybulk::counts %>% distinct(sample) %>% mutate(source = "AU")
#' tidybulk::counts %>% right_join(annotation)
#'
#' @export
right_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
	UseMethod("right_join")
}

#' @export
right_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
																 ...)
{
	dplyr::right_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}

#' @export
right_join.tidybulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
															 ...)
{
	x %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::right_join(y, by = by, copy = copy, suffix = suffix, ...) %>%

		# Attach attributes
		reattach_internals(x) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}


#' Full join datasets
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
#' @examples
#'`%>%` = magrittr::`%>%`
#' annotation = tidybulk::counts %>% distinct(sample) %>% mutate(source = "AU")
#' tidybulk::counts %>% full_join(annotation)
#'
#' @export
full_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
	UseMethod("full_join")
}

#' @export
full_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
																 ...)
{
	dplyr::full_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}

#' @export
full_join.tidybulk <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
															 ...)
{
	x %>%
		drop_class(c("tidybulk", "tt")) %>%
		dplyr::full_join(y, by = by, copy = copy, suffix = suffix, ...) %>%

		# Attach attributes
		reattach_internals(x) %>%

		# Add class
		add_class("tt") %>%
		add_class("tidybulk")

}

#' @importFrom dplyr do
#' @export
dplyr::do
