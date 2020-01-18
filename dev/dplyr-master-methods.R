
#' Arrange rows by column values
#'
#' See \code{dpyr::\link[dpyr:arrange]{arrange}} for details.
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
#' \Sexpr[stage=render,results=Rd]{dplyr:::methods_rd("arrange")}.
#' @export
#' @param .data A data frame, data frame extension (e.g. a tibble), or a
#'   lazy data frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for
#'   more details.
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Variables, or functions or
#'   variables. Use [desc()] to sort a variable in descending order.
#' @family single table verbs
#' @examples
#' arrange(mtcars, cyl, disp)
#' arrange(mtcars, desc(disp))
#'
#' # grouped arrange ignores groups
#' by_cyl <- mtcars %>% group_by(cyl)
#' by_cyl %>% arrange(desc(wt))
#' # Unless you specifically ask:
#' by_cyl %>% arrange(desc(wt), .by_group = TRUE)
arrange <- function(.data, ..., .by_group = FALSE) {
  UseMethod("arrange")
}

#' @param .by_group If `TRUE`, will sort first by grouping variable. Applies to
#'   grouped data frames only.
#' @rdname arrange
#' @export
#'
############# START ADDED TTBULK ###################################

arrange.default <- function(.data, ..., .by_group = FALSE) {

  dplyr::arrange(.data, ..., .by_group = .by_group)

}

#' @export
arrange.ttBulk <- function(.data, ..., .by_group = FALSE) {

  .data %>%
    drop_class(c("ttBulk", "tt")) %>%
    dplyr::arrange(.data, ..., .by_group = .by_group) %>%

    # Attach attributes
    add_attr(.data %>% attr("parameters"), "parameters") %>%

    # Add class
    add_class("tt") %>%
    add_class("ttBulk")

}

############# END ADDED TTBULK #####################################

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
#' one <- mtcars[1:4, ]
#' two <- mtcars[11:14, ]
#'
#' # You can supply data frames as arguments:
#' bind_rows(one, two)
#'
#' # The contents of lists are spliced automatically:
#' bind_rows(list(one, two))
#' bind_rows(split(mtcars, mtcars$cyl))
#' bind_rows(list(one, two), list(two, one))
#'
#'
#' # In addition to data frames, you can supply vectors. In the rows
#' # direction, the vectors represent rows and should have inner
#' # names:
#' bind_rows(
#'   c(a = 1, b = 2),
#'   c(a = 3, b = 4)
#' )
#'
#' # You can mix vectors and data frames:
#' bind_rows(
#'   c(a = 1, b = 2),
#'   tibble(a = 3:4, b = 5:6),
#'   c(a = 7, b = 8)
#' )
#'
#'
#' # Note that for historical reasons, lists containing vectors are
#' # always treated as data frames. Thus their vectors are treated as
#' # columns rather than rows, and their inner names are ignored:
#' ll <- list(
#'   a = c(A = 1, B = 2),
#'   b = c(A = 3, B = 4)
#' )
#' bind_rows(ll)
#'
#' # You can circumvent that behaviour with explicit splicing:
#' bind_rows(!!!ll)
#'
#'
#' # When you supply a column name with the `.id` argument, a new
#' # column is created to link each row to its original data frame
#' bind_rows(list(one, two), .id = "id")
#' bind_rows(list(a = one, b = two), .id = "id")
#' bind_rows("group 1" = one, "group 2" = two, .id = "groups")
#'
#' # Columns don't need to match when row-binding
#' bind_rows(data.frame(x = 1:3), data.frame(y = 1:4))
#' \dontrun{
#' # Rows do need to match when column-binding
#' bind_cols(data.frame(x = 1), data.frame(y = 1:2))
#' }
#'
#' bind_cols(one, two)
#' bind_cols(list(one, two))
#' @name bind
NULL

#' @export
#' @rdname bind
#' @export
#'
############# START ADDED TTBULK #####################################

bind_rows <- function(..., .id = NULL) {
  UseMethod("bind_rows")
}

#' @export
bind_rows.default <-  function(..., .id = NULL)
{
  dplyr::bind_rows(..., .id = .id)
}

#' @export
bind_rows.ttBulk <- function(..., .id = NULL)
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


  dplyr::bind_rows(..., .id = .id) %>%

    # Attach attributes
    add_attr(par, "parameters")

}

############# END ADDED TTBULK #####################################

#' @export
#' @rdname bind
############# START ADDED TTBULK #####################################

bind_cols <- function(..., .id = NULL) {
  UseMethod("bind_cols")
}

#' @export
bind_cols.default <-  function(..., .id = NULL)
{
  dplyr::bind_cols(..., .id = .id)
}

#' @export
bind_cols.ttBulk <- function(..., .id = NULL)
{

  tts = dplyr:::flatten_bindable(rlang::dots_values(...))

  dplyr::bind_cols(..., .id = .id) %>%

    # Attach attributes
    add_attr(tts[[1]] %>% attr("parameters") , "parameters")

}

############# END ADDED TTBULK #####################################
############# START ADDED TTBULK #####################################

#' @importFrom dplyr arrange_all
#' @export
dplyr::arrange_all

#' @importFrom dplyr arrange_at
#' @export
dplyr::arrange_at

#' @importFrom dplyr arrange_if
#' @export
dplyr::arrange_if

############# END ADDED TTBULK #####################################
############# START ADDED TTBULK #####################################

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
############# END ADDED TTBULK #####################################

############# START ADDED TTBULK #####################################

#' @importFrom dplyr distinct_all
#' @export
dplyr::distinct_all

#' @importFrom dplyr distinct_at
#' @export
dplyr::distinct_at

#' @importFrom dplyr distinct_if
#' @export
dplyr::distinct_if

############# END ADDED TTBULK #####################################

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
#' ```
#' starwars %>% filter(mass > mean(mass, na.rm = TRUE))
#' ```
#'
#' With the grouped equivalent:
#'
#' ```
#' starwars %>% group_by(gender) %>% filter(mass > mean(mass, na.rm = TRUE))
#' ```
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
#' \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("filter")}.
#' @seealso [filter_all()], [filter_if()] and [filter_at()].
#' @export
#' @examples
#' filter(starwars, species == "Human")
#' filter(starwars, mass > 1000)
#'
#' # Multiple criteria
#' filter(starwars, hair_color == "none" & eye_color == "black")
#' filter(starwars, hair_color == "none" | eye_color == "black")
#'
#' # Multiple arguments are equivalent to and
#' filter(starwars, hair_color == "none", eye_color == "black")
#'
#'
#' # The filtering operation may yield different results on grouped
#' # tibbles because the expressions are computed within groups.
#' #
#' # The following filters rows where `mass` is greater than the
#' # global average:
#' starwars %>% filter(mass > mean(mass, na.rm = TRUE))
#'
#' # Whereas this keeps rows with `mass` greater than the gender
#' # average:
#' starwars %>% group_by(gender) %>% filter(mass > mean(mass, na.rm = TRUE))
#'
#'
#' # Refer to column names stored as strings with the `.data` pronoun:
#' vars <- c("mass", "height")
#' cond <- c(80, 150)
#' starwars %>%
#'   filter(
#'     .data[[vars[[1]]]] > cond[[1]],
#'     .data[[vars[[2]]]] > cond[[2]]
#'   )
#' # Learn more in ?dplyr_tidy_eval
############# START ADDED TTBULK #####################################
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
filter.ttBulk <- function (.data, ..., .preserve = FALSE)
{
  .data %>%
    drop_class(c("ttBulk", "tt")) %>%
    dplyr::filter(.data, ..., .preserve = .preserve) %>%

    # Attach attributes
    add_attr(.data %>% attr("parameters"), "parameters") %>%

    # Add class
    add_class("tt") %>%
    add_class("ttBulk")

}
############# END ADDED TTBULK #####################################


############# START ADDED TTBULK #####################################

#' @importFrom dplyr filter_all
#' @export
dplyr::filter_all

#' @importFrom dplyr filter_at
#' @export
dplyr::filter_at

#' @importFrom dplyr filter_if
#' @export
dplyr::filter_if

############# END ADDED TTBULK #####################################

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
#' * `group_by()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("group_by")}.
#' * `ungroup()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("ungroup")}.
#' @export
#' @examples
#' by_cyl <- mtcars %>% group_by(cyl)
#'
#' # grouping doesn't change how the data looks (apart from listing
#' # how it's grouped):
#' by_cyl
#'
#' # It changes how it acts with the other dplyr verbs:
#' by_cyl %>% summarise(
#'   disp = mean(disp),
#'   hp = mean(hp)
#' )
#' by_cyl %>% filter(disp == max(disp))
#'
#' # Each call to summarise() removes a layer of grouping
#' by_vs_am <- mtcars %>% group_by(vs, am)
#' by_vs <- by_vs_am %>% summarise(n = n())
#' by_vs
#' by_vs %>% summarise(n = sum(n))
#'
#' # To removing grouping, use ungroup
#' by_vs %>%
#'   ungroup() %>%
#'   summarise(n = sum(n))
#'
#' # You can group by expressions: this is just short-hand for
#' # a mutate() followed by a group_by()
#' mtcars %>% group_by(vsam = vs + am)
#'
#' # By default, group_by() overrides existing grouping
#' by_cyl %>%
#'   group_by(vs, am) %>%
#'   group_vars()
#'
#' # Use add = TRUE to instead append
#' by_cyl %>%
#'   group_by(vs, am, .add = TRUE) %>%
#'   group_vars()
#'
#'
#' # when factors are involved, groups can be empty
#' tbl <- tibble(
#'   x = 1:10,
#'   y = factor(rep(c("a", "c"), each  = 5), levels = c("a", "b", "c"))
#' )
#' tbl %>%
#'   group_by(y) %>%
#'   group_rows()
#'
############# START ADDED TTBULK #####################################
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
filter.ttBulk <- function (.data, ..., .preserve = FALSE)
{
  .data %>%
    drop_class(c("ttBulk", "tt")) %>%
    dplyr::filter(.data, ..., .preserve = .preserve) %>%

    # Attach attributes
    add_attr(.data %>% attr("parameters"), "parameters") %>%

    # Add class
    add_class("tt") %>%
    add_class("ttBulk")

}
############# END ADDED TTBULK #####################################
group_by <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)) {
  UseMethod("group_by")
}

#' @export
group_by.data.frame <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)) {
  groups <- group_by_prepare(.data, ..., .add = .add)
  grouped_df(groups$data, groups$group_names, .drop)
}

#' @rdname group_by
#' @export
#' @param x A [tbl()]
ungroup <- function(x, ...) {
  UseMethod("ungroup")
}

#' @export
ungroup.grouped_df <- function(x, ...) {
  if (missing(...)) {
    as_tibble(x)
  } else {
    old_groups <- group_vars(x)
    to_remove <- tidyselect::vars_select(names(x), ...)

    new_groups <- setdiff(old_groups, to_remove)
    group_by(x, !!!syms(new_groups))
  }
}

#' @export
ungroup.rowwise_df <- function(x, ...) {
  ellipsis::check_dots_empty()
  as_tibble(x)
}

#' @export
ungroup.data.frame <- function(x, ...) {
  ellipsis::check_dots_empty()
  x
}

#' Group by a selection of variables
#'
#' These [scoped] variants of [group_by()] group a data frame by a
#' selection of variables. Like [group_by()], they have optional
#' [mutate] semantics.
#'
#' @family grouping functions
#' @inheritParams scoped
#' @inheritParams group_by
#' @param .add See [group_by()]
#'
#' @export
#'
#' @section Grouping variables:
#'
#' Existing grouping variables are maintained, even if not included in
#' the selection.
#'
#' @examples
#' # Group a data frame by all variables:
#' group_by_all(mtcars)
#'
#' # Group by variables selected with a predicate:
#' group_by_if(iris, is.factor)
#'
#' # Group by variables selected by name:
#' group_by_at(mtcars, vars(vs, am))
#'
#' # Like group_by(), the scoped variants have optional mutate
#' # semantics. This provide a shortcut for group_by() + mutate():
#' d <- tibble(x=c(1,1,2,2), y=c(1,2,1,2))
#' group_by_all(d, as.factor)
#' group_by_if(iris, is.factor, as.character)
group_by_all <- function(.tbl, .funs = list(), ..., .add = FALSE, .drop = group_by_drop_default(.tbl)) {
  funs <- manip_all(.tbl, .funs, enquo(.funs), caller_env(), ...)
  if (!length(funs)) {
    funs <- syms(tbl_vars(.tbl))
  }
  .group_by_static_drop(.tbl, !!!funs, .add = .add, .drop = .drop)
}
#' @rdname group_by_all
#' @export
group_by_at <- function(.tbl, .vars, .funs = list(), ..., .add = FALSE, .drop = group_by_drop_default(.tbl)) {
  funs <- manip_at(.tbl, .vars, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  if (!length(funs)) {
    funs <- tbl_at_syms(.tbl, .vars, .include_group_vars = TRUE)
  }
  .group_by_static_drop(.tbl, !!!funs, .add = .add, .drop = .drop)
}
#' @rdname group_by_all
#' @export
group_by_if <- function(.tbl, .predicate, .funs = list(), ..., .add = FALSE, .drop = group_by_drop_default(.tbl)) {
  funs <- manip_if(.tbl, .predicate, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  if (!length(funs)) {
    funs <- tbl_if_syms(.tbl, .predicate, .include_group_vars = TRUE)
  }
  .group_by_static_drop(.tbl, !!!funs, .add = .add, .drop = .drop)
}
#' Summarise multiple columns
#'
#' @description
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' The [scoped] variants of [summarise()] make it easy to apply the same
#' transformation to multiple variables.
#' There are three variants.
#'  * `summarise_all()` affects every variable
#'  * `summarise_at()` affects variables selected with a character vector or
#'   vars()
#'  * `summarise_if()` affects variables selected with a predicate function
#'
#' @inheritParams scoped
#' @param .cols This argument has been renamed to `.vars` to fit
#'   dplyr's terminology and is deprecated.
#' @return A data frame. By default, the newly created columns have the shortest
#'   names needed to uniquely identify the output. To force inclusion of a name,
#'   even when not needed, name the input (see examples for details).
#' @seealso [The other scoped verbs][scoped], [vars()]
#'
#' @section Grouping variables:
#'
#' If applied on a grouped tibble, these operations are *not* applied
#' to the grouping variables. The behaviour depends on whether the
#' selection is **implicit** (`all` and `if` selections) or
#' **explicit** (`at` selections).
#'
#' * Grouping variables covered by explicit selections in
#'   `summarise_at()` are always an error. Add `-group_cols()` to the
#'   [vars()] selection to avoid this:
#'
#'   ```
#'   data %>%
#'     summarise_at(vars(-group_cols(), ...), myoperation)
#'   ```
#'
#'   Or remove `group_vars()` from the character vector of column names:
#'
#'   ```
#'   nms <- setdiff(nms, group_vars(data))
#'   data %>% summarise_at(nms, myoperation)
#'   ```
#'
#' * Grouping variables covered by implicit selections are silently
#'   ignored by `summarise_all()` and `summarise_if()`.
#'
#' @section Naming:
#'
#' The names of the new columns are derived from the names of the
#' input variables and the names of the functions.
#'
#' - if there is only one unnamed function (i.e. if `.funs` is an unnamed list
#'   of length one),
#'   the names of the input variables are used to name the new columns;
#'
#' - for `_at` functions, if there is only one unnamed variable (i.e.,
#'   if `.vars` is of the form `vars(a_single_column)`) and `.funs` has length
#'   greater than one,
#'   the names of the functions are used to name the new columns;
#'
#' - otherwise, the new names are created by
#'   concatenating the names of the input variables and the names of the
#'   functions, separated with an underscore `"_"`.
#'
#' The `.funs` argument can be a named or unnamed list.
#' If a function is unnamed and the name cannot be derived automatically,
#' a name of the form "fn#" is used.
#' Similarly, [vars()] accepts named and unnamed arguments.
#' If a variable in `.vars` is named, a new column by that name will be created.
#'
#' Name collisions in the new columns are disambiguated using a unique suffix.
#'
#' @section Life cycle:
#'
#' The functions are maturing, because the naming scheme and the
#' disambiguation algorithm are subject to change in dplyr 0.9.0.
#'
#' @examples
#' by_species <- iris %>%
#'   group_by(Species)
#'
#'
#' # The _at() variants directly support strings:
#' starwars %>%
#'   summarise_at(c("height", "mass"), mean, na.rm = TRUE)
#'
#' # You can also supply selection helpers to _at() functions but you have
#' # to quote them with vars():
#' starwars %>%
#'   summarise_at(vars(height:mass), mean, na.rm = TRUE)
#'
#' # The _if() variants apply a predicate function (a function that
#' # returns TRUE or FALSE) to determine the relevant subset of
#' # columns. Here we apply mean() to the numeric columns:
#' starwars %>%
#'   summarise_if(is.numeric, mean, na.rm = TRUE)
#'
#' # If you want to apply multiple transformations, pass a list of
#' # functions. When there are multiple functions, they create new
#' # variables instead of modifying the variables in place:
#' by_species %>%
#'   summarise_all(list(min, max))
#'
#' # Note how the new variables include the function name, in order to
#' # keep things distinct. Passing purrr-style lambdas often creates
#' # better default names:
#' by_species %>%
#'   summarise_all(list(~min(.), ~max(.)))
#'
#' # When that's not good enough, you can also supply the names explicitly:
#' by_species %>%
#'   summarise_all(list(min = min, max = max))
#'
#' # When there's only one function in the list, it modifies existing
#' # variables in place. Give it a name to create new variables instead:
#' by_species %>% summarise_all(list(med = median))
#' by_species %>% summarise_all(list(Q3 = quantile), probs = 0.75)
#' @export
summarise_all <- function(.tbl, .funs, ...) {
  funs <- manip_all(.tbl, .funs, enquo(.funs), caller_env(), ...)
  summarise(.tbl, !!!funs)
}
#' @rdname summarise_all
#' @export
summarise_if <- function(.tbl, .predicate, .funs, ...) {
  funs <- manip_if(.tbl, .predicate, .funs, enquo(.funs), caller_env(), ...)
  summarise(.tbl, !!!funs)
}
#' @rdname summarise_all
#' @export
summarise_at <- function(.tbl, .vars, .funs, ..., .cols = NULL) {
  .vars <- check_dot_cols(.vars, .cols)
  funs <- manip_at(.tbl, .vars, .funs, enquo(.funs), caller_env(), ...)
  summarise(.tbl, !!!funs)
}

#' @rdname summarise_all
#' @export
summarize_all <- summarise_all
#' @rdname summarise_all
#' @export
summarize_if <- summarise_if
#' @rdname summarise_all
#' @export
summarize_at <- summarise_at

#' Mutate multiple columns
#'
#' @description
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' The [scoped] variants of [mutate()] and [transmute()] make it easy to apply
#' the same transformation to multiple variables. There are three variants:
#'  * _all affects every variable
#'  * _at affects variables selected with a character vector or vars()
#'  * _if affects variables selected with a predicate function:
#'
#' @inheritParams scoped
#' @inheritParams summarise_all
#' @return A data frame. By default, the newly created columns have the shortest
#'   names needed to uniquely identify the output. To force inclusion of a name,
#'   even when not needed, name the input (see examples for details).
#' @seealso [The other scoped verbs][scoped], [vars()]
#'
#' @section Grouping variables:
#'
#' If applied on a grouped tibble, these operations are *not* applied
#' to the grouping variables. The behaviour depends on whether the
#' selection is **implicit** (`all` and `if` selections) or
#' **explicit** (`at` selections).
#'
#' * Grouping variables covered by explicit selections in
#'   `mutate_at()` and `transmute_at()` are always an error. Add
#'   `-group_cols()` to the [vars()] selection to avoid this:
#'
#'   ```
#'   data %>% mutate_at(vars(-group_cols(), ...), myoperation)
#'   ```
#'
#'   Or remove `group_vars()` from the character vector of column names:
#'
#'   ```
#'   nms <- setdiff(nms, group_vars(data))
#'   data %>% mutate_at(vars, myoperation)
#'   ```
#'
#' * Grouping variables covered by implicit selections are ignored by
#'   `mutate_all()`, `transmute_all()`, `mutate_if()`, and
#'   `transmute_if()`.
#'
#' @inheritSection summarise_all Naming
#' @inheritSection summarise_all Life cycle
#'
#' @examples
#' iris <- as_tibble(iris)
#'
#' # All variants can be passed functions and additional arguments,
#' # purrr-style. The _at() variants directly support strings. Here
#' # we'll scale the variables `height` and `mass`:
#' scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
#' starwars %>% mutate_at(c("height", "mass"), scale2)
#'
#' # You can pass additional arguments to the function:
#' starwars %>% mutate_at(c("height", "mass"), scale2, na.rm = TRUE)
#'
#' # You can also pass formulas to create functions on the spot, purrr-style:
#' starwars %>% mutate_at(c("height", "mass"), ~scale2(., na.rm = TRUE))
#'
#' # You can also supply selection helpers to _at() functions but you have
#' # to quote them with vars():
#' iris %>% mutate_at(vars(matches("Sepal")), log)
#'
#' # The _if() variants apply a predicate function (a function that
#' # returns TRUE or FALSE) to determine the relevant subset of
#' # columns. Here we divide all the numeric columns by 100:
#' starwars %>% mutate_if(is.numeric, scale2, na.rm = TRUE)
#'
#' # mutate_if() is particularly useful for transforming variables from
#' # one type to another
#' iris %>% mutate_if(is.factor, as.character)
#' iris %>% mutate_if(is.double, as.integer)
#'
#'
#' # Multiple transformations ----------------------------------------
#'
#' # If you want to apply multiple transformations, pass a list of
#' # functions. When there are multiple functions, they create new
#' # variables instead of modifying the variables in place:
#' iris %>% mutate_if(is.numeric, list(scale2, log))
#'
#' # The list can contain purrr-style formulas:
#' iris %>% mutate_if(is.numeric, list(~scale2(.), ~log(.)))
#'
#' # Note how the new variables include the function name, in order to
#' # keep things distinct. The default names are not always helpful
#' # but you can also supply explicit names:
#' iris %>% mutate_if(is.numeric, list(scale = scale2, log = log))
#'
#' # When there's only one function in the list, it modifies existing
#' # variables in place. Give it a name to instead create new variables:
#' iris %>% mutate_if(is.numeric, list(scale2))
#' iris %>% mutate_if(is.numeric, list(scale = scale2))
#' @export
mutate_all <- function(.tbl, .funs, ...) {
  check_grouped(.tbl, "mutate", "all", alt = TRUE)
  funs <- manip_all(.tbl, .funs, enquo(.funs), caller_env(), ...)
  mutate(.tbl, !!!funs)
}
#' @rdname mutate_all
#' @export
mutate_if <- function(.tbl, .predicate, .funs, ...) {
  check_grouped(.tbl, "mutate", "if")
  funs <- manip_if(.tbl, .predicate, .funs, enquo(.funs), caller_env(), ...)
  mutate(.tbl, !!!funs)
}
#' @rdname mutate_all
#' @export
mutate_at <- function(.tbl, .vars, .funs, ..., .cols = NULL) {
  .vars <- check_dot_cols(.vars, .cols)
  funs <- manip_at(.tbl, .vars, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  mutate(.tbl, !!!funs)
}

#' Count observations by group
#'
#' @description
#' `count()` lets you quickly count the unique values of a variable, or
#' unique combinations of multiple variables. `tally()` is a lower-level
#' helper that works on already grouped data. Both can perform weighted counts,
#' if you give `wt` the name of a variable to weight by.
#'
#' `count()` and `tally()` are shortcuts for `summarise()`; `add_count()`
#' and `add_tally()` perform the identical operations but use `mutate()`
#' instead of `summarise()` so that they add a new column with group-wise
#' counts.
#'
#' @param x A data frame, data frame extension (e.g. a tibble), or a
#'   lazy data frame (e.g. from dbplyr or dtplyr).
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Variables to group by.
#' @param wt <[`tidy-eval`][dplyr_tidy_eval]> If omitted, will count the number of rows.
#'   If specified, will perform a "weighted" tally by summing the (non-missing)
#'   values of variable `wt`.
#'
#'   If omitted, and column `n` exists, it will automatically be used as a
#'   weighting variable, although you will have to specify `name` to provide
#'   a new name for the output.
#' @param sort If `TRUE` will sort output in descending order of `n`.
#' @param name The name of the new column in the output.
#'
#'   If omitted, it will default to `n`. If there's already a column called `n`,
#'   it will error, and require you to specify the name.
#' @param .drop For `count()`: if `FALSE` will include counts for empty groups
#'   (i.e. for levels of factors that don't exist in the data). Deprecated for
#'   `add_count()` since it didn't actually affect the output.
#' @return
#' An object of the same type as `.data`. `count()` and `add_count()`
#' group transiently, so the output has the same groups as the input.
#' @export
#' @examples
#' # count() is a convenient way to get a sense of the distribution of
#' # values in a dataset
#' starwars %>% count(species)
#' starwars %>% count(species, sort = TRUE)
#'
#' # use the `wt` argument to perform a weighted count. This is useful
#' # when the data has already been aggregated once
#' df <- tribble(
#'   ~name,    ~gender,   ~runs,
#'   "Max",    "male",       10,
#'   "Sandra", "female",      1,
#'   "Susan",  "female",      4
#' )
#' # counts rows:
#' df %>% count(gender)
#' # counts runs:
#' df %>% count(gender, wt = runs)
#'
#' # tally() is a lower-level function that assumes you've done the grouping
#' starwars %>% tally()
#' starwars %>% group_by(species) %>% tally()
#'
#' # both count() and tally() have add_ variants that work like
#' # mutate() instead of summarise
#' df %>% add_count(gender, wt = runs)
#' df %>% add_tally(wt = runs)
tally <- function(x, wt = NULL, sort = FALSE, name = NULL) {
  n <- tally_n(x, {{ wt }})
  name <- check_name(x, name)
  out <- summarise(x, !!name := !!n)

  if (sort) {
    arrange(out, desc(!!sym(name)))
  } else {
    out
  }
}

#' @rdname tally
#' @export
add_tally <- function(x, wt = NULL, sort = FALSE, name = NULL) {
  n <- tally_n(x, {{ wt }})
  name <- check_name(x, name)
  out <- mutate(x, !!name := !!n)

  if (sort) {
    arrange(out, desc(!!sym(name)))
  } else {
    out
  }
}

#' @export
#' @rdname tally
count <- function(x, ..., wt = NULL, sort = FALSE, name = NULL, .drop = group_by_drop_default(x)) {

  if (!missing(...)) {
    out <- group_by(x, ..., .add = TRUE, .drop = .drop)
  } else {
    out <- x
  }

  out <- tally(out, wt = !!enquo(wt), sort = sort, name = name)
  dplyr_reconstruct(out, x)
}

#' @rdname tally
#' @export
add_count <- function(x, ..., wt = NULL, sort = FALSE, name = NULL, .drop = deprecated()) {
  if (!missing(.drop)) {
    lifecycle::deprecate_warn("1.0.0", "add_count(.drop = )")
  }

  if (!missing(...)) {
    out <- group_by(x, ..., .add = TRUE)
  } else {
    out <- x
  }
  out <- add_tally(out, wt = !!enquo(wt), sort = sort, name = name)

  name <- check_name(x, name)
  x[[name]] <- out[[name]]
  x
}

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
#' \Sexpr[stage=render,results=Rd]{dplyr:::methods_rd("rename")}.
#' @family single table verbs
#' @export
#' @examples
#' iris <- as_tibble(iris) # so it prints a little nicer
#' rename(iris, petal_length = Petal.Length)
#' @export
rename <- function(.data, ...) {
  UseMethod("rename")
}

#' @export
rename.data.frame <- function(.data, ...) {
  loc <- tidyselect::eval_rename(expr(c(...)), .data)
  # eval_rename() only returns changes
  names <- names(.data)
  names[loc] <- names(loc)

  set_names(.data, names)
}
#' Group input by rows
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("questioning")}
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
#' @param data Input data frame.
#' @export
#' @examples
#' df <- expand.grid(x = 1:3, y = 3:1)
#' df_done <- df %>% rowwise() %>% do(i = seq(.$x, .$y))
#' df_done
#' df_done %>% summarise(n = length(i))
rowwise <- function(data) {
  abort_if_not(is.data.frame(data))

  structure(data, class = c("rowwise_df", "tbl_df", "tbl", "data.frame"))
}

# Do ---------------------------------------------------------------------------

#' @export
do.rowwise_df <- function(.data, ...) {
  # Create ungroup version of data frame suitable for subsetting
  group_data <- ungroup(.data)

  args <- enquos(...)
  named <- named_args(args)

  # Create new environment, inheriting from parent, with an active binding
  # for . that resolves to the current subset. `_i` is found in environment
  # of this function because of usual scoping rules.
  mask <- new_data_mask(new_environment())
  current_row <- function() lapply(group_data[`_i`, , drop = FALSE], "[[", 1)
  env_bind_do_pronouns(mask, current_row)

  n <- nrow(.data)
  m <- length(args)

  out <- replicate(m, vector("list", n), simplify = FALSE)
  names(out) <- names(args)
  p <- progress_estimated(n * m, min_time = 2)

  for (`_i` in seq_len(n)) {
    for (j in seq_len(m)) {
      out[[j]][`_i`] <- list(eval_tidy(args[[j]], mask))
      p$tick()$print()
    }
  }

  if (!named) {
    label_output_dataframe(NULL, out, groups(.data), group_by_drop_default(.data))
  } else {
    label_output_list(NULL, out, groups(.data))
  }
}
#' Sample n rows from a table
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("retired")}
#' `sample_n()` and `sample_frac()` have been retired in favour of
#' [slice_sample()]. While they will not be deprecated in the near future,
#' retirement means that we will only perform critical bug fixes, so we recommend
#' moving to the newer alternative.
#'
#' These functions were retired because we realised it was more convenient to
#' have two mutually exclusive arguments to one function, rather than two
#' separate functions. This also made it to clean up a few other smaller
#' design issues with `sample_n()`/`sample_frac`:
#'
#' * The connection to `slice()` was not obvious.
#' * The name of the first argument, `tbl`, is inconsistent with other
#'   single table verbs which use `.data`.
#' * The `size` argument uses tidy evaluation, which is surprising and
#'   undocumented.
#' * It was easier to remove the deprecated `.env` argument.
#' * `...` was in a suboptimal position.
#'
#' @keywords internal
#' @param tbl A data.frame.
#' @param size <[`tidy-select`][dplyr_tidy_select]>
#'   For `sample_n()`, the number of rows to select.
#'   For `sample_frac()`, the fraction of rows to select.
#'   If `tbl` is grouped, `size` applies to each group.
#' @param replace Sample with or without replacement?
#' @param weight <[`tidy-select`][dplyr_tidy_select]> Sampling weights.
#'   This must evaluate to a vector of non-negative numbers the same length as
#'   the input. Weights are automatically standardised to sum to 1.
#' @param .env DEPRECATED.
#' @param ... ignored
#' @examples
#' by_cyl <- mtcars %>% group_by(cyl)
#'
#' # sample_n() -> slice_sample() ----------------------------------------------
#' sample_n(mtcars, 10)
#' sample_n(mtcars, 50, replace = TRUE)
#' sample_n(mtcars, 10, weight = mpg)
#'
#' # Changes:
#' # * explicitly name the `n` argument,
#' # * the `weight` argument is now `weight_by`.
#'
#' slice_sample(mtcars, n = 10)
#' slice_sample(mtcars, n = 50, replace = TRUE)
#' slice_sample(mtcars, n = 10, weight_by = mpg)
#'
#' # Note that sample_n() would error if n was bigger than the group size
#' # slice_sample() will just use the available rows for consistency with
#' # the other slice helpers like slice_head()
#'
#' # sample_frac() -> slice_sample() -------------------------------------------
#' sample_frac(mtcars)
#' sample_frac(mtcars, replace = TRUE)
#'
#' # Changes:
#' # * use prop = 1 to randomly sample all rows
#'
#' slice_sample(mtcars, prop = 1)
#' slice_sample(mtcars, prop = 1, replace = TRUE)
#'
#' @export
sample_n <- function(tbl, size, replace = FALSE, weight = NULL, .env = NULL, ...) {
  UseMethod("sample_n")
}

#' @export
sample_n.default <- function(tbl, size, replace = FALSE, weight = NULL,
                             .env = parent.frame(), ...) {
  bad_args("tbl", "must be a data frame, not {friendly_type_of(tbl)}")
}

#' @export
sample_n.data.frame <- function(tbl, size, replace = FALSE,
                                weight = NULL, .env = NULL, ...) {
  if (!is_null(.env)) {
    inform("`.env` is deprecated and no longer has any effect")
  }

  size <- enquo(size)
  weight <- enquo(weight)

  slice(tbl, sample.int(n(), check_size(!!size, n(), replace = replace), replace = replace, prob = !!weight))
}

#' @rdname sample_n
#' @export
sample_frac <- function(tbl, size = 1, replace = FALSE, weight = NULL, .env = NULL, ...) {
  UseMethod("sample_frac")
}

#' @export
sample_frac.default <- function(tbl, size = 1, replace = FALSE, weight = NULL,
                                .env = parent.frame(), ...) {
  bad_args("tbl", "must be a data frame, not {friendly_type_of(tbl)}")
}

#' @export
sample_frac.data.frame <- function(tbl, size = 1, replace = FALSE,
                                   weight = NULL, .env = NULL, ...) {
  if (!is_null(.env)) {
    inform("`.env` is deprecated and no longer has any effect")
  }

  size <- enquo(size)
  weight <- enquo(weight)

  slice(tbl, sample.int(n(), round(n() * check_frac(!!size, replace = replace)), replace = replace, prob = !!weight))
}


