
#' Arrange rows by column values
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
arrange.data.frame <- function(.data, ..., .by_group = FALSE) {
  if (missing(...)) {
    return(.data)
  }

  loc <- arrange_rows(.data, ..., .by_group = .by_group)
  dplyr_row_slice(.data, loc)
}


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
bind_rows <- function(..., .id = NULL) {
  dots <- dots_values(...)
  if (length(dots) == 1 && is.list(dots[[1]]) && !is.data.frame(dots[[1]])) {
    dots <- dots[[1]]
  }
  dataframe_ish <- function(.x) {
    is.data.frame(.x) || (vec_is(.x) && !is.null(names(.x)))
  }
  dots <- keep(
    flatten_if(dots, function(.x) is.list(.x) && !is.data.frame(.x)),
    function(.x) !is.null(.x)
  )

  dots <- keep(dots, function(.x) !is.null(.x))
  dots <- flatten_if(dots, function(.x) is.list(.x) && !dataframe_ish(.x))

  if (!is_null(.id)) {
    if (!(is_string(.id))) {
      bad_args(".id", "must be a scalar string, ",
        "not {friendly_type_of(.id)} of length {length(.id)}"
      )
    }
    if (!all(have_name(dots) | map_lgl(dots, is_empty))) {
      dots <- compact(dots)
      names(dots) <- seq_along(dots)
    }
  }
  if (!is.null(names(dots)) && !all(map_lgl(dots, dataframe_ish))) {
    dots <- list(as_tibble(dots))
  }

  for (i in seq_along(dots)) {
    .x <- dots[[i]]
    if (!is.data.frame(.x) && !vec_is(.x)) {
      abort(glue("Argument {i} must be a data frame or a named atomic vector"))
    }

    if (is.null(names(.x))) {
      abort(glue("Argument {i} must have names"))
    }
  }

  dots <- map(dots, function(.x) if(is.data.frame(.x)) .x else tibble(!!!as.list(.x)))
  result <- vec_rbind(!!!dots, .names_to = .id)
  if (length(dots) && is_tibble(first <- dots[[1L]])) {
    if (is_grouped_df(first)) {
      result <- grouped_df(result, group_vars(first), group_by_drop_default(first))
    } else {
      class(result) <- class(first)
    }
  }
  result
}

#' @export
#' @rdname bind
bind_cols <- function(...) {
  dots <- dots_values(...)
  not_null <- function(.x) !is.null(.x)
  dots <- keep(dots, not_null)

  # nothing to bind, return a dummy tibble
  if (!length(dots)) {
    return(tibble())
  }

  # Before things are squashed, we need
  # some information about the "first" data frame
  if (is.data.frame(dots[[1]]) || !is.list(dots[[1]])) {
    first <- dots[[1]]
  } else {
    first <- dots[[1]][[1]]
  }

  dots <- squash_if(dots, function(.x) is.list(.x) && !is.data.frame(.x))
  dots <- keep(dots, not_null)
  if (!length(dots)) {
    return(tibble())
  }

  res <- vec_cbind(!!!dots)
  if (length(dots)) {
    if (is_grouped_df(first)) {
      res <- grouped_df(res, group_vars(first), group_by_drop_default(first))
    } else if(inherits(first, "rowwise_df")){
      res <- rowwise(res)
    } else if(is_tibble(first) || !is.data.frame(first)) {
      res <- as_tibble(res)
    }
  }
  res
}

#' Arrange rows by a selection of variables
#'
#' These [scoped] variants of [arrange()] sort a data frame by a
#' selection of variables. Like [arrange()], you can modify the
#' variables before ordering with the `.funs` argument.
#'
#' @inheritParams scoped
#' @inheritParams arrange
#'
#' @section Grouping variables:
#'
#' The grouping variables that are part of the selection participate
#' in the sorting of the data frame.
#'
#' @export
#' @examples
#' df <- as_tibble(mtcars)
#' df
#' arrange_all(df)
#'
#' # You can supply a function that will be applied before taking the
#' # ordering of the variables. The variables of the sorted tibble
#' # keep their original values.
#' arrange_all(df, desc)
#' arrange_all(df, list(~desc(.)))
arrange_all <- function(.tbl, .funs = list(), ..., .by_group = FALSE) {
  funs <- manip_all(.tbl, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  if (!length(funs)) {
    funs <- syms(tbl_vars(.tbl))
  }
  arrange(.tbl, !!!funs, .by_group = .by_group)
}
#' @rdname arrange_all
#' @export
arrange_at <- function(.tbl, .vars, .funs = list(), ..., .by_group = FALSE) {
  funs <- manip_at(.tbl, .vars, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  if (!length(funs)) {
    funs <- tbl_at_syms(.tbl, .vars, .include_group_vars = TRUE)
  }
  arrange(.tbl, !!!funs, .by_group = .by_group)
}
#' @rdname arrange_all
#' @export
arrange_if <- function(.tbl, .predicate, .funs = list(), ..., .by_group = FALSE) {
  funs <- manip_if(.tbl, .predicate, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  if (!length(funs)) {
    funs <- tbl_if_syms(.tbl, .predicate, .include_group_vars = TRUE)
  }
  arrange(.tbl, !!!funs, .by_group = .by_group)
}
#' Select distinct rows by a selection of variables
#'
#' These [scoped] variants of [distinct()] extract distinct rows by a
#' selection of variables. Like `distinct()`, you can modify the
#' variables before ordering with the `.funs` argument.
#'
#' @param .keep_all If `TRUE`, keep all variables in `.data`.
#'   If a combination of `...` is not distinct, this keeps the
#'   first row of values.
#' @inheritParams scoped
#' @export
#'
#' @section Grouping variables:
#'
#' The grouping variables that are part of the selection are taken
#' into account to determine distinct rows.
#'
#' @examples
#' df <- tibble(x = rep(2:5, each = 2) / 2, y = rep(2:3, each = 4) / 2)
#' df
#' distinct_all(df)
#' distinct_at(df, vars(x,y))
#' distinct_if(df, is.numeric)
#'
#' # You can supply a function that will be applied before extracting the distinct values
#' # The variables of the sorted tibble keep their original values.
#' distinct_all(df, round)
#' arrange_all(df, list(~round(.)))
distinct_all <- function(.tbl, .funs = list(), ..., .keep_all = FALSE) {
  funs <- manip_all(.tbl, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  if (!length(funs)) {
    funs <- syms(tbl_vars(.tbl))
  }
  distinct(.tbl, !!!funs, .keep_all = .keep_all)
}
#' @rdname distinct_all
#' @export
distinct_at <- function(.tbl, .vars, .funs = list(), ..., .keep_all = FALSE) {
  funs <- manip_at(.tbl, .vars, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  if (!length(funs)) {
    funs <- tbl_at_syms(.tbl, .vars, .include_group_vars = TRUE)
  }
  distinct(.tbl, !!!funs, .keep_all = .keep_all)
}
#' @rdname distinct_all
#' @export
distinct_if <- function(.tbl, .predicate, .funs = list(), ..., .keep_all = FALSE) {
  funs <- manip_if(.tbl, .predicate, .funs, enquo(.funs), caller_env(), .include_group_vars = TRUE, ...)
  if (!length(funs)) {
    funs <- tbl_if_syms(.tbl, .predicate, .include_group_vars = TRUE)
  }
  distinct(.tbl, !!!funs, .keep_all = .keep_all)
}
#' Filter within a selection of variables
#'
#' These [scoped] filtering verbs apply a predicate expression to a
#' selection of variables. The predicate expression should be quoted
#' with [all_vars()] or [any_vars()] and should mention the pronoun
#' `.` to refer to variables.
#'
#' @inheritParams scoped
#' @param .vars_predicate A quoted predicate expression as returned by
#'   [all_vars()] or [any_vars()].
#'
#'   Can also be a function or purrr-like formula. In this case, the
#'   intersection of the results is taken by default and there's
#'   currently no way to request the union.
#' @param .preserve when `FALSE` (the default), the grouping structure
#'   is recalculated based on the resulting data, otherwise it is kept as is.
#' @export
#'
#' @section Grouping variables:
#'
#' The grouping variables that are part of the selection are taken
#' into account to determine filtered rows.
#'
#' @examples
#' # While filter() accepts expressions with specific variables, the
#' # scoped filter verbs take an expression with the pronoun `.` and
#' # replicate it over all variables. This expression should be quoted
#' # with all_vars() or any_vars():
#' all_vars(is.na(.))
#' any_vars(is.na(.))
#'
#'
#' # You can take the intersection of the replicated expressions:
#' filter_all(mtcars, all_vars(. > 150))
#'
#' # Or the union:
#' filter_all(mtcars, any_vars(. > 150))
#'
#'
#' # You can vary the selection of columns on which to apply the
#' # predicate. filter_at() takes a vars() specification:
#' filter_at(mtcars, vars(starts_with("d")), any_vars((. %% 2) == 0))
#'
#' # And filter_if() selects variables with a predicate function:
#' filter_if(mtcars, ~ all(floor(.) == .), all_vars(. != 0))
#'
#'
#' # We're working on a new syntax to allow functions instead,
#' # including purrr-like lambda functions. This is already
#' # operational, but there's currently no way to specify the union of
#' # the predicate results:
#' mtcars %>% filter_at(vars(hp, vs), ~ . %% 2 == 0)
filter_all <- function(.tbl, .vars_predicate, .preserve = FALSE) {
  syms <- syms(tbl_vars(.tbl))
  pred <- apply_filter_syms(.vars_predicate, syms, .tbl)
  filter(.tbl, !!pred, .preserve = .preserve)
}
#' @rdname filter_all
#' @export
filter_if <- function(.tbl, .predicate, .vars_predicate, .preserve = FALSE) {
  syms <- tbl_if_syms(.tbl, .predicate, .include_group_vars = TRUE)
  pred <- apply_filter_syms(.vars_predicate, syms, .tbl)
  filter(.tbl, !!pred, .preserve = .preserve)
}
#' @rdname filter_all
#' @export
filter_at <- function(.tbl, .vars, .vars_predicate, .preserve = FALSE) {
  syms <- tbl_at_syms(.tbl, .vars, .include_group_vars = TRUE)
  pred <- apply_filter_syms(.vars_predicate, syms, .tbl)
  filter(.tbl, !!pred, .preserve = .preserve)
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


