

#' Get matrix from tibble
#'
#'
#'
#'
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames The column name of the input tibble that will become the rownames of the output matrix
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#' library(tibble)
#' tibble(.feature = "CD3G", count=1) |> as_matrix(rownames=.feature)
#'
#' @export
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  
  # Fix NOTEs
  . = NULL
  
  rownames = enquo(rownames)
  tbl %>%
    
    # Through warning if data frame is not numerical beside the rownames column (if present)
    ifelse_pipe(
      do_check &&
        tbl %>%
        # If rownames defined eliminate it from the data frame
        ifelse_pipe(!quo_is_null(rownames), ~ .x[,-1], ~ .x) |>
        dplyr::summarise_all(class) |>
        tidyr::gather(variable, class) |>
        pull(class) |>
        unique() %>%
        `%in%`(c("numeric", "integer")) |> not() |> any(),
      ~ {
        warning("tidybulk says: there are NON-numerical columns, the matrix will NOT be numerical")
        .x
      }
    ) |>
    as.data.frame() |>
    
    # Deal with rownames column if present
    ifelse_pipe(
      !quo_is_null(rownames),
      ~ .x |>
        magrittr::set_rownames(tbl |> pull(!!rownames)) |>
        select(-1)
    ) |>
    
    # Convert to matrix
    as.matrix()
}