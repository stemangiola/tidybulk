

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
  df <- tbl
  # Through warning if data frame is not numerical beside the rownames column (if present)
  check_df <- df
  if (!quo_is_null(rownames)) {
    check_df <- check_df[,-1]
  }
  if (do_check && check_df %>% dplyr::summarise_all(class) %>% tidyr::gather(variable, class) %>% pull(class) %>% unique() %>% `%in%`(c("numeric", "integer")) %>% not() %>% any()) {
    warning("tidybulk says: there are NON-numerical columns, the matrix will NOT be numerical")
  }
  df <- as.data.frame(df)
  # Deal with rownames column if present
  if (!quo_is_null(rownames)) {
    df <- df %>% magrittr::set_rownames(tbl %>% pull(!!rownames)) %>% select(-1)
  }
  # Convert to matrix
  as.matrix(df)
}