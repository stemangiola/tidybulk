#' log10_reverse_trans
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description it perform log scaling and reverse the axis. Useful to plot negative log probabilities. To not be used directly but with ggplot (e.g. scale_y_continuous(trans = "log10_reverse") )
#'
#' @importFrom scales trans_new
#' @importFrom scales log_breaks
#'
#' @return A scales object
#'
#' @examples
#'
#' library(ggplot2)
#' library(tibble)
#'
#' tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) |>
#'  ggplot(aes(fold_change , pvalue)) +
#'  geom_point() +
#'  scale_y_continuous(trans = "log10_reverse")
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. https://ggplot2.tidyverse.org
#'
#' @export
log10_reverse_trans <- function() {
  trans <- function(x) -log10(x)
  inv <- function(x) 10^(-x)

  trans_new("log10_reverse", trans, inv, log_breaks(base = 10))
}

#' logit scale
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description it perform logit scaling with right axis formatting. To not be used directly but with ggplot (e.g. scale_y_continuous(trans = "log10_reverse") )
#'
#' @importFrom scales label_scientific
#' @importFrom scales extended_breaks
#' @importFrom stats qlogis plogis
#'
#' @return A scales object
#'
#' @examples
#'
#' library(ggplot2)
#' library(tibble)
#'
#' tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) |>
#'  ggplot(aes(fold_change , pvalue)) +
#'  geom_point() +
#'  scale_y_continuous(trans = "log10_reverse")
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. https://ggplot2.tidyverse.org
#'
#' @export
logit_trans <- function(){

  check_and_install_packages("functional")

  trans <- qlogis
  inv <- plogis

  trans_new("logit",
            transform = trans,
            inverse = inv,
            breaks = functional::Compose(trans, extended_breaks(), inv),
            format = label_scientific(digits = 2)
  )
}
