
#' Log10 reverse transformation for ggplot2
#'
#' Creates a transformation that applies -log10(x) to data, useful for
#' visualizing p-values or other values where smaller values should be
#' displayed larger.
#'
#' @return A transformation object that can be used with ggplot2's scale functions
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' # Example usage with p-values
#' ggplot(data, aes(x = pvalue)) +
#'   geom_histogram() +
#'   scale_x_continuous(trans = log10_reverse_trans())
#' }
#'
#' @importFrom scales trans_new
#' @importFrom scales label_scientific
#' @export
log10_reverse_trans <- function() {
  trans <- function(x) -log10(x)
  inv <- function(x) 10^(-x)
  
  # Custom breaks function for even spacing in transformed space
  breaks <- function(x) {
    # Get the range in transformed space
    trans_range <- trans(x)
    # Remove infinite values
    trans_range <- trans_range[is.finite(trans_range)]
    if (length(trans_range) == 0) return(numeric(0))
    
    # Create evenly spaced breaks in transformed space
    min_trans <- min(trans_range)
    max_trans <- max(trans_range)
    even_breaks <- seq(min_trans, max_trans, length.out = 5)
    # Convert back to original space for labeling
    result <- inv(even_breaks)
    # Force the result to be returned
    return(result)
  }
  
  # Custom format function for scientific notation
  format <- function(x) {
    scales::label_scientific(digits = 2)(x)
  }

  scales::trans_new("log10_reverse", trans, inv, breaks, format = format)
}

#' scale_y_log10_reverse
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description A wrapper function that provides evenly spaced ticks with scientific notation for log10 reverse transformed y-axis. This is particularly useful for volcano plots and other plots showing p-values. The function applies a log10 transformation and reverses the axis to better visualize p-values without hard transforming the data, while maintaining the original p-value scale for interpretation. This allows you to see the full range of p-values with proper scaling while keeping the original values readable.
#'
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom scales label_scientific
#'
#' @param breaks Number of breaks to display (default: 5)
#' @param digits Number of digits for scientific notation (default: 2)
#' @param ... Additional arguments passed to scale_y_continuous
#'
#' @return A ggplot2 scale object
#'
#' @examples
#'
#' library(ggplot2)
#' library(tibble)
#'
#' # Create test data
#' test_data <- tibble(
#'   pvalue = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5),
#'   fold_change = 1:6
#' )
#'
#' # Use the wrapper function
#' test_data |>
#'   ggplot(aes(fold_change, pvalue)) +
#'   geom_point() +
#'   scale_y_log10_reverse()
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. https://ggplot2.tidyverse.org
#'
#' @export
scale_y_log10_reverse <- function(breaks = 5, digits = 2, ...) {
  
  # Function to create evenly spaced breaks for log10_reverse transformation
  create_even_breaks <- function(pvalues, n = breaks) {
    # Remove zeros and get finite values
    pvalues <- pvalues[pvalues > 0 & is.finite(pvalues)]
    
    if (length(pvalues) == 0) return(numeric(0))
    
    # Transform to log10 space
    trans_values <- -log10(pvalues)
    
    # Create evenly spaced breaks in transformed space
    min_trans <- min(trans_values)
    max_trans <- max(trans_values)
    even_breaks <- seq(min_trans, max_trans, length.out = n)
    
    # Convert back to original space
    10^(-even_breaks)
  }
  
  # Create a custom breaks function that will be called by ggplot2
  breaks_func <- function(x) {
    create_even_breaks(x, n = breaks)
  }
  
  # Return the scale with custom breaks and scientific notation
  scale_y_continuous(
    trans = "log10_reverse",
    breaks = breaks_func,
    labels = label_scientific(digits = digits),
    ...
  )
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
