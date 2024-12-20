#' Standard Error
#'
#' Calculates standard error on a vector
#'
#' @param x a numeric vector. `NA` will be automatically ignored.
#'
#' @returns Standard error calculation
#'
#' @examples
#' a <- c(1:10)
#' SE(a)
#'
#' a <- c(1:10, NA, NA)
#' SE(a)
#'
#' a <- c(NA, NA, NA)
#' SE(a)
#'
#' @export
SE <- function(x) {
  x <- as.vector(na.omit(x))
  sd(x)/sqrt(length(x))
}
