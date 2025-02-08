#' `dist` to `data.frame` long format conversion
#'
#' A function to change the data class from `dist` to `data.frame` in long format.
#'
#' @param inDist A matrix of class dist, as the output of `ape::dist.dna`.
#' @param Convert100 If TRUE converts proportion distance (or divergence) input to percentages (`FASLE` is the default).
#'
#' @return The distance data in class 'data.frame' long format
#'
#' @examples
#' library(ape)
#' test.dist <- dist.dna(test.DNA)
#'
#' Dist.df <- Dist2DF.Wide(test.dist)
#' Dist.df
#'
#' Dist.df.100 <- Dist2DF.Wide(test.dist, Convert100 = TRUE)
#' Dist.df.100
#'
#' Dist.df.100x <- Dist2DF.Wide(test.dist*100) # equivalent to `Convert100 = TRUE
#' Dist.df.100x
#'
#' @export
Dist2DF.Long <- function(inDist,
                         Convert100 = FALSE) {
  if (class(inDist) != "dist") stop("Wrong input format for distance data -- class() of input is not 'dist'")
  if (Convert100 != TRUE & Convert100 != FALSE) stop("Wrong input for 'Convert100' -- Options are 'TRUE', 'FALSE'")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  if (isTRUE(length(attr(inDist, "variance"))>0) == FALSE & Convert100 == FALSE)
    data.frame(
      x1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
      x2 = rep(B[-length(B)], (length(B)-1):1),
      Distance = as.vector(inDist))
  else
    if (isTRUE(length(attr(inDist, "variance"))>0 ) == TRUE & Convert100 == FALSE){
      data.frame(
        x1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
        x2 = rep(B[-length(B)], (length(B)-1):1),
        Distance = as.vector(inDist),
        Var = attr(inDist, "variance"))
    } else
      if (isTRUE(length(attr(inDist, "variance"))>0) == FALSE  & Convert100 == TRUE){
        data.frame(
          x1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
          x2 = rep(B[-length(B)], (length(B)-1):1),
          Distance = as.vector(inDist * 100))
      }  else
        if (isTRUE(length(attr(inDist, "variance"))>0)== TRUE & Convert100 == TRUE){
          data.frame(
            x1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
            x2 = rep(B[-length(B)], (length(B)-1):1),
            Distance = as.vector(inDist * 100),
            Var = (attr(inDist, "variance")) * 100)
        }
}
