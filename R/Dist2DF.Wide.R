#' 'dist' to 'data.frame' wide format conversion
#'
#' A function to change the data class from dist to data.frame in wide format.
#'
#' @param inDist A matrix of class `dist`, as the output of `ape::dist.dna`.
#' @param Convert100 If `TRUE` converts proportion distance (or divergence) input to percentages (`FASLE` is the default).
#'
#' @returns The distance data in class 'data.frame' wide format.
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
Dist2DF.Wide <- function(inDist,
                         Convert100 = FALSE) {
  if (class(inDist) != "dist") stop("Wrong input format for distance data -- class() of input is not 'dist'")
  if (Convert100 != TRUE & Convert100 != FALSE) stop("Wrong input for 'Convert100' -- Options are 'TRUE', 'FALSE'")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  if (Convert100 == FALSE){
    inDistLong <- data.frame(
      row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
      col = rep(B[-length(B)], (length(B)-1):1),
      Distance = as.vector(inDist))
    outW <- reshape(inDistLong,  ## Change the long format table into a classic wide table
                    idvar = "row",
                    timevar = "col",
                    direction = "wide")
    colnames(outW) <- gsub("Distance.","",colnames(outW))
    names(outW)[1] <- NA
    outW

  } else
    if (Convert100 == TRUE){
      inDistLong <- data.frame(
        row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
        col = rep(B[-length(B)], (length(B)-1):1),
        Distance = as.vector(inDist * 100))
      outW <- reshape(inDistLong,  ## Change the long format table into a classic wide table
                      idvar = "row",
                      timevar = "col",
                      direction = "wide")
      colnames(outW) <- gsub("Distance.","",colnames(outW))
      names(outW)[1] <- NA
      outW
    }
}
