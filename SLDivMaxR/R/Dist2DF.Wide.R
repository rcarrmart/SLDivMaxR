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
#' library("ape")
#' COI.data <- read.FASTA(file="COI.fasta")
#' Dist <- dist.dna(COI.data, model = "K80")
#'
#' # Wide matrix-like dataframe of proportional distances
#' DistW1 <- Dist2DF.Wide(Dist)
#' head(DistW1)
#'
#' # Wide matrix-like dataframe of distances in percentages
#' DistW2 <- Dist2DF.Wide(Dist, Convert100 = TRUE)
#' head(DistW2)
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
