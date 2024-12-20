#' Mean Distances by Group
#'
#' A function that applies `hclust` and `cutree` to a matrix of class `"dist"` to cluster samples by a  maximum distance threshold or to a set number of groups.
#'
#' @param inDust A `"dist"` object matrix, output of `ape::dist.dna` or alike matrix
#' @param Convert100 Converts proportions to percentages.
#' @param get.group Provides discreet grouping using maximum divergence threshold between sequences.
#' @param n.groups Groups low-diverging sequences to a specified number of groups, regardless their within group divergences.
#'
#' @returns A data.frame with the mean distance between groups group distances from a distance matrix produced by `ape::dist.dna`
#'      Note: `Dist2DF.Groups` should not be used as a species delimitation analysis as it does not explicitly tests a hypothesis.
#'
#' @export
Dist2DF.Groups <- function(inDist,
                           Convert100 = FALSE,
                           get.group = NULL,
                           n.groups = NULL) {
  if (class(inDist) != "dist") stop("Wrong input format for distance data -- class() of input is not 'dist'")
  if (Convert100 != TRUE & Convert100 != FALSE) stop("Wrong input for 'Convert100' -- Options are 'TRUE', 'FALSE'")
  if (is.null(get.group) == FALSE & class(get.group) != "numeric") stop("Wrong input for 'get.group' -- must be numeric")
  if (is.null(n.groups) == FALSE & class(n.groups) != "numeric") stop("Wrong input for 'n.groups' -- must be numeric")
  if (is.null(get.group) == TRUE & is.null(n.groups) == TRUE) stop("No imput for neither `get.group` or `n.groups`")
  if (is.null(get.group) == FALSE & is.null(n.groups) == FALSE)
    stop("Imput entered for both `get.group` or `n.groups`, choose one.")
  if (is.null(get.group) == FALSE | is.null(n.groups) == FALSE) {
    x1 <- hclust(inDist, "complete")
    x2 <- cutree(x1, k = n.groups, h = get.group)
    data.frame(Specimen = names(x2),
               Group = paste("group", sprintf(paste0("%0", nchar(max(x2)),".0f"),as.numeric(x2)), sep = "_"))
  }
  if (is.null(get.group) == FALSE | is.null(n.groups) == FALSE & Convert100 == TRUE) {
    x1 <- hclust((inDist*100), "complete")
    x2 <- cutree(x1, k = n.groups, h = get.group)
    data.frame(Specimen = names(x2),
               Group = paste("group", sprintf(paste0("%0", nchar(max(x2)),".0f"),as.numeric(x2)), sep = "_"))
  }
}
