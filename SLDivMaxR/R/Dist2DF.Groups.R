#' Mean Distances between Group
#'
#' A function that clusters samples of a distance matrix (`class == 'dist'`) by a maximum distance threshold or to a set number of groups.
#'
#' @param inDust A `"dist"` object matrix, output of `ape::dist.dna` or alike matrix
#' @param Convert100 Converts proportions to percentages.
#' @param get.group Provides discreet grouping using maximum divergence threshold between sequences.
#' @param n.groups Groups low-diverging sequences to a specified number of groups, regardless their within group divergences.
#' @param method Options to use pairwise comparisons ("PW", default), or using the `hclust` and `cutree` approach ("CT"). If using `n.groups`, then `method` should be set to "CT"
#' @param G.Name A user-defined character to name the groupings. Default is `"group"``.
#'
#' @returns A data.frame with the mean distance between groups group distances from a distance matrix produced by `ape::dist.dna`
#'      Note: `Dist2DF.Groups` should not be used as a species delimitation analysis as it does not explicitly tests a hypothesis.
#'
#' @export
Dist2DF.Groups <- function(inDist,
                           Convert100 = FALSE,
                           get.group = NULL,
                           n.groups = NULL,
                           method = "PW",
                           G.Name = "group") {
  if (class(inDist) != "dist") stop("Wrong input format for distance data -- class() of input is not 'dist'")
  if (Convert100 != TRUE & Convert100 != FALSE) stop("Wrong input for 'Convert100' -- Options are 'TRUE', 'FALSE'")
  if (is.null(get.group) == FALSE & class(get.group) != "numeric") stop("Wrong input for 'get.group' -- must be numeric")
  if (is.null(n.groups) == FALSE & class(n.groups) != "numeric") stop("Wrong input for 'n.groups' -- must be numeric")
  if (is.null(get.group) == TRUE & is.null(n.groups) == TRUE) stop("No imput for neither `get.group` or `n.groups`")
  if (is.null(get.group) == FALSE & is.null(n.groups) == FALSE)
    stop("Imput entered for both `get.group` or `n.groups`, choose one.")
  if (method != "CT" & method != "PW")
    stop("Invalid imput entered for `method`. Options are `CT` or `PW`.")
  if (method == "PW" & is.null(n.groups) == FALSE)
    stop("Invalid imput entered for `method`. If `n.groups == TRUE`, option for `method` should be `CT`.")

  if (Convert100) {
    inDist <- inDist * 100
  }
  if (method == "CT"){
  if (is.null(get.group) == FALSE | is.null(n.groups) == FALSE) {
    x1 <- hclust(inDist, "complete")
    x2 <- cutree(x1, k = n.groups, h = get.group)
    data.frame(Specimen = names(x2),
               Group = paste(G.Name, sprintf(paste0("%0", nchar(max(x2)),".0f"),as.numeric(x2)), sep = "_"))
  }
    } else
    if (method == "PW"){

      distMatrix <- as.matrix(inDist)
      diag(distMatrix) <- NA

      GPs <- rep(NA, length(rownames(distMatrix)))
      gID <- 0
      for (i in 1:(nrow(distMatrix) - 1)) {
        if (is.na(GPs[i])) {
          gID  <- gID  + 1
          GPs[i] <- gID
          for (j in (i + 1):nrow(distMatrix)) {
            if (distMatrix[i, j] <= get.group) {
              GPs[j] <- gID
            }
          }
        }
      }
      maxGP <- max(GPs, na.rm = TRUE)
      data.frame(Specimen = rownames(distMatrix),
                 Group = paste(G.Name, sprintf(paste0("%0", nchar(max(maxGP)),".0f"), as.numeric(GPs)), sep = "_"))
    }
}
