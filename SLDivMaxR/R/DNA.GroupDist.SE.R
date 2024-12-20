#' Divergence Standard Error
#'
#' Calculates standard error of specified between and within group genetic divergences.
#'
#'
#' @param inDNA the imput aligned sequences in `.fasta` format.
#' @param Dist.type specifies if the distance output should be in proportions (`"prop"`), percentages (`"percentage"`), or number of nucleotides that are different (`"nucleotide"`).
#' @param Model specifies the evolutionary model to use. Same as in `ape::dist.dna`. Options are: `"RAW"`, `"JC69"`, `"K80"` (default), `"F81"`, `"K81"`, `"F84"`, `"T92"`, `"TN93"`, `"GG95"`, `"LOGDET"`, `"BH87"`, `"PARALIN"`, `"N"`, `"TS"`, `"TV"`, `"INDEL"`, and `"INDELBLOCK"`. Note: if `Model` is `"N"` and `Dist.type` is set for `"prop"` or `"percentage"`, then an error will be return.
#' @param GAMMA same as `gamma` in `ape::dist.dna()`.
#' @param PW.deletion specifies if Pairwise deletion should be considered (`TURE`, default) or not (`FALSE`).
#' @param by.group a vector to group sequences on different categories (e.g., populations, species, lineages, etc.). The vector should follow the same order as the sequence alignment. The vector can be the `.csv` output of species delimitation analyses such as `ASAP` and `ABGD`.
#' @param within.group creates a `data.frame` with each within group standard error, and number of saples ("N") by group. A vector with discreet groupings should be imported in `by.group`. Groups with N = 1 will return `NA`'s.
#'
#' @returns A data.frame with the standard errors of the mean distance between or within groups distances.
#'
#' @export
DNA.GroupDist.SE <- function(inDNA,
                             Dist.type = "percent",
                             Out.Format = "long",
                             Model = "K80",
                             GAMMA = FALSE,
                             PW.deletion = TRUE,
                             by.group,
                             within.group = FALSE) {

  MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", "TN93", "GG95", "LOGDET", "BH87", "PARALIN",
              "N", "TS", "TV", "INDEL", "INDELBLOCK")
  imod <- pmatch(toupper(Model), MODELS)
  if (is.na(imod)) stop(paste("'Model' must be one of:", paste("\"", MODELS,"\"", sep = "", collapse = " "),
                              ". See ape::dist.dna for more information.", sep = ""))

  DT <- c("percent", "prop", "nucleotide")
  iDT <- pmatch((Dist.type), DT)
  if (is.na(iDT)) stop("Wrong 'Dist.type' input. Options are: 'prop', 'percent', 'nucleotide'.")

  if (is.null(by.group) == TRUE) stop("Missing 'by.group' argument")
  if (length(by.group) != length(inDNA)) stop("'by.group' and 'inDNA' are not of the same length")

  if (is.null(by.group) == FALSE ) {
    if (Dist.type == "prop" & Model != "N" & Model != "RAW"){
      inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = FALSE,
                              as.matrix = FALSE, gamma = GAMMA)} else
                                if (Dist.type == "percent" & Model != "N" & Model != "RAW") {
                                  inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = FALSE,
                                                          as.matrix = FALSE, gamma = GAMMA) * 100
                                } else
                                  if (Dist.type == "nucleotide"){
                                    inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = ("N"), variance = FALSE,
                                                            as.matrix = FALSE, gamma = GAMMA)}
  }

  G <- factor(by.group, exclude = NULL)

  grow <- G[as.dist(row(as.matrix(inDist)))]
  gcol <- G[as.dist(col(as.matrix(inDist)))]

  fx <- as.numeric(grow) >= as.numeric(gcol)
  x1 <- ifelse(fx, grow, gcol)
  x2 <- ifelse(!fx, grow, gcol)
  n <- table(G)
  tx <- matrix(TRUE, nlevels(G), nlevels(G))
  diag(tx) <- n > 1
  tx[upper.tri(tx)] <- FALSE
  out <- matrix(NA, nlevels(G), nlevels(G))
  tmp <- tapply(inDist, list(x1, x2), SE)
  out[(tx)] <- tmp[!is.na(tmp)]
  rownames(out) <- colnames(out) <- levels(G)
  class(out) <- c("matix")

  if (within.group == TRUE) {
    data.frame(Group = names(n), N = as.vector(n), Within.Group.SE = as.vector(diag(out)))
  } else

    if (Out.Format == "long" & within.group == FALSE) {
      Dout <- as.dist(out)
      x <- Dist2DF.Long(Dout)
      names(x) <- c("x1", "x2", "Standard_Error")
      x
    } else

      if (Out.Format == "wide" & within.group == FALSE) {
        Dout <- as.dist(out)
        Dist2DF.Wide(Dout)}
}
