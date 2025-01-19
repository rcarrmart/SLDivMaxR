#' Single locus divergence data.frame
#'
#' A function that applies `ape::dist.dna()` but changes `"dist"` class output from `ape::dist.dna()` (or similar data) to `data.frame`.
#'
#' @param inDNA the input aligned sequences in `.fasta` format.
#'
#' @param Dist.type specifies if the distance output should be in proportions (`"prop"`), percentages (`"percentage"`), or number of nucleotides that are different (`"nucleotide"`).
#'
#' @param Model specifies the evolutionary model to use. Same as in `ape::dist.dna`. Options are: `"RAW"`, `"JC69"`, `"K80"` (default), `"F81"`, `"K81"`, `"F84"`, `"T92"`, `"TN93"`, `"GG95"`, `"LOGDET"`, `"BH87"`, `"PARALIN"`, `"N"`, `"TS"`, `"TV"`, `"INDEL"`, and `"INDELBLOCK"`. Note: if `Model` is `"N"` and `Dist.type` is set for `"prop"` or `"percentage"`, then an error will be return.
#'
#' @param GAMMA same as `gamma` in `ape::dist.dna()`.
#'
#' @param PW.deletion specifies if Pairwise deletion should be considered (`TURE`, default) or not (`FALSE`).
#'
#' @param Var if `TRUE` output is the calculated variation (`FALSE` is default)
#'
#' @param by.group a vector to group sequences on different categories (e.g., populations, species, lineages, etc.). The vector should follow the same order as the sequence alignment, by applying `vegan::meandist()`. The vector can be the `.csv` output of species delimitation analyses such as `ASAP` and `ABGD`.
#'
#' @param `group.summary` applies `vegan::meandist()` to formulate a data.frame with within groups mean distance, between groups mean distance and overall distance. A vector with discreet groupings should be imported in `by.group`.
#'
#' @param `within.group` applies `vegan::meandist()` to create a `data.frame` with each within group mean divergence, and number of samples ("N") by group. A vector with discreet groupings should be imported in `by.group`. Groupls with N = 1 will return `NA`'s.
#'
#' @returns A dataframe with the mean distance between groups or within group distances.
#'          Note: `DNA.Dist.df` should not be used as a species delimitation analysis as it does not explicitly tests a hypothesis.
#'
#' @examples
#' library("ape")
#' #UNder construction!
#'
#' @export
DNA.Dist.df <- function(inDNA,
                        Dist.type = "percent",
                        Out.Format = "long",
                        Model = "K80",
                        GAMMA = FALSE,
                        PW.deletion = TRUE,
                        Var = FALSE,
                        by.group = NULL,
                        group.summary = FALSE,
                        within.group = FALSE) {

  MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", "TN93", "GG95", "LOGDET", "BH87", "PARALIN", "N", "TS",
              "TV", "INDEL", "INDELBLOCK")
  imod <- pmatch(toupper(Model), MODELS)
  if (is.na(imod)) stop(paste("'Model' must be one of:", paste("\"", MODELS,"\"", sep = "", collapse = " "),
                              ". See ape::dist.dna for more information.", sep = ""))

  DT <- c("percent", "prop", "nucleotide")
  iDT <- pmatch((Dist.type), DT)
  if (is.na(iDT)) stop("Wrong 'Dist.type' input. Options are: 'prop', 'percent', 'nucleotide'.")

  if (Dist.type != "nucleotide" & Model == "N")
    stop("Change 'Dist.type' to 'nucleotide' if input for Model is 'N'.")

  dist.fx <- function(model, distType) {
    if (distType == "percent" && model != "N") {
      dist.mtx <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = model, variance = Var,
                                as.matrix = FALSE, gamma = GAMMA) * 100
    } else if (distType == "prop" && model != "N") {
      dist.mtx <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = model, variance = Var,
                                as.matrix = FALSE, gamma = GAMMA)
    } else if (distType == "nucleotide" | model == "N") {
      dist.mtx <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = "N", variance = Var,
                                as.matrix = FALSE, gamma = GAMMA)
    }
    return(dist.mtx)
  }

  out.form <- function(dist.mtx, format) {
    if (format == "long") {
      Dist2DF.Long(dist.mtx)
    } else if (format == "wide") {
      Dist2DF.Wide(dist.mtx)
    }
  }

  if (is.null(by.group)) {
    dist.mtx <- dist.fx(MODELS[imod], Dist.type)
    out.form(dist.mtx, Out.Format)

  } else if (isTRUE(group.summary)) {  # Group summary
    dist.mtx <- dist.fx(MODELS[imod], Dist.type)
    md <- vegan::meandist(dist.mtx, by.group)
    a <- summary(md)
    return(data.frame("." = c("Within Groups", "Between Groups", "Overall"),
                      Average = c(a$W, a$B, a$D)))

  } else if (isTRUE(within.group)) { # Within-group
    dist.mtx <- dist.fx(MODELS[imod], Dist.type)
    md <- vegan::meandist(dist.mtx, by.group)
    return(data.frame(Group = names(attributes(md)$n),
                      N = as.vector(attributes(md)$n),
                      Within.Distance = as.vector(diag(md))))

  } else {                    # Group-wise without summary or within-group
    dist.mtx <- dist.fx(MODELS[imod], Dist.type)
    md <- vegan::meandist(dist.mtx, by.group)
    md <- as.dist(md)
    out.form(md, Out.Format)
  }
}
