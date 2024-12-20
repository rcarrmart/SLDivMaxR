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

  if (is.null(by.group) == TRUE ) {
    if (isTRUE(group.summary) == TRUE) stop("No grouping in `by.group` argument included")

    if (Out.Format == "long" ) {

      if (Dist.type == "prop" & Model != "N") { ## From 'ape' package
        inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                as.matrix = FALSE, gamma = GAMMA)
        Dist2DF.Long(inDist)
      } else
        if (Dist.type == "percent" & Model != "N") {
          inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                   as.matrix = FALSE, gamma = GAMMA) * 100
          Dist2DF.Long(inDist2)
        } else
          if (Dist.type == "nucleotide") {
            inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = "N", variance = Var,
                                     as.matrix = FALSE, gamma = GAMMA)
            Dist2DF.Long(inDist2)}
    }else

      if (Out.Format == "wide") {
        if (Dist.type == "prop" & Model != "N") { ## From 'ape' package
          inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                  as.matrix = FALSE, gamma = GAMMA)
          Dist2DF.Wide(inDist)
        } else
          if (Dist.type == "percent" & Model != "N") {
            inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                     as.matrix = FALSE, gamma = GAMMA) * 100
            Dist2DF.Wide(inDist2)
          } else
            if (Dist.type == "nucleotide") {
              inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = "N", variance = Var,
                                       as.matrix = FALSE, gamma = GAMMA)
              Dist2DF.Wide(inDist2)
            } else
              stop("Wrong 'Out.Format' input. Options are: 'long', 'wide'.")}
  } else

    if (is.null(by.group) == FALSE & isTRUE(group.summary) == FALSE & within.group == FALSE) {
      if (Out.Format == "long" ) {

        if (Dist.type == "prop" & Model != "N") { ## From 'ape' package
          inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                  as.matrix = FALSE, gamma = GAMMA)
          md <- vegan::meandist(inDist, by.group)
          md <- as.dist(md)
          Dist2DF.Long(md)
        } else
          if (Dist.type == "percent" & Model != "N") {
            inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                     as.matrix = FALSE, gamma = GAMMA) * 100
            md <- vegan::meandist(inDist2, by.group)
            md <- as.dist(md)
            Dist2DF.Long(md)
          } else
            if (Dist.type == "nucleotide") {
              inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = "N", variance = Var,
                                       as.matrix = FALSE, gamma = GAMMA)
              md <- vegan::meandist(inDist2, by.group)
              md <- as.dist(md)
              Dist2DF.Long(md)}
      } else

        if (Out.Format == "wide") {
          if (Dist.type == "prop" & Model != "N") { ## From 'ape' package
            inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                    as.matrix = FALSE, gamma = GAMMA)
            md <- vegan::meandist(inDist, by.group)
            md <- as.dist(md)
            Dist2DF.Wide(md)
          } else
            if (Dist.type == "percent" & Model != "N") {
              inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                       as.matrix = FALSE, gamma = GAMMA) * 100
              md <- vegan::meandist(inDist2, by.group)
              md <- as.dist(md)
              Dist2DF.Wide(md)
            } else
              if (Dist.type == "nucleotide") {
                inDist3 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = "N", variance = Var,
                                         as.matrix = FALSE, gamma = GAMMA)
                md <- vegan::meandist(inDist3, by.group)
                md <- as.dist(md)
                Dist2DF.Wide(md)
              } else
                stop("Wrong 'Out.Format' input. Options are: 'long', 'wide'.")}
    }else

      if (isTRUE(group.summary) == TRUE & is.null(by.group) == FALSE & within.group == FALSE) {
        if (Dist.type == "prop" & Model != "N") { ## From 'ape' package
          inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                  as.matrix = FALSE, gamma = GAMMA)
          md <- vegan::meandist(inDist, by.group)
          a <- summary(md)
          data.frame("."= c("Within Groups", "Between Groups", "Overall"), Average = c(a$W,a$B,a$D))
        } else
          if (Dist.type == "percent" & Model != "N"){
            inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                     as.matrix = FALSE, gamma = GAMMA) * 100
            md <- vegan::meandist(inDist2, by.group)
            a <- summary(md)
            data.frame("."= c("Within Groups", "Between Groups", "Overall"), Average = c(a$W,a$B,a$D))
          } else
            if (Dist.type == "nucleotide") {
              inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = "N", variance = Var,
                                       as.matrix = FALSE, gamma = GAMMA)
              md <- vegan::meandist(inDist2, by.group)
              a <- summary(md)
              data.frame("."= c("Within Groups", "Between Groups", "Overall"), Average = c(a$W,a$B,a$D))}
      } else

        if (within.group == TRUE & is.null(by.group) == FALSE) {
          if (Dist.type == "prop" & Model != "N") { ## From 'ape' package
            inDist <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                    as.matrix = FALSE, gamma = GAMMA)
            md <- vegan::meandist(inDist, by.group)
            data.frame(Group = names(attributes(md)$n), N = as.vector(attributes(md)$n),
                       Within.Distance = as.vector(diag(md)))
          } else
            if (Dist.type == "percent" & Model != "N") {
              inDist2 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = (MODELS[imod]), variance = Var,
                                       as.matrix = FALSE, gamma = GAMMA) * 100
              md <- vegan::meandist(inDist2, by.group)
              data.frame(Group = names(attributes(md)$n), N = as.vector(attributes(md)$n),
                         Within.Distance = as.vector(diag(md)))
            } else
              if (Dist.type == "nucleotide") {
                inDist3 <- ape::dist.dna(inDNA, pairwise.deletion = PW.deletion, model = "N", variance = Var,
                                         as.matrix = FALSE, gamma = GAMMA)
                md <- vegan::meandist(inDist3, by.group)
                data.frame(Group = names(attributes(md)$n), N = as.vector(attributes(md)$n),
                           Within.Distance = as.vector(diag(md)))
              } else
                if (within.group == TRUE & is.null(by.group) == TRUE) stop("Missing `by.group` argument") }
}
