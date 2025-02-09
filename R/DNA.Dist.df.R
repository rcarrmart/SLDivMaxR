#' Single locus divergence data.frame
#'
#' A function that applies `ape::dist.dna()` but changes `"dist"` class output from `ape::dist.dna()` (or similar data) to `data.frame`.
#'
#' @param inDNA The input aligned sequences in fasta format.
#'
#' @param Dist.type Specifies the distance output format. Options are: `"prop"`: Proportions; `"percentage"`: Percentages; and `"nucleotide"`: Number of nucleotides that are different.
#'
#' @param Model Specifies the evolutionary model to use, as in `ape::dist.dna`. Options include: `"RAW"`, `"JC69"`, `"K80"` (default), `"F81"`, `"K81"`, `"F84"`, `"T92"`, `"TN93"`, `"GG95"`, `"LOGDET"`, `"BH87"`, `"PARALIN"`, `"N"`, `"TS"`, `"TV"`, `"INDEL"`, and `"INDELBLOCK"`. Note: If `Model = "N"` and `Dist.type` is set to `"prop"` or `"percentage"`, an error will be returned.
#'
#' @param GAMMA Same as the `gamma` parameter in `ape::dist.dna()`.
#'
#' @param PW.deletion Specifies whether pairwise deletion should be considered (`TRUE`, default) or not (`FALSE`).
#'
#' @param Var If `TRUE`, the output will include calculated variation. Default is `FALSE`.
#'
#' @param by.group A vector for grouping sequences into different categories (e.g., populations, species, lineages). The order should match the sequence alignment and be used with `vegan::meandist()`.
#'
#' @param group.summary Applies `vegan::meandist()` to generate a `data.frame` with within-group mean distance, between-group mean distance, and overall distance. A vector with discrete groupings should be provided in `by.group`.
#'
#' @param within.group Applies `vegan::meandist()` to create a `data.frame` with each group’s mean divergence and the number of samples ("N") by group. A vector with discrete groupings should be provided in `by.group`. Groups with `N = 1` will return `NA`s.
#'
#' @returns A dataframe with the mean distance between groups or within group distances.
#'          Note: `DNA.Dist.df` should not be used as a species delimitation analysis as it does not explicitly tests a hypothesis.
#'
#' @examples
#' # Distance matrix class == data.frame of all specimens
#' DNA.Dist.df(test.DNA, Out.Format = "wide")
#'
#' # same as above but with JC96 model
#' DNA.Dist.df(test.DNA, Out.Format = "wide", Model = "JC69")
#'
#' # as above, but long pair-wise comparison format
#' DNA.Dist.df(test.DNA, Out.Format = "long", Model = "JC69")
#'
#' ####
#' # Defining groupings
#' ####
#'
#' #create data.frame with specimen and grouping (e.g., species) IDs
#' Species.info <- data.frame(Specimen = paste0("Ind",LETTERS[1:9]),
#'                            Species = paste0("Sp",c(rep(1,3),rep(2,3),rep(3,3))),
#'                            Seq.names = names(test.DNA))
#'
#' Species.info  # Check that they match.
#' # The funtion uses the order of the `by.group` vector to group the sequences
#'
#' DNA.Dist.df(test.DNA, Out.Format = "long", Model = "JC69", by.group = Species.info$Species)
#'
#' # Reorganize the output
#' Comparison <- DNA.Dist.df(test.DNA, Out.Format = "long", Model = "JC69",
#'                           by.group = Species.info$Species)
#' Comparison <- data.frame(Pairs = paste(Comparison$x1, Comparison$x2, sep = "_vs_"),
#'                          Distance = Comparison$Distance)
#' Comparison
#'
#' #wide format
#' WDist <- DNA.Dist.df(test.DNA, Out.Format = "wide", Model = "JC69",
#'                       by.group = Species.info$Species)
#' WDist
#'
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
  if (is.na(imod))  {
                    cat("Model provided:", Model, "\n")
                    stop(paste("'Model' must be one of:", paste("\"", MODELS,"\"", sep = "", collapse = " "),
                              ". See ape::dist.dna for more information.", sep = ""))
                    }

  DT <- c("percent", "prop", "nucleotide")
  iDT <- pmatch((Dist.type), DT)
  if (is.na(iDT)) stop("Wrong 'Dist.type' input. Options are: 'prop', 'percent', 'nucleotide'.")

  if (Dist.type != "nucleotide" & Model == "N")
    stop("Change 'Dist.type' to 'nucleotide' if input for Model is 'N'.")
  if (!is.null(by.group) && length(by.group) != length(inDNA)) stop("`by.group` vector and the DNA list are of different lengths.")
  if (is.na(by.group) == TRUE) stop("Missing data on `by.group` vector")

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
