#' Group Divergence Standard Error
#'
#' Calculates standard error of specified between and within group genetic divergences.
#'
#'
#' @param inDNA The input aligned sequences in fasta format.
#' @param Dist.type Specifies the distance output format. Options are `"prop"`: Proportions; `"percentage"`: Percentages; and `"nucleotide"`: Number of nucleotides that differ.
#' @param Model Specifies the evolutionary model to use, as in ape::dist.dna. Options include: `"RAW"`, `"JC69"`, `"K80"` (default), `"F81"`, `"K81"`, `"F84"`, `"T92"`, `"TN93"`, `"GG95"`, `"LOGDET"`, `"BH87"`, `"PARALIN"`, `"N"`, `"TS"`, `"TV"`, `"INDEL"`, and `"INDELBLOCK"`. Note: If `Model = "N"` and `Dist.type` is set to `"prop"` or `"percentage"`, an error will be returned.
#' @param GAMMA Same as `gamma` in `ape::dist.dna()`.
#' @param PW.deletion Specifies whether pairwise deletion should be considered (TRUE, default) or not (FALSE).
#' @param by.group A vector for grouping sequences into different categories (e.g., populations, species, lineages). The vector should follow the same order as the sequence alignment.
#' @param within.group Creates a data.frame with the within-group standard error and the number of samples ("N") by group. A vector with discrete groupings should be provided in by.group. Groups with N = 1 will return NAs.
#'
#' @returns A data.frame with the standard errors of the mean distance between or within groups distances.
#'
#' @examples
#' # Defining groupings
#'
#' #create data.frame with specimen and grouping (e.g., species) IDs
#' Species.info <- data.frame(Specimen = paste0("Ind",LETTERS[1:9]),
#'                             Species = paste0("Sp",c(rep(1,3),rep(2,3),rep(3,3))),
#'                             Seq.names = names(test.DNA))
#'
#' Species.info  # Check that they match.
#' # The funtion uses the order of the `by.group` vector to group the sequences
#'
#' DNA.GroupDist.SE(test.DNA, Out.Format = "long", Model = "JC69",
#'                  by.group = Species.info$Species)
#'
# Reorganize the output
#' CompSE<- DNA.GroupDist.SE(test.DNA, Out.Format = "long", Model = "JC69",
#'                             by.group = Species.info$Species)
#' CompSE <- data.frame(Pairs = paste(CompSE$x1, CompSE$x2, sep = "_vs_"),
#'                      SE = CompSE$Standard_Error)
#' CompSE
#'
#' #wide format
#' WDistSE <- DNA.GroupDist.SE(test.DNA, Out.Format = "wide", Model = "JC69",
#'                             by.group = Species.info$Species)
#' WDistSE
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
