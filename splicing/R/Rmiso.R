
# TODO: better API

MISO <- function(geneStructure, gene=1L, reads, readLength,
                 noIterations=5000, noBurnIn=500, noLag=10,
                 hyperparameters=rep(1,noIso(geneStructure)[1]),
                 paired=reads$paired, fragmentProb=NULL, fragmentStart,
                 normalMean, normalVar, numDevs) {

  if (!paired) {
    res <- .Call("R_splicing_miso", geneStructure, as.integer(gene), reads,
                 as.integer(readLength), as.integer(noIterations),
                 as.integer(noBurnIn), as.integer(noLag),
                 as.double(hyperparameters),
                 PACKAGE="splicing")
  } else {
    if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }
    res <- .Call("R_splicing_miso_paired", geneStructure, as.integer(gene),
                 reads, as.integer(readLength), as.integer(noIterations),
                 as.integer(noBurnIn), as.integer(noLag),
                 as.double(hyperparameters), fragmentProb,
                 as.integer(fragmentStart), as.double(normalMean),
                 as.double(normalVar), as.double(numDevs),
                 PACKAGE="splicing")
  }

  res
}

MISO.Trinity <- function(matchMatrix, fragmentLength=NULL, isoLength,
                         readLength, noIterations=5000, noBurnIn=500,
                         noLag=10, hyperparameters=rep(1, length(isoLength)),
                         paired=FALSE, fragmentProb=NULL, fragmentStart,
                         normalMean, normalVar, numDevs) {
  if (!paired) {
    res <- .Call("R_splicing_miso_trinity", matchMatrix,
                 as.integer(isoLength), as.integer(readLength),
                 as.integer(noIterations), as.integer(noBurnIn),
                 as.integer(noLag), as.double(hyperparameters),
                 PACKAGE="splicing")
  } else {
    if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }
    res <- .Call("R_splicing_miso_paired_trinity",
                 matchMatrix, fragmentLength,
                 as.integer(isoLength), as.integer(readLength),
                 as.integer(noIterations), as.integer(noBurnIn),
                 as.integer(noLag), as.double(hyperparameters), fragmentProb,
                 as.integer(fragmentStart), as.double(normalMean),
                 as.double(normalVar), as.double(numDevs),
                 PACKAGE="splicing")
  }

  res
}

postMean <- function(misoResult)
  UseMethod("postMean")

postMean.MISOresult <- function(misoResult) {
  if (!inherits(misoResult, "MISOresult")) {
    stop("Does not seem to be a MISO result")
  }
  rowMeans(misoResult$samples)
}

logLik.MISOresult <- function(object, ...) {
  object$logLik
}

runData <- function(misoResult)
  UseMethod("runData")

runData.MISOresult <- function(misoResult) {
  misoResult$runData
}
