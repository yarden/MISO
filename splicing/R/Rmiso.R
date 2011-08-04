
MISO <- function(geneStructure, gene=1L, reads, readLength,
                 noIterations=5000, noBurnIn=500, noLag=10,
                 hyperparameters=rep(1,noIso(geneStructure)[1]),
                 paired=reads$paired, insertProb=NULL, insertStart,
                 normalMean, normalVar, numDevs) {

  if (!paired) {
    res <- .Call("R_splicing_miso", geneStructure, as.integer(gene), reads,
                 as.integer(readLength), as.integer(noIterations),
                 as.integer(noBurnIn), as.integer(noLag),
                 as.double(hyperparameters),
                 PACKAGE="splicing")
  } else {
    if (!is.null(insertProb)) { insertProb <- as.double(insertProb) }
    res <- .Call("R_splicing_miso_paired", geneStructure, as.integer(gene),
                 reads, as.integer(readLength), as.integer(noIterations),
                 as.integer(noBurnIn), as.integer(noLag),
                 as.double(hyperparameters), insertProb,
                 as.integer(insertStart), as.double(normalMean),
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
