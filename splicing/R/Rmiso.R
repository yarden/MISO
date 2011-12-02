
MISO <- function(geneStructure, gene=1L, reads,
                 readLength=getReadLength(reads),
                 overHang=1L, noChains=2, noIterations=5000,
                 maxIterations=100000, noBurnIn=500,
                 noLag=10, start=c("auto", "uniform", "random", "given",
                             "linear"), startPsi=NULL, 
                 stopCond=c("convMean", "fixedno"),
                 hyperparameters=rep(1, noIso(geneStructure)[gene]),
                 paired=isPaired(reads), fragmentProb=NULL, fragmentStart=0L,
                 normalMean, normalVar, numDevs) {

  if (length(readLength) != 1) {
    stop("Variable read length is currently not supported")
  }

  useq <- unique(eachSeqName(reads))
  if (length(useq) != 1) {
    stop("Reads from multiple sequences")
  }
  if (useq != seqIds(geneStructure)[gene]) {
    stop("Sequence ids no not match")
  }
  
  start <- switch(match.arg(start), "auto"=0L, "uniform"=1L, "random"=2L,
                  "given"=3L, "linear"=4L)
  stopCond <- switch(match.arg(stopCond), "fixedno"=0L, convMean=1L);
  
  if (!is.null(startPsi)) {
    startPsi <- structure(as.numeric(startPsi), dim=dim(startPsi))
  }
  
  if (!paired) {
    res <- .Call("R_splicing_miso", geneStructure, as.integer(gene), reads,
                 as.integer(readLength), as.integer(noChains),
                 as.integer(noIterations), as.integer(maxIterations),
                 as.integer(noBurnIn), as.integer(noLag),
                 as.double(hyperparameters),
                 as.integer(overHang), start, startPsi, stopCond,
                 PACKAGE="splicing")
  } else {
    if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }
    res <- .Call("R_splicing_miso_paired", geneStructure, as.integer(gene),
                 reads, as.integer(readLength), as.integer(noChains),
                 as.integer(noIterations), as.integer(maxIterations),
                 as.integer(noBurnIn), as.integer(noLag),
                 as.double(hyperparameters), fragmentProb,
                 as.integer(fragmentStart), as.double(normalMean),
                 as.double(normalVar), as.double(numDevs),
                 as.integer(overHang), start, startPsi, stopCond,
                 PACKAGE="splicing")
  }

  res$geneStructure <- selectGenes(geneStructure, gene)
  res
}

## Do all genes, reads are kept on the disk

runMISO <- function(geneStructure, readsfile,
                    results = c("return", "files", "Rfiles"),
                    resultDir = ".", overWrite = FALSE, readLength = NULL,
                    verbose = TRUE, snowCluster = NULL,
                    seed = NULL, ...) {

  run(runFunction=MISO, writeFunction=writeMISO, geneStructure=geneStructure,
      readsfile=readsfile, results=results, resultDir=resultDir,
      overWrite=overWrite, readLength=readLength, verbose=verbose,
      snowCluster=snowCluster, seed=seed, ...)
}     

print.MISOresult <- function(x, ...) {
  ci <- t(confint(x))-postMean(x)
  ci <- apply(abs(ci), 1, max)
  cat(sep="", "MISO ", geneIds(x$geneStructure), ", ",
      runData(x)$noIso, " i: ",
      paste(collapse=", ", sep="", round(postMean(x), 2), " (+-",
            round(ci,2), ")"), "\n")
}

writeMISO <- function(misoResult, file)
  UseMethod("writeMISO")

writeMISO.MISOresult <- function(misoResult, file) {
  if (!inherits(misoResult, "MISOresult")) {
    stop("This does not seem to be a MISO result")
  }
  .Call("R_splicing_writemiso", misoResult, as.character(file),
        PACKAGE="splicing")
  invisible(NULL)
}

readMISO <- function(file) {
  lines <- readLines(file)
  lines <- lines[ lines != "" ]
  headers <- c(grep("^\\[[a-zA-Z0-9_-]*\\]$", lines), length(lines)+1)
  records <- mapply(headers[-length(headers)], headers[-1],
                    FUN=function(f, t) lines[(f+1):(t-1)])
  names(records) <- sub("[", "", fixed=TRUE,
                        sub("]", "", fixed=TRUE,
                            lines[headers[-length(headers)]]))
  
  ## runData
  header <- strsplit(records[["runData"]], ":")
  head1 <- sapply(header, "[[", 1)
  head2 <- lapply(header, "[[", 2)
  names(head2) <- head1
  runData <- list(noIso=as.numeric(head2[['noIso']]),
                  noIters=as.numeric(head2[['noIters']]),
                  noBurnIn=as.numeric(head2[['noBurnIn']]),
                  noLag=as.numeric(head2[['noLag']]),
                  noAccepted=as.numeric(head2[['noAccepted']]),
                  noRejected=as.numeric(head2[['noRejected']]))

  ## samples
  tc <- textConnection(records[["samples"]])
  samples <- t(matrix(scan(tc, what=0.0, quiet=TRUE), nc=runData$noIso,
                      byrow=TRUE))
  close(tc)

  ## logLik
  tc <- textConnection(records[["logLik"]])
  logLik <- scan(tc, what=0.0, quiet=TRUE)
  close(tc)

  ## matchMatrix
  tc <- textConnection(records[["matchMatrix"]])
  matchMatrix <- t(matrix(scan(tc, what=0.0, quiet=TRUE), nc=runData$noIso,
                          byrow=TRUE))
  close(tc)

  ## classTemplates
  tc <- textConnection(records[["classTemplates"]])
  classTemplates <- matrix(scan(tc, what=0.0, quiet=TRUE), nr=runData$noIso)
  close(tc)  

  ## classCounts
  tc <- textConnection(records[["classCounts"]])
  classCounts <- scan(tc, what=0.0, quiet=TRUE)
  close(tc)

  ## geneStructure
  tmp <- tempfile()
  cat(records[["geneStructure"]], sep="\n", file=tmp)
  geneStructure <- readGFF3(tmp)
  
  res <- list(samples=samples, logLik=logLik, matchMatrix=matchMatrix,
              classTemplates=classTemplates, classCounts=classCounts,
              runData=runData, geneStructure=geneStructure)
  class(res) <- "MISOresult"
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

confint.MISOresult <- function(misoResult, level=0.95) {
  if (level <= 0 || level >= 1) {
    stop("`level' must be greater than 0 and less than 1")
  }
  alpha <- 1-level
  lw.idx <- round(alpha/2 * ncol(misoResult$samples))
  up.idx <- round((1-alpha/2) * ncol(misoResult$samples))
  ordered <- apply(misoResult$samples, 1, sort)
  rbind(ordered[lw.idx,], ordered[up.idx,])
}
