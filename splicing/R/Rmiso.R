
MISO <- function(geneStructure, gene=1L, reads,
                 readLength=getReadLength(reads),
                 overHang=1L, noChains=1, noIterations=5000, noBurnIn=500,
                 noLag=10, start=c("auto", "uniform", "random", "given"),
                 startPsi=NULL, startAlpha=NULL,
                 hyperparameters=rep(1, noIso(geneStructure)[gene]),
                 paired=reads$paired, fragmentProb=NULL, fragmentStart=0L,
                 normalMean, normalVar, numDevs) {

  if (length(readLength) != 1) {
    stop("Variable read length is currently not supported")
  }

  start <- switch(match.arg(start), "auto"=0L, "uniform"=1L, "random"=2L,
                  "given"=3L)
  if (!is.null(startPsi)) {
    startPsi <- structure(as.numeric(startPsi), dim=dim(startPsi))
  }
  if (!is.null(startAlpha)) {
    startAlpha <- structure(as.numeric(startAlpha), dim=dim(startAlpha))
  }
  
  if (!paired) {
    res <- .Call("R_splicing_miso", geneStructure, as.integer(gene), reads,
                 as.integer(readLength), as.integer(noChains),
                 as.integer(noIterations), as.integer(noBurnIn),
                 as.integer(noLag), as.double(hyperparameters),
                 as.integer(overHang), start, startPsi, startAlpha,
                 PACKAGE="splicing")
  } else {
    if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }
    res <- .Call("R_splicing_miso_paired", geneStructure, as.integer(gene),
                 reads, as.integer(readLength), as.integer(noIterations),
                 as.integer(noBurnIn), as.integer(noLag),
                 as.double(hyperparameters), fragmentProb,
                 as.integer(fragmentStart), as.double(normalMean),
                 as.double(normalVar), as.double(numDevs),
                 as.integer(overHang),
                 PACKAGE="splicing")
  }

  res$geneStructure <- selectGenes(geneStructure, gene)
  res
}

## Do all genes, reads are kept on the disk

runMISO <- function(geneStructure, readsfile, 
                    results = c("return", "files", "Rfiles"),
                    resultDir = ".", readLength = NULL,
                    verbose = TRUE, snowCluster = NULL,
                    seed = NULL, ...) {

  results <- match.arg(results)

  if (is.character(geneStructure)) {
    geneStructure <- readGFF3(geneStructure)
  }
  if (!isGFF3(geneStructure)) {
    stop("Invalid gene structure(s), must be a GFF3 object or a filename")
  }

  ids <- geneIds(geneStructure)
  
  get.region <- function(gene) {
    start <- geneStructure$start[geneStructure$gid[gene]+1]
    end <- geneStructure$end[geneStructure$gid[gene]+1]
    seqid <- geneStructure$seqid[gene]
    seqid_str <- geneStructure$seqid_str[seqid+1]
    paste(sep="", seqid_str, ":", start, "-", end)
  }
  
  run.return <- function(misoResult, gene) {
    misoResult
  }

  run.files <- function(misoResult, gene) {
    fname <- paste(sep="", resultDir, "/", ids[gene], ".miso")
    writeMISO(misoResult, fname)
    fname
  }

  run.Rfiles <- function(misoResult, gene) {
    fname <- paste(sep="", resultDir, "/", ids[gene], ".Rdata")
    save(misoResult, file=fname)
    fname
  }
  
  funcs <- list("return"=run.return, "files"=run.files, "Rfiles"=run.Rfiles)

  myclusterExport <- function(cl, list) {
    for (name in list) {
      clusterCall(cl, assign, name, get(name))
    }
  }

  if (!is.null(snowCluster)) {
    require(snow)
    if (!is.null(seed)) {
      clusterCall(snowCluster, set.seed, seed)
    }
    myclusterExport(snowCluster, c("get.region", "readsfile",
                                   "readLength", "verbose", "ids",
                                   "geneStructure", "funcs", "results"))
    myapply <- function(...) parLapply(cl=snowCluster, ...)
  } else {
    myapply <- lapply
  }

  res <- myapply(seq_along(ids), function(x) {
    region <- get.region(x)
    reads <- readSAM(readsfile, region=region)
    rl <- if (is.null(readLength)) {
      getReadLength(reads)
    } else {
      readLength
    }
    if (verbose) { message("Running gene # ", x, ", ", ids[x]) }
    res <- MISO(geneStructure, gene=x, reads=reads, readLength=rl, ...)
    funcs[[results]](res, x)
  })
  
  names(res) <- ids
  res
}

print.MISOresult <- function(x, ...) {
  cat(sep="", "MISO ", geneIds(x$geneStructure), ", ",
      runData(x)$noIso, " i: ",
      paste(collapse=" ", round(postMean(x), 2)), "\n")
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
