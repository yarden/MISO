
matchIso <- function(geneStructure, gene=1, reads, overHang=1L,
                     readLength=0L, paired=isPaired(reads),
                     fragmentProb=NULL, fragmentStart=0L, normalMean=NA,
                     normalVar=NA, numDevs=4) {

  if (paired) {
    getreadlength <- function(cig) {
      sum(as.numeric(strsplit(gsub("[0-9]+[ABCDEFGHIJKLNOPQRSTUVWXYZ]", "",
                                   cig), "M")[[1]]))
    }
    if (readLength==0) { readLength=getreadlength(reads$cigar[1]) }
    if (!is.null(fragmentProb)) { fragmentProb <- as.numeric(fragmentProb) }
    .Call("R_splicing_matchIso_paired", geneStructure, as.integer(gene),
          as.integer(reads$position), as.character(reads$cigar),
          as.integer(readLength), as.integer(overHang), fragmentProb,
          as.integer(fragmentStart), as.numeric(normalMean),
          as.numeric(normalVar), as.numeric(numDevs),
          PACKAGE="splicing")
  } else {
    .Call("R_splicing_matchIso", geneStructure, as.integer(gene),
          as.integer(reads$position), as.character(reads$cigar),
          as.integer(overHang), as.integer(readLength),
          PACKAGE="splicing")
  }
}

solveIso <- function(geneStructure, gene=1L, reads,
                     readLength=getReadLength(reads),
                     overHang=1L, scale=TRUE,
                     paired=isPaired(reads), fast=FALSE, fragmentProb=NULL,
                     fragmentStart=0L, normalMean=NA, normalVar=NA,
                     numDevs=4) {

  if (length(readLength) != 1) {
    stop("Variable read length is currently not supported")
  }
  
  if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }

  if (!paired) { 
    res <- .Call("R_splicing_solve_gene", geneStructure, as.integer(gene),
                 as.integer(readLength), as.integer(reads$position),
                 as.character(reads$cigar), as.integer(overHang),
                 as.logical(scale), PACKAGE="splicing")
  } else {
    res <- .Call("R_splicing_solve_gene_paired", geneStructure,
                 as.integer(gene), as.integer(readLength),
                 as.integer(overHang), as.logical(scale),
                 as.integer(reads$position), as.character(reads$cigar),
                 as.logical(fast), fragmentProb, as.integer(fragmentStart),
                 as.double(normalMean), as.double(normalVar),
                 as.double(numDevs), PACKAGE="splicing")
  }

  res$geneStructure <- selectGenes(geneStructure, gene)
  res$classTemplates <- ifelse(res$assignment > 0, 1, 0)
  res  
}

writeLinear <- function(linResult, file) {
  .Call("R_splicing_writelinear", linResult, as.character(file),
        PACKAGE="splicing")
  invisible(NULL)
}

runLinear <- function(geneStructure, readsfile,
                      results=c("return", "files", "Rfiles"),
                      resultDir=".", overWrite=FALSE, readLength=NULL,
                      verbose=TRUE, snowCluster=NULL, seed=NULL,
                      dropBadCigar=FALSE, ...) {

  run(runFunction=solveIso, writeFunction=writeLinear,
      geneStructure=geneStructure,
      readsfile=readsfile, results=results, resultDir=resultDir,
      overWrite=overWrite, readLength=readLength, verbose=verbose,
      snowCluster=snowCluster, seed=seed, dropBadCigar=dropBadCigar, ...)
}

## Linear bias, to be rewritten...

solveIsoLinBias <- function(geneStructure, gene=1L, reads,
                            readLength=getReadLength(reads), overHang=1L,
                            scale=TRUE, paired=isPaired(reads), fast=FALSE,
                            fragmentProb=NULL, fragmentStart=0L,
                            normalMean=NA, normalVar=NA, numDevs=4) {

  require(nnls)
  
  if (length(readLength) != 1) {
    stop("Variable read length is currently not supported")
  }
  if (!is.null(fragmentProb)) {
    fragmentProb <- as.double(fragmentProb)
  }

  ## Assignment matrix one
  am0 <- assignmentMatrix(geneStructure=geneStructure, gene=gene,
                          readLength=readLength, overHang=overHang,
                          bias=0, paired=paired, fast=fast,
                          fragmentProb=fragmentProb,
                          fragmentStart=fragmentStart,
                          normalMean=normalMean, normalVar=normalVar,
                          numDevs=numDevs)

  ## Assignment matrix two
  am1 <- assignmentMatrix(geneStructure=geneStructure, gene=gene,
                          readLength=readLength, overHang=overHang,
                          bias=1, paired=paired, fast=fast,
                          fragmentProb=fragmentProb,
                          fragmentStart=fragmentStart,
                          normalMean=normalMean, normalVar=normalVar,
                          numDevs=numDevs)

  ## This shouldn't happen, but check anyway
  if (!all(colnames(am0)==colnames(am1))) {
    stop("Assignment matrix columns do not match")
  }
  
  ## Match the reads
  mat <- matchIso(geneStructure=geneStructure, gene=gene, reads=reads,
                  overHang=overHang, readLength=readLength,
                  paired=paired, fragmentProb=fragmentProb,
                  fragmentStart=fragmentStart, normalMean=normalMean,
                  normalVar=normalVar, numDevs=numDevs)
  if (paired) { mat <- mat[[1]] }    
  matvec <- table(apply((mat!=0)+0L, 2, paste,
                        collapse=""))[colnames(am0)]
  is.na(matvec) <- 0

  ## Test identifyability
  sol0 <- nnls(t(0*am1 + am0), cbind(matvec))
  sol1 <- nnls(t(1*am1 + am0), cbind(matvec))
  if (abs(sum(sol0$residuals^2) - sum(sol1$residuals^2)) < 1e-8) {
    stop("Linear bias parameter is not identifyable")
  }
  
  f <- function(x, param) {
    am <- x * param$am1 + param$am0
    sol <- nnls(t(am), cbind(param$matvec))
    res <- sum(sol$residuals^2)
    res
  }

  optsol <- optimize(f=f, c(0,1), tol=.Machine$double.eps,
                     param=list(am0=am0, am1=am1, matvec=matvec))

  ## Final solution
  n <- function(x) x/sum(x)
  am <- optsol$minimum * am1 + am0
  finsol <- nnls(t(am), cbind(matvec))
  list(a=optsol$minimum, expression=n(finsol$x),
       residuals=finsol$residuals, finsol=finsol)  
}

# Quadratic bias, to be rewritten

solveIsoQuadBias <- function(geneStructure, gene=1L, reads,
                             readLength=getReadLength(reads), overHang=1L,
                             scale=TRUE, paired=isPaired(reads),
                             fast=FALSE, fragmentProb=NULL,
                             fragmentStart=0L, normalMean=NA,
                             normalVar=NA, numDevs=4,
                             lower=-Inf, upper=Inf, ...) {

  require(nnls)

  if (length(readLength) != 1) {
    stop("Variable read length is currently not supported")
  }
  if (!is.null(fragmentProb)) {
    fragmentProb <- as.double(fragmentProb)
  }

  ## Assignment matrix one
  am0 <- assignmentMatrix(geneStructure=geneStructure, gene=gene,
                          readLength=readLength, overHang=overHang,
                          bias=0, paired=paired, fast=fast,
                          fragmentProb=fragmentProb,
                          fragmentStart=fragmentStart,
                          normalMean=normalMean, normalVar=normalVar,
                          numDevs=numDevs)

  ## Assignment matrix two
  am1 <- assignmentMatrix(geneStructure=geneStructure, gene=gene,
                          readLength=readLength, overHang=overHang,
                          bias=1, paired=paired, fast=fast,
                          fragmentProb=fragmentProb,
                          fragmentStart=fragmentStart,
                          normalMean=normalMean, normalVar=normalVar,
                          numDevs=numDevs)

  ## Assignment matrix three
  am2 <- assignmentMatrix(geneStructure=geneStructure, gene=gene,
                          readLength=readLength, overHang=overHang,
                          bias=2, paired=paired, fast=fast,
                          fragmentProb=fragmentProb,
                          fragmentStart=fragmentStart,
                          normalMean=normalMean, normalVar=normalVar,
                          numDevs=numDevs)

  ## This shouldn't happen, but check anyway
  if (!all(colnames(am0)==colnames(am1)) ||
      !all(colnames(am0)==colnames(am2))) {
    stop("Assignment matrix columns do not match")
  }

  ## Match the reads
  mat <- matchIso(geneStructure=geneStructure, gene=gene, reads=reads,
                  overHang=overHang, readLength=readLength,
                  paired=paired, fragmentProb=fragmentProb,
                  fragmentStart=fragmentStart, normalMean=normalMean,
                  normalVar=normalVar, numDevs=numDevs)
  if (paired) { mat <- mat[[1]] }    
  matvec <- table(apply((mat!=0)+0L, 2, paste,
                        collapse=""))[colnames(am0)]
  is.na(matvec) <- 0

  ## Test identifyability. Fix 'a' to two selected values and
  ## then check the residuals. If they are the same for both
  ## fixed 'a' values, then the model is not idendifyable
  f <- function(x, param) {
    am <- x * param$am2 + param$b * param$am1 + param$am0
    sol <- nnls(t(am), cbind(param$matvec))
    res <- sum(sol$residuals^2)
    res
  }
  sol0 <- optimize(f=f, c(-100,100), tol=.Machine$double.eps,
                   param=list(am0=am0, am1=am1, am2=am2, b=0, matvec=matvec))
  sol1 <- optimize(f=f, c(-100,100), tol=.Machine$double.eps,
                   param=list(am0=am0, am1=am1, am2=am2, b=1, matvec=matvec))
  if (abs(sol0$objective - sol1$objective) < 1e-5) {
    stop("Quadratic bias parameters are not identifyable")
  }

  ## OK, ready to optimize quadratic parameters
  fn <- function(x) {
    am <- x[1] * am2 + x[2] * am1 + am0
    sol <- nnls(t(am), cbind(matvec))
    res <- sum(sol$residuals^2)
    res
  }

  optsol <- optim(c(0,0), fn, lower=lower, upper=upper, ...)

  ## Final solution
  n <- function(x) x/sum(x)
  am <- optsol$par[1] * am2 + optsol$par[2] * am1 + am0
  finsol <- nnls(t(am), cbind(matvec))
  list(a=optsol$par[1], b=optsol$par[2], expression=n(finsol$x),
       residuals=finsol$residuals, finsol=finsol)
}
                             
                             
