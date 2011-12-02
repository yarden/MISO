
matchIso <- function(geneStructure, gene=1, reads, overHang=1L,
                     readLength=0L, paired=isPaired(reads),
                     fragmentProb=NULL, fragmentStart=0L, normalMean,
                     normalVar, numDevs) {

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
                     overHang=1L, paired=isPaired(reads), fragmentProb=NULL,
                     fragmentStart=0L, normalMean, normalVar, numDevs) {

  if (length(readLength) != 1) {
    stop("Variable read length is currently not supported")
  }
  
  if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }

  if (!paired) { 
    .Call("R_splicing_solve_gene", geneStructure, as.integer(gene),
          as.integer(readLength), as.integer(reads$position),
          as.character(reads$cigar), as.integer(overHang),
          PACKAGE="splicing")
  } else {
    .Call("R_splicing_solve_gene_paired", geneStructure, as.integer(gene),
          as.integer(readLength), as.integer(overHang),
          as.integer(reads$position), as.character(reads$cigar), fragmentProb,
          as.integer(fragmentStart), as.double(normalMean),
          as.double(normalVar), as.double(numDevs),
          PACKAGE="splicing")
  }
}

## TODO: writeFunction

runLinear <- function(geneStructure, readsfile,
                      results=c("return", "files", "Rfiles"),
                      resultDir=".", overWrite=FALSE, readLength=NULL,
                      verbose=TRUE, snowCluster=NULL, seed=NULL, ...) {

  run(runFunction=solveIso, writeFunction=NULL, geneStructure=geneStructure,
      readsfile=readsfile, results=results, resultDir=resultDir,
      overWrite=overWrite, readLength=readLength, verbose=verbose,
      snowCluster=snowCluster, seed=seed, ...)
}
