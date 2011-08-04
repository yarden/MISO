
matchIso <- function(geneStructure, gene=1, reads) {

  res <- .Call("R_splicing_matchIso", geneStructure, as.integer(gene),
               as.integer(reads$position), as.character(reads$cigar),
               PACKAGE="splicing")

  paired <- attr(reads, "paired")
  if (!is.null(paired) && paired) {
    res <- (res[,seq(1,ncol(res),by=2)] & res[,seq(2,ncol(res),by=2)]) + 0
    res <- apply(res, 2, paste, collapse="")
    res <- rep(res, each=2)
    res
  } else {
    apply(res, 2, paste, collapse="")
  }
}

solveIsoOne <- function(geneStructure, gene=1L, reads, readLength=33L,
                        paired=FALSE, insertProb=NULL, insertStart=0L,
                        normalMean, normalVar, numDev) {

  if (!is.null(insertProb)) { insertProb <- as.double(insertProb) }

  if (!paired) { 
    .Call("R_splicing_solve_gene", geneStructure, as.integer(gene),
          as.integer(readLength), as.integer(reads$position),
          as.character(reads$cigar),
          PACKAGE="splicing")
  } else {
    .Call("R_splicing_solve_gene_paired", geneStructure, as.integer(gene),
          as.integer(readLength), as.integer(reads$position),
          as.character(reads$cigar), insertProb,
          as.integer(insertStart), as.double(normalMean),
          as.double(normalVar), as.double(numDev),
          PACKAGE="splicing")
  }
}
