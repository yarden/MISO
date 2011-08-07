
matchIso <- function(geneStructure, gene=1, reads, paired=reads$paired,
                     insertProb=NULL, insertStart=0L,
                     normalMean, normalVar, numDevs) {

  if (paired) {
    getreadlength <- function(cig) {
      sum(as.numeric(strsplit(gsub("[0-9]+[ABCDEFGHIJKLNOPQRSTUVWXYZ]", "",
                                   cig), "M")[[1]]))
    }
    readLength=getreadlength(reads$cigar[1])
    if (!is.null(insertProb)) { insertProb <- as.numeric(insertProb) }
    .Call("R_splicing_matchIso_paired", geneStructure, as.integer(gene),
          as.integer(reads$position), as.character(reads$cigar),
          as.integer(readLength), insertProb, as.integer(insertStart),
          as.numeric(normalMean), as.numeric(normalVar),
          as.numeric(numDevs),
          PACKAGE="splicing")
  } else {
    .Call("R_splicing_matchIso", geneStructure, as.integer(gene),
          as.integer(reads$position), as.character(reads$cigar),
          PACKAGE="splicing")
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
