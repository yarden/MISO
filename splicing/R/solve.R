
matchIso <- function(geneStructure, gene=1, reads, paired=reads$paired,
                     fragmentProb=NULL, fragmentStart=0L,
                     normalMean, normalVar, numDevs) {

  if (paired) {
    getreadlength <- function(cig) {
      sum(as.numeric(strsplit(gsub("[0-9]+[ABCDEFGHIJKLNOPQRSTUVWXYZ]", "",
                                   cig), "M")[[1]]))
    }
    readLength=getreadlength(reads$cigar[1])
    if (!is.null(fragmentProb)) { fragmentProb <- as.numeric(fragmentProb) }
    .Call("R_splicing_matchIso_paired", geneStructure, as.integer(gene),
          as.integer(reads$position), as.character(reads$cigar),
          as.integer(readLength), fragmentProb, as.integer(fragmentStart),
          as.numeric(normalMean), as.numeric(normalVar),
          as.numeric(numDevs),
          PACKAGE="splicing")
  } else {
    .Call("R_splicing_matchIso", geneStructure, as.integer(gene),
          as.integer(reads$position), as.character(reads$cigar),
          PACKAGE="splicing")
  }
}

solveIso <- function(geneStructure, gene=1L, reads, readLength=33L,
                     paired=FALSE, fragmentProb=NULL, fragmentStart=0L,
                     normalMean, normalVar, numDevs) {

  if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }

  if (!paired) { 
    .Call("R_splicing_solve_gene", geneStructure, as.integer(gene),
          as.integer(readLength), as.integer(reads$position),
          as.character(reads$cigar),
          PACKAGE="splicing")
  } else {
    .Call("R_splicing_solve_gene_paired", geneStructure, as.integer(gene),
          as.integer(readLength), as.integer(reads$position),
          as.character(reads$cigar), fragmentProb,
          as.integer(fragmentStart), as.double(normalMean),
          as.double(normalVar), as.double(numDevs),
          PACKAGE="splicing")
  }
}
