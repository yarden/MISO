
assignmentMatrix <- function(geneStructure, gene=1L, readLength=33L,
                             paired=FALSE, fragmentProb=NULL, fragmentStart=0L,
                             normalMean, normalVar, numDevs) {

  if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }
  
  if (!paired) { 
    res <- .Call("R_splicing_assignment_matrix", geneStructure,
                 as.integer(gene), as.integer(readLength),
                 PACKAGE="splicing")
  } else {
    res <- .Call("R_splicing_paired_assignment_matrix", geneStructure,
                 as.integer(gene), as.integer(readLength),
                 fragmentProb, as.integer(fragmentStart),
                 as.double(normalMean), as.double(normalVar),
                 as.double(numDevs),
                 PACKAGE="splicing")
  }
  
  colnames(res) <- apply(res, 2, function(x) {
    paste(ifelse(x==0, '0', '1'), collapse="")
  })
  res
}

condIso <- function(geneStructure, gene=1, readLength,
                    type=c("relative", "absolute"),
                    norm=c("2","1","inf"), paired=FALSE,
                    fragmentProb=NULL, fragmentStart=0L, normalMean,
                    normalVar, numDevs) {

  type <- switch(match.arg(type), "relative"=0, "absolute"=1)
  norm <- switch(match.arg(norm), "2"=0, "1"=1, "inf"=2)
  if (!is.null(fragmentProb)) { fragmentProb <- as.numeric(fragmentProb) }
  if (missing(normalMean)) { normalMean <- 0 }
  if (missing(normalVar)) { normalVar <- 0 }
  if (missing(numDevs)) { numDevs <- 0 }

  .Call("R_splicing_gene_complexity", geneStructure, as.integer(gene),
        as.integer(readLength), as.integer(type), as.integer(norm),
        as.logical(paired), fragmentProb, as.integer(fragmentStart),
        as.double(normalMean), as.double(normalVar), as.double(numDevs),
        PACKAGE="splicing")
} 
