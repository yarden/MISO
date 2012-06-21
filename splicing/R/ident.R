
assignmentMatrix <- function(geneStructure, gene=1L, readLength,
                             overHang=1L, paired=FALSE, fast=FALSE,
                             fragmentProb=NULL, fragmentStart=0L,
                             normalMean=NA, normalVar=NA, numDevs=4) {

  if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }
  
  if (!paired) { 
    res <- .Call("R_splicing_assignment_matrix", geneStructure,
                 as.integer(gene), as.integer(readLength),
                 as.integer(overHang),
                 PACKAGE="splicing")
  } else {
    res <- .Call("R_splicing_paired_assignment_matrix", geneStructure,
                 as.integer(gene), as.integer(readLength),
                 as.integer(overHang), as.logical(fast), fragmentProb,
                 as.integer(fragmentStart), as.double(normalMean),
                 as.double(normalVar), as.double(numDevs),
                 PACKAGE="splicing")
  }
  
  colnames(res) <- apply(res, 2, function(x) {
    paste(ifelse(x==0, '0', '1'), collapse="")
  })
  res
}

geneComplexity <- function(geneStructure, gene=1, readLength, overHang=1L,
                           type=c("relative", "absolute"),
                           norm=c("2","1","inf"), paired=FALSE, fast=FALSE,
                           fragmentProb=NULL, fragmentStart=0L, normalMean=NA,
                           normalVar=NA, numDevs=4) {

  type <- switch(match.arg(type), "relative"=0, "absolute"=1)
  norm <- switch(match.arg(norm), "2"=0, "1"=1, "inf"=2)
  if (!is.null(fragmentProb)) { fragmentProb <- as.numeric(fragmentProb) }

  .Call("R_splicing_gene_complexity", geneStructure, as.integer(gene),
        as.integer(readLength), as.integer(overHang), as.integer(type),
        as.integer(norm), as.logical(paired), as.logical(fast), fragmentProb,
        as.integer(fragmentStart), as.double(normalMean),
        as.double(normalVar), as.double(numDevs),
        PACKAGE="splicing")
} 
