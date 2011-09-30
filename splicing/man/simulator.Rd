\name{Read simulator}
\alias{simulateReads}
\alias{realPsi}
\title{Simulating RNASeq reads from real or artifical genes}
\description{
  TODO
}
\usage{
simulateReads (geneStructure, gene = 1, expression, noReads, readLength, 
    paired = FALSE, chrname = seqIds(geneStructure)[gene],
    chrlen = geneLength(geneStructure)[gene], 
    qname = "read-\%s", fragmentProb = NULL, fragmentStart = 0L, 
    normalMean, normalVar, numDevs)
realPsi(geneStructure, gene = 1, reads, readLength=getReadLength(reads))
}
\arguments{
  \item{geneStructure}{A GFF3 object containing one of more genes.}
  \item{gene}{The numeric id of the gene for which the reads are
    generated. If the GFF3 object contains a single gene, then this must
    be 1.}
  \item{expression}{A numeric vector giving the expression levels
    for the isoforms.}
  \item{noReads}{The number of reads to generate. For paired-end reads
    this means the number of read pairs.}
  \item{readLength}{Numeric scalar, the read length.}
  \item{paired}{Logical scalar, whether to generate paired reads.}
  \item{chrname}{The sequence (chromosome) name.}
  \item{chrlen}{TODO}
  \item{qname}{TODO}
  \item{fragmentProb}{TODO}
  \item{fragmentStart}{TODO}
  \item{normalMean}{TODO}
  \item{normalVar}{TODO}
  \item{numDevs}{TODO}
  \item{reads}{TODO}
}
\details{
  TODO
}
\value{
  \code{simulateReads} returns a SAM object.

  \code{realPsi} returns a numeric vector.
}
% \references{}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
\seealso{\code{\link{createGene}} for creating artificial genes, that
  are useful for simulations.
}
\examples{
## A gene with an alternative exon
mygene <- createGene(exons=list(c(1,1000), c(2001,3000), c(4001,5000)),
                     isoforms=list(c(1,2,3), c(1,3)))
reads <- simulateReads(mygene, expression=c(2/3,1/3), noReads=200,
                       readLength=40)

realPsi(mygene, reads=reads)

## Simple expression level estimator
solveIso(mygene, reads=reads)$expression
}
\keyword{misc}
