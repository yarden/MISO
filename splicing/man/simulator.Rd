\name{Read simulator}
\alias{genReadsForGene}
\alias{genReads}
\title{Simulating RNASeq reads from real or artifical genes}
\usage{
genReadsForGene (geneStructure, expression,
     gene = geneIds(geneStructure)[1], coverage = 10, readLength = 33,
     paired = FALSE, insertDistr = NULL, normalInsert = c(250, 15))
genReads (geneStructure, expression, genes = geneIds (geneStructure))
}
\arguments{
  \item{geneStructure}{A GFF3 object containing one of more genes.}
  \item{expression}{A numeric vector for \code{genReadsForGene}, a list
    of numeric vectors for \code{genReads}. The list can be named, then
    the names should correspond to gene IDs in \code{geneStructure}. If
    the list is not named, then it is assumed to have the same order as
    the genes in \code{geneStructure}.}
  \item{gene}{If generating reads for a single gene with
    \code{genReadsForGene}, then this argument should be a single gene
    ID or a number, that defines the gene for which reads are
    generated. By default the first gene is used.}
  \item{coverage}{Coverage, a numeric scalar. This specifies the number
    of reads to generate. It is calculated as the number of read that
    cover each base the given number of times, on average. }
  \item{readLength}{Numeric scalar, the read length.}
  \item{paired}{Logical scalar, whether to generate paired reads.}
  \item{insertDistr}{Numeric vector, the distribution of insert
    length. The first value is the probability of having insert length
    one, etc. If this vector is not given or \code{NULL}, then
    \code{normalInsert} is used to calculate the insert length
    distribution. This argument is ignored if \code{paired} is
    \code{FALSE}.}
  \item{normalInsert}{A numeric vector of length two. This is an
    alternative to the explicit \code{insertDistr} form. The first value
    is the mean, the second the standard deviation of the (distretized)
    normal distribution that is used as the insert length distribution.
    This argument is ignored if \code{paired} is FALSE, or
    \code{insertDistr} is given (and not \code{NULL}).}
  \item{genes}{The numbers or ids of the genes for which reads are
    generated. By default all genes are used.}
}
\details{
  TODO
}
\value{
  Both functions return a SAM object.
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
reads <- genReadsForGene(mygene, expression=c(2/3,1/3))

## Simple expression level estimator
solveIso(mygene, reads)
}
\keyword{misc}
