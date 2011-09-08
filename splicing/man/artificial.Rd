\name{Artificial genes}
\alias{createGene}
\alias{mergeGenes}
\title{Create artificial gene structures for simulation purposes}
\usage{
createGene (exons, isoforms, id = "insilicogene",
     seqid = "seq1", source = "protein_coding", strand = "+")
mergeGenes (\dots, mode=c("disjunct"))
}
\arguments{
  \item{exons}{A list of numeric vectors, each of length two. This list
    gives all constitutive and alternative exons, and an exon is given
    by its start and end coordinates (inclusive). The exons are
    numbered according to their order in this list.}
  \item{isoforms}{A list of isoforms, each isoform is given by a numeric
    vector, a list of exons, where the numbers refer to the order of the
    exons in the \code{exons} argument.}
  \item{id}{A character scalar, the id of the gene.}
  \item{seqid}{Character scalar, the sequence id (i.e. usually
    chromosome name) of the gene.}
  \item{source}{Character scalar, it corresponds to the \code{source}
    field in the GFF3 format.}
  \item{strand}{Character scalar, the strand of the gene, usually
    \code{+} or \code{-}.}
  \item{\dots}{The genes to merge. If you have your genes in a list,
    then you can use \code{do.call}, see the example below.}
  \item{mode}{Character scalar, defines how the genes are merged. Right
    now only \code{disjunct} is supported, which means that the genes
    are artificially made disjunct first, by increasing the coordinates
    of all exons. Obviously, this is useful for artificial genes only.}
}
\details{
  \code{createGene} creates an artificial gene structure, that can be
  used in simulations, when studying properties of methods estimating
  isoform expression.

  \code{mergeGenes} merges gene structures from two or more GFF objects.
}
\value{
  A GFF3 object.
}
% \references{}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
% \note{}
\seealso{
  \code{\link{gff3-class}} for the description of the GFF3 class.
}
\examples{
## A simple gene, with a skipped exon
mygene <- createGene (exons=list(c(1,1000), c(2001,3000), c(4001,5000)),
                      isoforms=list(c(1,2,3), c(1,3)))
noIso(mygene)
geneIds(mygene)
getIso(mygene)

## Create ten genes, with an increasingly bigger skipped exon,
## and merge them
genelist <- lapply(seq(10,100,by=10), function(x) {
  createGene(exons=list(c(1,1000), c(2001, 2000+x), c(4001,5000)),
             isoforms=list(c(1,2,3), c(1,3)), id=sprintf("gene-\%d", x))
})
genes <- do.call(mergeGenes, genelist)
noIso(genes)
geneIds(genes)
isoLength(genes)
}
\keyword{misc}
