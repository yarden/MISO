\name{GFF3 objects}
\alias{gff3-class}
\alias{GFF3}
\alias{readGFF3}
\alias{writeGFF3}
\alias{noGenes}
\alias{noGenes.gff3}
\alias{seqIds}
\alias{seqIds.gff3}
\alias{geneIds}
\alias{geneIds.gff3}
\alias{selectGenes}
\alias{selectGenes.gff3}
\alias{noIso}
\alias{noIso.gff3}
\alias{noExons}
\alias{noExons.gff3}
\alias{getIso}
\alias{getIso.gff3}
\alias{isGFF3}
\alias{geneTypes}
\alias{geneTypes.gff3}
\alias{geneLength}
\alias{geneLength.gff3}
\alias{isoLength}
\alias{isoLength.gff3}
\alias{totalExonLength}
\alias{totalExonLength.gff3}
\alias{getSpecies}
\alias{getSpecies.gff3}
\title{GFF3 objects, containing gene structures}
\description{
  Isoform information for single or multiple genes is stored in GFF3
  objects.
}
\usage{
isGFF3 (object)

\method{noGenes}{gff3} (gff3)
\method{noIso}{gff3} (gff3)
\method{noExons}{gff3} (gff3)

\method{geneLength}{gff3} (gff3)
\method{isoLength}{gff3} (gff3)
\method{totalExonLength}{gff3} (gff3, overlap=TRUE)

\method{seqIds}{gff3} (gff3)
\method{geneIds}{gff3} (gff3)
\method{geneTypes}{gff3} (gff3)
\method{getIso}{gff3} (gff3, collapse=FALSE)
\method{getSpecies}{gff3} (gff3)

\method{selectGenes}{gff3} (gff3, idx)

readGFF3 (file)
writeGFF3 (gff3, file)
}
\arguments{
  \item{object}{An R object to be inspected.}
  \item{gff3}{A GFF3 object, i.e. an object of class \code{gff3}.}
  \item{collapse}{Logical scalar, whether to return all isoforms in a
    single vector, or by gene. The default is \code{FALSE}, the results
    are returned gene by gene.}
  \item{overlap}{Logical scalar, whether to count overlapping regions
     multiple times.}
  \item{idx}{A numeric or character vector, the genes to be selected.}
  \item{file}{The file or connection the GFF data is read from or
    written to. This argument is passed to \code{read.table} or
    \code{write.table}.}
}
\details{
  A GFF3 object, i.e. an object with class \code{gff3} contains a number
  of gene structures, from the same organism.

  In this manual page, we list the various operations of GFF3 objects.

  \code{isGFF3} decides whether its argument is a GFF3 object, it checks
  the \code{class} attribute of the supplied object.

  \code{noGenes} gives the number of genes. \code{noIso} gives the
  number of isoforms, \code{noExons} the number of different exons
  for each gene.

  \code{geneLength} returns the length of the unspliced gene, in number
  of bases. This is calculated using all exons of the gene, both
  constitutive and alternative.

  \code{isoLength} gives the length of the different mature isoforms,
  for all genes in the GFF3 object.

  \code{totalExonLength} calculates the total length of exonic regions
  in the genes. Overlapping regions are by default counted multiple
  times.

  \code{seqIds} returns the sequence ID for each gene.
  
  \code{geneIds} gives the ids of the genes.

  \code{getIso} gives the ids of the transcripts, separately for each
  gene.

  \code{getSpecies} returns the scientific name of the species
  of which the genes belong to, if available.
  
  \code{selectGenes} can be used to select a subset of genes from a GFF3
  object, based on their gene ids, or numbers.

  \code{readGFF3} reads a GFF3 file and returns a GFF3 object, while
  \code{writeGFF3} does the opposite, it writes a GFF3 object to a GFF3
  file.
}
\value{
  \code{isGFF3} returns a logical scalar.

  \code{noGenes} returns a numeric scalar.

  \code{noIso}, \code{noExons} and \code{geneLength} return a named
  numeric vector. 

  \code{isoLength} returns either a named list of named numeric vectors
  (if \code{collapse} is \code{FALSE}), or a named numeric vector.

  \code{totalExonLength} returns a named numeric vector, the names the
  the gene ids.

  \code{seqIds}, \code{geneIds} and \code{geneTypes} return a character
  vector.

  \code{getIso} returns a named list of character vectors.

  \code{getSpecies} returns a character scalar.
  
  \code{selectGenes} returns a GFF3 object.

  \code{readGFF3} returns a GFF3 object.

  \code{writeGFF3} returns \code{NULL}, invisibly.
}
% \references{}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
% \note{}
% \seealso{}
\examples{
## Create a GFF with artificial genes
genelist <- lapply(seq(20,100,by=20), function(x) {
  createGene(exons=list(c(1,1000), c(2001, 2000+x), c(4001,5000)),
             isoforms=list(c(1,2,3), c(1,3)), id=sprintf("gene-\%d", x))
})
genes <- do.call(merge, genelist)

## Query some properties
noGenes(genes)
noIso(genes)
noExons(genes)
geneLength(genes)
isoLength(genes)
geneIds(genes)
geneTypes(genes)
getIso(genes)

## Create a subset
twogenes <- selectGenes(genes, c(3,1))
noGenes(twogenes)
geneIds(twogenes)

## A gene that has two exons with a common part
mygene <- createGene(exons=list(c(1,1000), c(501,1500), c(2001,3000)),
                     isoforms=list(c(1,3), c(2,3)), id="ogene")
noExons(mygene)
isoLength(mygene)
totalExonLength(mygene)
totalExonLength(mygene, overlap=FALSE)
}
\keyword{misc}
