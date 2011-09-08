\name{GFF3 objects}
\alias{gff3-class}
\alias{readGFF3}
\alias{writeGFF3}
\alias{noGenes}
\alias{seqIds}
\alias{geneIds}
\alias{selectGenes}
\alias{noIso}
\alias{noExons}
\alias{getIso}
\alias{isGFF3}
\alias{geneTypes}
\alias{geneLength}
\alias{isoLength}
\alias{totalExonLength}
\title{GFF3 objects, containing gene structures}
\usage{
isGFF3 (object)

noGenes (gff3)
noIso (gff3)
noExons (gff3)

geneLength (gff3)
isoLength (gff3, collapse=FALSE)
totalExonLength (gff3, overlap=TRUE)

seqIds (gff3)
geneIds (gff3)
geneTypes (gff3)
getIso (gff3)

selectGenes (gff3, idx)

readGFF3 (file, sep = "\t", col.names = c("seqid", "source", 
     "type", "start", "end", "score", "strand", "phase", "attributes"),
     stringsAsFactors = FALSE, ...)
writeGFF3 (gff3, file, sep = "\t", quote = FALSE, row.names = FALSE,
     col.names = FALSE, ...)
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
  \item{sep}{Column separator, defaults to a single tabulator.}
  \item{col.names}{The column names to assign to the GFF data, when
    reading it from a file. Do not change the default, unless you
    really know what you are doing.}
  \item{stringsAsFactors}{Whether to read the data as factors, again,
    changing the default is not recommended.}
  \item{quote}{Whether to quote fields in the output file, it is not
    suggested to change the default, which is \code{FALSE}.}
  \item{row.names}{Whether to write row names to the file, leave this on
    \code{FALSE}, otherwise the file will not conform to the GFF3
     format.}
}
\details{
  A GFF3 object, i.e. an object with class \code{gff3} contains a number
  of gene structures, from the same organism.

  Currently it is implemented as a data frame with two extra indices, to
  be able to find genes and transcripts quickly in the table.

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

  \code{selectGenes} can be used to select a subset of genes from a GFF3
  object, based on their gene ids, or numbers.

  \code{readGFF3} reads a GFF3 file and returns a GFF3 object, while
  \code{writeGFF3} does the opposite, it writes a GFF3 object to a GFF3
  file.
}
\value{
  \code{isGFF3} returns a logical scalar.

  \code{noGenes} returns a numeric scalar

  \code{noIso}, \code{noExons} and \code{geneLength} return a named
  numeric vector. 

  \code{isoLength} returns either a named list of named numeric vectors
  (if \code{collapse} is \code{FALSE}), or a named numeric vector.

  \code{totalExonLength} returns a named numeric vector, the names the
  the gene ids.

  \code{seqIds}, \code{geneIds} and \code{geneTypes} return a character
  vector.

  \code{getIso} return a named list of character vectors.

  \code{selectGenes} returns a GFF3 object.

  \code{readGFF3} returns a GFF3 object.

  \code{writeGFF3} returns whatever \code{\link{write.table}} returns.
}
% \references{}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
% \note{}
\seealso{\code{\link{gtf2gff3}} for converting GTF files to GFF3 format
  and \code{\link{gff32gtf}} for the opposite.
}
\examples{
## Create a GFF with artificial genes
genelist <- lapply(seq(20,100,by=20), function(x) {
  createGene(exons=list(c(1,1000), c(2001, 2000+x), c(4001,5000)),
             isoforms=list(c(1,2,3), c(1,3)), id=sprintf("gene-\%d", x))
})
genes <- do.call(mergeGenes, genelist)

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
