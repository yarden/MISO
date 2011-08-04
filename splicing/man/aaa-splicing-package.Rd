\name{splicing-package}
\alias{splicing-package}
\alias{splicing}
\docType{package}
\title{
  Package to study alternative splicing.
}
\description{
  This package contains functions to quantitate alternatively spliced
  isoforms, right now with a simple linear deconvolution. It also has
  additional functions to create artificial gene structures, simulate
  reads from real or artificial genes, calculate the complexity of gene
  structures and convert between GTF and GFF3 files.
}
% \details{}
\author{ Gabor Csardi \email{gcsardi@stat.harvard.edu} }
\references{
  The \code{\link{listGTFs}} and \code{\link{downloadGTF}} functions can
  be used to download gene information from Ensembl, for various
  species.

  Gene information can be stored in various formats, the GTF and GFF3
  formats are the ones that are most commonly used. 
  The \code{\link{gtf2gff3}} and \code{\link{gff32gtf}} functions
  convert between the GTF and GFF3 format.

  \code{\link{readGFF3}}
  \code{\link{writeGFF3}}
  \code{\link{writeGTF}}

  \code{\link{noGenes}}
  \code{\link{geneIds}}
  \code{\link{selectGenes}}
  \code{\link{noIso}}
  \code{\link{noExons}}
  \code{\link{getIso}}
  \code{\link{isGFF3}}
  \code{\link{geneTypes}}
  \code{\link{geneLength}}
  \code{\link{isoLength}}

  \code{\link{createGene}}
  \code{\link{mergeGenes}}
  
  \code{\link{geneReadsForGene}}
  \code{\link{genReads}}
  
  \code{\link{writeSAM}}
  \code{\link{sam2bam}}

  \code{\link{condIso}}

  \code{\link{miso}}

  \code{\link{plotIso}}
  \code{\link{plotIsoSize}}
  \code{\link{plotIsoPFD}}

  \code{\link{matchIso}}

  \code{\link{condIso}}
  
  \code{\link{solveIso}}
}
\keyword{package}
% \seealso{}
% \examples{}
