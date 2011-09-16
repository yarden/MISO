\name{splicing-package}
\alias{splicing-package}
\alias{splicing}
\docType{package}
\title{
  Package to study alternative splicing.
}
\description{
  This package contains functions to quantitate alternatively spliced
  isoforms. It implements a simple linear deconvolution and the Mixture
  of Isoforms (MISO) method. It also has additional functions to create
  artificial gene structures, simulate reads from real or artificial
  genes, calculate the complexity of gene structures. It support simple
  visualization of gene structures.
}
\details{
  We list the most important functions of the splicing package
  here. Please see their manual pages for the details and additional
  functions.
  
  The \code{\link{listEnsemblGenes}} and \code{\link{getEnsemblGenes}}
  functions can be used to download gene information from Ensembl, for
  various species. 

  The package uses GFF3 objects for storing gene structures.
  \code{\link{readGFF3}} reads a GFF3 file, and \code{\link{writeGFF3}}
  stored gene structure information in a GFF3 files.

  GFF3 objects have the following methods:
  \itemize{
    \item \code{\link{isGFF3}}
    \item \code{\link{writeGFF3}}
    \item \code{\link{noGenes}}
    \item \code{\link{geneIds}}
    \item \code{\link{selectGenes}}
    \item \code{\link{noIso}}
    \item \code{\link{noExons}}
    \item \code{\link{getIso}}
    \item \code{\link{isGFF3}}
    \item \code{\link{geneTypes}}
    \item \code{\link{geneLength}}
    \item \code{\link{isoLength}}
    \item \code{\link{totalExonLength}}
    \item \code{\link{getSpecies}}
    \item \code{\link{mergeGFF3}}
  }
  
  The \code{\link{createGene}} function creates (possibly artificial)
  gene structures.
  
  \code{simulateReads} simulates short reads for a given gene
  structure. \code{\link{realPsi}} gives the expression levels of the
  various isoforms for a gene with simulated reads.

  \code{\link{rpk}} calculates the reads-per-kilobase measure for a gene
  and a set of reads.

  SAM objects are used for storing short reads. Their operations:
  \itemize{
    \item \code{\link{noReads}}
    \item \code{\link{isPaired}}
    \item \code{\link{isMixed}}
    \item \code{\link{noPairs}}
    \item \code{\link{noSingles}}
    \item \code{\link{toTable}}
    \item \code{\link{readLength}}
    \item \code{\link{seqNames}}
    \item \code{\link{selectReads}}
    \item \code{\link{filterReads}}
    \item \code{\link{getStrand}}
    \item \code{\link{getAttributes}}
    \item \code{\link{listAttributes}}
    \item \code{\link{getAttribute}}
    \item \code{\link{getIsoform}}
    \item \code{\link{writeSAM}}
  }
  
  \code{\link{condIso}} calculates the complexity of a gene structure.

  \code{\link{MISO}} runs the Mixture of Isoforms algorithm for a given
  gene and set of reads, and returns a MISOresult object.

  For visualizing gene structures, see the \code{\link{plotIso}},
  \code{\link{plotIsoSize}} and \code{\link{plotIsoPDF}} functions.

  \code{\link{matchIso}} matches a set of reads to all isoforms of a
  gene. 

  \code{\link{solveIso}} quantitates the isoform expression levels for a
  gene, using a simple (and fast!) linear deconvolution algorithm.
}
\references{
  Katz Y, Wang ET, Airoldi EM, Burge CB: Analysis and design of RNA
  sequencing experiments for identifying isoform regulation.
  Nat Methods. 2010 Dec ; 7(12): 1009-15
}
\author{ Gabor Csardi \email{gcsardi@stat.harvard.edu} }
\keyword{package}
% \seealso{}
% \examples{}
