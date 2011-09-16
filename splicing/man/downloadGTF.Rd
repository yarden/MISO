\name{Getting genes from Ensembl}
\alias{listGTFs}
\alias{downloadGTF}
\title{
  List and download gene structures from Ensembl.
}
\description{
  \code{downloadGTF} downloads gene structures from Ensembl.
}
\usage{
listEnsemblGenes(baseURL = "ftp://ftp.ensembl.org/pub/current_gtf/")
getEnsemblGenes(speciesName,
      baseURL = "ftp://ftp.ensembl.org/pub/current_gtf/",
      verbose = TRUE)
}
\arguments{
  \item{baseURL}{The URL of the Ensembl server.}
  \item{speciesName}{The scientific name of the species,
    e.g. \sQuote{Homo sapiens}. This is case insensitive.}
  \item{verbose}{Logical scalar, whether a progress bar for the
    download progress.}
}
\details{
  \code{listEnsemblGenes} lists the species for which gene structures
  are available from Ensembl.

  \code{getEnsemblGenes} downloads the gene structures for a given
  species from Ensembl.
}
\value{
  \code{listGTFs} returns a character vector.

  \code{downloadGTF} returns a GFF3 object.
}
% \references{}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
% \note{}
\seealso{\code{\link{writeGFF3}}}
\examples{
\dontrun{
## List all organisms at the Ensembl web site
listEnsemblGenes()

## Download yeast and save it as a GFF3 file
yeast <- getEnsemblGenes("Saccharomyces cerevisiae")
writeGFF(yeast, "yeast.gff3")

## Explore the downloaded genes
yeast
noGenes(yeast)
geneIds(yeast)[1:10]
}
}
\keyword{misc}
