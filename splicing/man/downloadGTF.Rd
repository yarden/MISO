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
listGTFs(baseURL = "ftp://ftp.ensembl.org/pub/current_gtf/")
downloadGTF(speciesName,
      baseURL = "ftp://ftp.ensembl.org/pub/current_gtf/",
      verbose = FALSE)
}
\arguments{
  \item{baseURL}{The URL of the Ensembl server.}
  \item{speciesName}{The scientific name of the species,
    e.g. \sQuote{Homo sapiens}. This is case insensitive.}
  \item{verbose}{Logical scalar, whether to show the download progress.}
}
% \details{}
\value{
  \code{listGTFs} returns a character vector.

  \code{downloadGTF} returns a data frame, that can be converted to a
  GFF3 object via \code{\link{gtf2gff2}}.
}
% \references{}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
% \note{}
\seealso{\code{\link{gtf2gff3}}}
\examples{
## List all organisms at the Ensembl web site
\dontrun{listGTFs()}

## Download yeast and convert it to GFF3
yeast <- downloadGTF("Saccharomyces cerevisiae")
yeastGFF <- gtf2gff3(yeast)
noGenes(yeastGFF)
geneIds(yeastGFF)[1:10]
}
\keyword{misc}
