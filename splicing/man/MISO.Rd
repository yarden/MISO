\name{MISO}
\alias{MISO}
\title{Mixture of Isoforms (MISO)}
\description{
  TODO
}
\usage{
MISO (geneStructure, gene = 1L, reads, readLength = readLength(reads),
    overHang = 1L, noIterations = 5000, noBurnIn = 500, noLag = 10,
    hyperparameters = rep(1, noIso(geneStructure)[gene]),
    paired = reads$paired, fragmentProb = NULL, 
    fragmentStart = 0L, normalMean, normalVar, numDevs) 
}
\arguments{
  \item{geneStructure}{}
  \item{gene}{}
  \item{reads}{}
  \item{readLength}{}
  \item{overHang}{}
  \item{noIterations}{}
  \item{noBurnIn}{}
  \item{noLag}{}
  \item{hyperparameters}{}
  \item{paired}{}
  \item{fragmentProb}{}
  \item{fragmentStart}{}
  \item{normalMean}{}
  \item{normalVar}{}
  \item{numDevs}{}
}
\details{
  TODO
}
\value{
  A \code{MISOresult} object.
}
\references{TODO}
\author{TODO}
\seealso{TODO}
\examples{
## TODO
}
\keyword{misc}
