\name{MISO}
\alias{MISO}
\alias{runMISO}
\alias{readMISO}
\alias{writeMISO}
\title{Mixture of Isoforms (MISO)}
\description{
  TODO
}
\usage{
MISO (geneStructure, gene = 1L, reads, readLength = getReadLength(reads),
    overHang = 1L, noIterations = 5000, noBurnIn = 500, noLag = 10,
    hyperparameters = rep(1, noIso(geneStructure)[gene]),
    paired = isPaired(reads), fragmentProb = NULL, 
    fragmentStart = 0L, normalMean, normalVar, numDevs)

runMISO(geneStructure, readsfile, results = c("return", "files",
    "Rfiles"), resultDir = ".", readLength = NULL, verbose = TRUE,
    snowCluster = NULL, seed = NULL, \dots)

writeMISO(misoResult, file)
readMISO(file)
}
\arguments{
  \item{geneStructure}{TODO}
  \item{gene}{TODO}
  \item{reads}{TODO}
  \item{readLength}{TODO}
  \item{overHang}{TODO}
  \item{noIterations}{TODO}
  \item{noBurnIn}{TODO}
  \item{noLag}{TODO}
  \item{hyperparameters}{TODO}
  \item{paired}{TODO}
  \item{fragmentProb}{TODO}
  \item{fragmentStart}{TODO}
  \item{normalMean}{TODO}
  \item{normalVar}{TODO}
  \item{numDevs}{TODO}
  \item{readsfile}{TODO}
  \item{results}{TODO}
  \item{resultDir}{TODO}
  \item{verbose}{TODO}
  \item{snowCluster}{TODO}
  \item{seed}{TODO}
  \item{\dots}{TODO}
  \item{misoResult}{TODO}
  \item{file}{TODO}
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
