\name{Read objects}
\alias{SAM}
\alias{splicingSAM}
\alias{getAttribute}
\alias{getAttribute.splicingSAM}
\alias{getAttributes}
\alias{getAttributes.splicingSAM}
\alias{getIsoform}
\alias{getIsoform.splicingSAM}
\alias{getReadLength}
\alias{getReadLength.splicingSAM}
\alias{getStrand}
\alias{getStrand.splicingSAM}
\alias{isMixed}
\alias{isMixed.splicingSAM}
\alias{isPaired}
\alias{isPaired.splicingSAM}
\alias{listAttributes}
\alias{listAttributes.splicingSAM}
\alias{noPairs}
\alias{noPairs.splicingSAM}
\alias{noReads}
\alias{noReads.splicingSAM}
\alias{noSingles}
\alias{noSingles.splicingSAM}
\alias{readSAM}
\alias{selectReads}
\alias{selectReads.splicingSAM}
\alias{seqNames}
\alias{seqNames.splicingSAM}
\alias{toTable}
\alias{toTable.splicingSAM}
\alias{writeSAM}
\alias{writeSAM.splicingSAM}
\alias{mergeReads}
\title{A set of short reads, corresponding to a SAM or BAM file}
\description{
  TODO
}
\usage{
readSAM(filename, region = NULL)
\method{writeSAM}{splicingSAM} (reads, conn)
\method{toTable}{splicingSAM} (reads, \dots)

\method{selectReads}{splicingSAM} (reads, idx)
mergeReads(\dots)

\method{seqNames}{splicingSAM} (reads)
\method{getReadLength}{splicingSAM} (reads)
\method{getStrand}{splicingSAM} (reads)

\method{noReads}{splicingSAM} (reads)
\method{noPairs}{splicingSAM} (reads)
\method{noSingles}{splicingSAM} (reads)
\method{isMixed}{splicingSAM} (reads)
\method{isPaired}{splicingSAM} (reads)

\method{getAttributes}{splicingSAM} (reads)
\method{listAttributes}{splicingSAM} (reads)
\method{getAttribute}{splicingSAM} (reads, attr)
\method{getIsoform}{splicingSAM} (reads)
}
\arguments{
  \item{filename}{TODO}
  \item{region}{TODO}
  \item{reads}{TODO}
  \item{conn}{TODO}
  \item{\dots}{TODO}
  \item{idx}{TODO}
  \item{attr}{TODO}
}
\details{
  TODO
}
\value{
  TODO
}
% \references{TODO}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
\seealso{\code{\link{filterReads}} for omitting some of the reads.}
\examples{
## TODO
}
\keyword{misc}
