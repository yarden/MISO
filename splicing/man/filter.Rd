\name{Filtering reads}
\alias{filterReads}
\alias{FILTER}
\title{Filtering reads based on various criteria}
\description{
  TODO
}
\usage{
\method{filterReads}{splicingSAM}(reads, \dots)
FILTER(name, \dots)
}
\arguments{
  \item{reads}{The reads to filter, a \code{SAM} object.}
  \item{\dots}{For \code{filterReads} the filters themselves, see
    details below.

    For \code{FILTER} the parameters of the filter.
  } 
  \item{name}{The name of the filter, see implemented filters below.}
}
\details{
  TODO
}
\value{
  \code{filterReads} returns the filtered SAM object.
}
% \references{}
\author{
  Gabor Csardi \email{gcsardi@stat.harvard.edu}
}
% \seealso{}
\examples{
## TODO
}
\keyword{misc}
