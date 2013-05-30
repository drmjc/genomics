\name{mirorder}
\alias{mirorder}
\title{Ordering permutation of mirbase ID's
Get the index that each ID within a vector of mirbase ID's would be at if \code{mirs}
was sorted/}
\usage{
  mirorder(mirs)
}
\arguments{
  \item{mirs}{a character vector of mirbase ID's (either
  precursor or mature), duplicates allowed, NA's will be
  removed.}
}
\value{
  the numeric index that each mir \emph{would} be in if
  mirs was sorted.
}
\description{
  Ordering permutation of mirbase ID's Get the index that
  each ID within a vector of mirbase ID's would be at if
  \code{mirs} was sorted/
}
\examples{
\dontrun{
# mirsort(affymir)
}
}
\author{
  Mark Cowley, 2008-06-24
}

