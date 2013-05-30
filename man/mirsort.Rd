\name{mirsort}
\alias{mirsort}
\title{Sort a vector of mirbase ID's}
\usage{
  mirsort(mirs)
}
\arguments{
  \item{mirs}{a character vector of mirbase ID's (either
  precursor or mature), duplicates allowed, NA's will be
  removed.}
}
\value{
  the mirs that have been sorted alpha-numerically, such
  that hsa-mir-2 comes before hsa-mir-199
}
\description{
  Sort a vector of mirbase ID's
}
\examples{
# not run
# mirsort(affymir)
}
\author{
  Mark Cowley, 2008-06-24
}

