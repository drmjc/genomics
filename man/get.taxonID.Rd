\name{get.taxonID}
\alias{get.taxonID}
\title{get.taxonID}
\usage{
  get.taxonID(taxid)
}
\arguments{
  \item{taxid}{species name, some aliases, or 2 letter
  code}
}
\value{
  numeric taxon ID
}
\description{
  return taxon ID for a species name or short code
}
\examples{
get.taxonID("mouse")
get.taxonID("mm")
get.taxonID("hs")
get.taxonID("rat")
get.taxonID("rn")
}
\author{
  Mark Cowley, 2013-05-30
}

