\name{order.refseq}
\alias{order.refseq}
\title{order RefSeq ID's by validation status}
\usage{
  order.refseq(x, na.last = TRUE, decreasing = FALSE)
}
\arguments{
  \item{x}{undocumented}

  \item{na.last}{see order}

  \item{decreasing}{see order}
}
\value{
  Undocumented return value
}
\description{
  return the order of a vector of RefSeq status's such that
  those that are VALIDATED come before those that are
  PREDICTED.
}
\details{
  Here are the status rankings: VALIDATED >> REVIEWED >>
  PROVISIONAL >> INFERRED >> PREDICTED = GENOME ANNOTATION
  = WGS 1 2 3 4 5 5 5 see ?order for equivalent method
  using normal integers.
}
\author{
  Mark Cowley, 10 April 2006
}

