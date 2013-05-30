\name{alignment2ucsc}
\alias{alignment2ucsc}
\title{Convert an alignment to a ucsc location.}
\usage{
  alignment2ucsc(x)
}
\arguments{
  \item{x}{an alignment object (see Details)}
}
\value{
  a ucsc string, eg \dQuote{chr13:40727709-40745894}
}
\description{
  alignment looks like: name chrom chromStart chromEnd
  strand 1 NM_011547 chr13 40727709 40745894 - 2 NM_010791
  chr11 101693600 101710444 - 3 NM_013684 chr17 15204872
  15222390 +
}
\author{
  Mark Cowley
}

