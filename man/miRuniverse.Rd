\name{miRuniverse}
\alias{miRuniverse}
\title{Obtain the miR universe
What are the known miR's in a given mirbase table.}
\usage{
  miRuniverse(mirbase)
}
\arguments{
  \item{mirbase}{a mirbase table}
}
\value{
  a vector of all known miR's
}
\description{
  the mirbase table can be native (contains Mature1 or
  Mature2 columns), or updated to Mature5p/Mature3p style,
  otherwise it can be a miRbase where the miR's are in the
  ID column
}
\author{
  Mark Cowley, 2009-01-09
}

