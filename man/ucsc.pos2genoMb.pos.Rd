\name{ucsc.pos2genoMb.pos}
\alias{ucsc.pos2genoMb.pos}
\title{chromosome position to genome position
Useful for Manhattan style plots where markers on chr 2 should be to the
right of those on chr1}
\usage{
  ucsc.pos2genoMb.pos(x, marker.info)
}
\arguments{
  \item{x}{a vector of ucsc-style genomic locations}

  \item{marker.info}{a table of marker locations, with at
  least these columns: Mb, GenoMb, chr}
}
\value{
  vector of genomic locations
}
\description{
  chromosome position to genome position Useful for
  Manhattan style plots where markers on chr 2 should be to
  the right of those on chr1
}
\author{
  Mark Cowley, 2011-08-02
}

