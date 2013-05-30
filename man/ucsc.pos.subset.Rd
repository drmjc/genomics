\name{ucsc.pos.subset}
\alias{ucsc.pos.subset}
\title{ucsc.pos.subset}
\usage{
  ucsc.pos.subset(intervals, start, stop, strand)
}
\arguments{
  \item{intervals}{vector of UCSC style intevals, eg
  "chr1:100-300"}

  \item{start}{vector of numerics defining the starting
  base within each interval (1-based)}

  \item{stop}{vector of numerics defining the ending base
  within each interval (1-based). NB stop can be longer
  than the original interval.}

  \item{strand}{either a charcter vector of "+"/"-", or
  {0,1}, or {-1,1}}
}
\value{
  a vector of UCSC style intervals, adjusted as detailed
  above
}
\description{
  For a set of UCSC style intervals, produce a new set of
  intervals that are subsets of those intervals. Strand
  must be provided, because intervals are strand-less, but
  the logic to subset an interval depends very much on the
  strand.
}
\details{
  eg: subset chr1:1000-1100 (+) to the interval that starts
  at 2-7 -> chr1:1001-1006 eg: subset chr1:1000-1100 (-) to
  the interval that starts at 2-7 -> chr1:1094-1099
}
\examples{
ucsc.pos.subset("chr1:1000-1100", 1,11, "+")
# [1] "chr1:1000-1010"
ucsc.pos.subset("chr1:1000-1100", 1,11, "-")
# [1] "chr1:1090-1100"
# (you can use > 1 entry too.)
}
\author{
  Mark Cowley, 2009-01-08
}

