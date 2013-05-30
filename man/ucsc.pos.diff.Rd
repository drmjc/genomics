\name{ucsc.pos.diff}
\alias{ucsc.pos.diff}
\title{Method to determine the distance between two genomic regions in "ucsc.pos"
format.}
\usage{
  ucsc.pos.diff(from, to,
    method = c("min", "start", "stop", "middle"))
}
\value{
  an integer, or Inf
}
\description{
  Method to determine the distance between two genomic
  regions in "ucsc.pos" format.
}
\examples{
ucsc.pos.diff("chr9:1000-2000", "chr9:1500-2100", "min")
# [1] 0
ucsc.pos.diff("chr9:1000-2000", "chr9:1500-2100", "start")
# [1] 500
ucsc.pos.diff("chr9:1000-2000", "chr9:1500-2100", "stop")
# [1] 100
ucsc.pos.diff("chr9:1000-2000", "chr9:1500-2100", "middle")
# [1] 300
ucsc.pos.diff("chr9:1000-2000", "chr9:2100-2200", "min")
# [1] 100
ucsc.pos.diff("chr9:1000-2000", "chr9:2100-2200", "start")
# [1] 1100
ucsc.pos.diff("chr9:1000-2000", "chr9:2100-2200", "stop")
# [1] 200
ucsc.pos.diff("chr9:1000-2000", "chrX:2100-2200", "stop")
# [1] Inf
}
\author{
  Mark Cowley, 18/4/07
}

