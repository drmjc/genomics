\name{range_merge_markers}
\alias{range_merge_markers}
\title{merge a range of markers}
\usage{
  range_merge_markers(x, marker.names, gaplength = 0,
    sep = "-", single.string = ", ")
}
\arguments{
  \item{x}{1:length(marker.names)}

  \item{marker.names}{The marker names (or words) to look
  up the indices supplied in x.}

  \item{gaplength}{allow N gaps between consecutive markers
  to be allowed. 0 implies no gaps; 1 implies if markers A
  and C are significant, then marker B is expected to be
  significant, and included in the range.}

  \item{sep}{Which character should be used to delimit a
  range of words?}

  \item{single.string}{if NULL, then a character vector is
  returned. If not NULL then it's value will be used to
  concatenate the ranges into a single string, thus ", "
  comma separates the ranges. For clarity, avoid using the
  same symbol as specified by sep.}
}
\value{
  A character vector with ranges identified by the first
  and last marker, and singletons identified by their name.
}
\description{
  Usful function to annotate the markers that are
  significant for a particular gene in a genome scan. First
  work out which markers are significant (either as
  indices, or as the marker names) then supply these
  markers to this function, along with the names of all
  markers, and you will be supplied with a char vector of
  these markers, and if a few markers in a row are
  significant, then just the start and end marker with a
  "-" (or sep) is returned.
}
\examples{
range_merge_markers(c(1:4,6), LETTERS)
# [1] "A-D" "F"
range_merge_markers(c(LETTERS[1:4], LETTERS[6]), LETTERS)
# [1] "A-D" "F"
}
\author{
  Mark Cowley, 22 march 2006
}

