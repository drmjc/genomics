\name{mirIDsplit}
\alias{mirIDsplit}
\title{Split the parts of a miRanda miRNA ID into the constituent parts}
\usage{
  mirIDsplit(mirs)
}
\arguments{
  \item{mirs}{a character vector of mirbase ID's (either
  precursor or mature), duplicates allowed, NA's will NOT
  be removed.}
}
\value{
  a 4 column data.frame with colnames: species, accession,
  ids, isoform
}
\description{
  Split the parts of a miRanda miRNA ID into the
  constituent parts
}
\examples{
mirIDsplit("hsa-let-7a-3")
# [1] "hsa" "let" "7a" "3"
}
\author{
  Mark Cowley, 2008-06-24
}

