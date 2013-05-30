\name{mirsplit}
\alias{mirsplit}
\title{Split the parts of a mirana miRNA ID into the consituent parts}
\usage{
  mirsplit(mirs)
}
\arguments{
  \item{mirs}{a character vector of mirbase ID's (either
  precursor or mature), duplicates allowed, NA's will NOT
  be removed.}
}
\value{
  a 4 column \code{data.frame} with colnames:
  \dQuote{species}, \dQuote{accession}, \dQuote{ids},
  \dQuote{isoform}
}
\description{
  Split the parts of a mirana miRNA ID into the consituent
  parts
}
\examples{
mirIDsplit("hsa-let-7a-3")
# species   accession   ids  isoform
#   "hsa"       "let"  "7a"      "3"
}
\author{
  Mark Cowley, 2008-06-24
}

