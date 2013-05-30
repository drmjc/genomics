\name{import.mirbase.dropped}
\alias{import.mirbase.dropped}
\title{Import the miRNA.dead file
Import the miRNA.dead file, and return the mir ID's that have been dropped.}
\usage{
  import.mirbase.dropped(file, species = "hsa")
}
\arguments{
  \item{file}{the path to \dQuote{miRNA.dead.gz}. It can be
  unzipped or gzipped}

  \item{species}{\code{NULL} means import all species,
  otherwise use a 3 letter species code, eg \dQuote{hsa}}
}
\description{
  Import the miRNA.dead file Import the miRNA.dead file,
  and return the mir ID's that have been dropped.
}
\author{
  Mark Cowley, 2009-01-13
}

