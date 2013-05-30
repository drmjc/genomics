\name{mirna.pre2mature}
\alias{mirna.pre2mature}
\title{Convert pre to mature microRNA's
Convert a list of precursor mir ID's to all possible mature
sequences/synonyms}
\usage{
  mirna.pre2mature(mirs, patterns = c("*", "-5p", "-3p"),
    major.only = FALSE)
}
\arguments{
  \item{mirs}{a character vector of precursor mirbase ID's}

  \item{patterns}{the possible patterns that make up a
  mature miRNA ID}

  \item{major.only}{return only the major form of the
  mature miRNA? (ie only 1 value)}
}
\value{
  a list with N elements, one for each mir in 'mirs'. Each
  element is a character vector of possible miRNA ID's.
}
\description{
  Convert pre to mature microRNA's Convert a list of
  precursor mir ID's to all possible mature
  sequences/synonyms
}
\author{
  Mark Cowley, isoD
}

