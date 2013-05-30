\name{miRNA.fix.miR.ids}
\alias{miRNA.fix.miR.ids}
\title{Fix the names of microRNA's}
\usage{
  miRNA.fix.miR.ids(miRids)
}
\arguments{
  \item{miRids}{a vector of miR ID's that need checking for
  validity.}
}
\value{
  a vector of miR ID's that have been fixed according to a
  set of simple rules.
}
\description{
  If there are any systematic naming problems with miRs,
  this method can fix them. Currently implemented rules:\cr
  replace -star -> *
}
\author{
  Mark Cowley, 2009-01-13
}

