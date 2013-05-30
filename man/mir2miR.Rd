\name{miR2mir}
\alias{miR2mir}
\title{Convert mature miR ID's into the immature, precursor form.}
\usage{
  miR2mir(miRs, patterns = c("\\*", "-5p", "-3p"))
}
\arguments{
  \item{miRs}{a character vector of mature miR ID's}

  \item{patterns}{the possible patterns on the end of the
  miR ID}
}
\value{
  a character vector of the mir ID corresponding to each
  miR ID, with the caveat that there may be multiple mir
  genes corresponding to each mature miR
}
\description{
  NB, this has no way of knowing which mir gene if there
  are more than one...
}
\author{
  Mark Cowley, 2008-06-25
}

