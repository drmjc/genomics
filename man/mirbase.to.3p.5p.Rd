\name{mirbase.to.3p.5p}
\alias{mirbase.to.3p.5p}
\title{Convert miR ID's to 3p/5p
The miRNA.xls file from miRbase uses the notation Mature1 and Mature2.}
\usage{
  mirbase.to.3p.5p(mirbase)
}
\arguments{
  \item{mirbase}{see \code{\link{import.miRbase}}. You must
  specify \code{predict.arm=TRUE}}
}
\value{
  a data.frame same dim as mirbase, but with mature1 and 2
  changed for 5p and 3p, and miR data shuttled into the
  correct columns depending on whether it's on the 5' or 3'
  arm. The Mature1_Arm/Mature2_Arm columns are dropped
  though.
}
\description{
  If a mir has 2 miRs then the Mature 1 and 2 are always
  the 5' and 3' respectively. BUT if there's only one known
  miR for a mir, then you don't know if it's a 5' or 3'.
  This function converts the mirbase file from Mature1/2 to
  Mature5p/3p
}
\author{
  Mark Cowley, 2009-01-09
}

