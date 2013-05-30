\name{import.knownGene2gene}
\alias{import.knownGene2gene}
\title{Download knownGeneToLocusLink from UCSC Table Browser to: Import the
knownToLocusLink file, downloaded from the UCSC Table Browser.}
\usage{
  import.knownGene2gene(taxid = TAXID.MOUSE,
    files = list(mouse = "/biomirror/ucsc/tables/mm8/knownToLocusLink.gz", human = "/biomirror/ucsc/tables/hg18/knownToLocusLink.gz", rat = "/biomirror/ucsc/tables/rn4/knownToLocusLink.gz"))
}
\arguments{
  \item{taxid}{the numerical taxonID}

  \item{files}{a list with your file name in it.}
}
\description{
  Download knownGeneToLocusLink from UCSC Table Browser to:
  Import the knownToLocusLink file, downloaded from the
  UCSC Table Browser.
}
\author{
  Mark Cowley, 30/3/07
}

