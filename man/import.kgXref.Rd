\name{import.kgXref}
\alias{import.kgXref}
\title{import the kgXref file downloaded manually from UCSC genome browser.}
\usage{
  import.kgXref(file = "/biomirror/hgdownload.cse.ucsc.edu/tables/kgXref.gz")
}
\arguments{
  \item{file}{path to kgXref file}
}
\value{
  see import.ucsc.file
}
\description{
  By default it looks in the
  /biomirror/hgdownload.cse.ucsc.edu/tables for a file
  called kgXref.gz #kgID mRNA spID spDisplayID geneSymbol
  <snip> AB001750 AB001750 O35551 RABE1_MOUSE AB001750
  <snip> AB001926 AB001926 P48722-2 P48722-2 AB001926
  <snip> AB003147 AB003147 O54770 O54770_MOUSE AB003147
  <snip> AB003503 AB003503 O88180 O88180_MOUSE AB003503
  <snip> <snip>
}
\author{
  Mark Cowley, 8 August 2006
}

