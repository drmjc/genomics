\name{import.knowngene}
\alias{import.knowngene}
\title{import the knownGene file downloaded manually from UCSC genome browser.}
\usage{
  import.knowngene(file = "/biomirror/hgdownload.cse.ucsc.edu/tables/knownGene")
}
\arguments{
  \item{file}{path to knownGene file}
}
\value{
  see import.ucsc.file
}
\description{
  By default it looks in the
  /biomirror/hgdownload.cse.ucsc.edu/tables for a file
  called knownGene #name chrom strand txStart txEnd
  cdsStart cdsEnd <snip> NM_001011874 chr1 - 3174829
  3630690 3176369 3630540 <snip> AK035064 chr1 + 3709446
  3783878 3709490 3783206 <snip> NM_011283 chr1 - 4304522
  4320770 4304977 4313203 <snip> NM_011441 chr1 - 4451317
  4456803 4452105 4453796 <snip> <snip>
}
\author{
  Mark Cowley, 8 August 2006
}

