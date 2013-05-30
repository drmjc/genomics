\name{import.unigene2gene}
\alias{import.unigene2gene}
\title{Import the gene2unigene file.}
\usage{
  import.unigene2gene(taxid = TAXID.MOUSE,
    file = "/biomirror/ftp.ncbi.nih.gov/gene/DATA/gene2unigene")
}
\arguments{
  \item{taxid}{numeric taxon id}

  \item{file}{path to gene2unigene file}
}
\value{
  \code{data.frame}
}
\description{
  Import the gene2unigene file.
}
\author{
  Mark Cowley, 30/3/07
}

