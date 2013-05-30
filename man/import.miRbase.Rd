\name{import.miRbase}
\alias{import.miRbase}
\title{Import miRBase records.
From the miRBase FTP site ([1]), there is a file called
\dQuote{miRNA.xls}, which contains lots of the useful information about each miRNA.}
\usage{
  import.miRbase(miRNA.xls.file, species.code = "hsa",
    gff = NULL, predict.arm = TRUE)
}
\arguments{
  \item{miRNA.xls.file}{the path to the miRNA.xls file}

  \item{species.code}{3 letter code such as \dQuote{hsa},
  \dQuote{mmu}, \dQuote{rno}}

  \item{gff}{if \code{NULL} then method will look for a
  file called \dQuote{genomes/<species code>.gff} in the
  same directory as the \code{miRNA.xls.file}; if
  \code{NULL}, then no location info will be inferred;
  otherwise provide the path to the \dQuote{<species
  code>.gff} file}

  \item{predict.arm}{logical: If \code{TRUE}, then predict
  which arm the mature miR exists on (3'/5'). Otherwise set
  to \code{FALSE}.}
}
\description{
  This function imports the \code{miRNA.xls} file, subsets
  the mir's to those from a particular species, and
  re-sorts the resulting file. It can also obtain the
  locations of the miRNA's, either automatically (ie by
  mirroring the entire \sQuote{sequences/CURRENT} directory
  from above FTP site), or by specifying the gff file.
}
\author{
  Mark Cowley, 2009-01-07
}
\references{
  [1]
  \url{ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT}
}

