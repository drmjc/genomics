\name{import.ncbi.file}
\alias{import.ncbi.file}
\title{import.ncbi.file}
\usage{
  import.ncbi.file(file, taxid, colnames, sortby.col)
}
\arguments{
  \item{file}{path to ncbi-formatted file}

  \item{taxid}{numeric taxon id}

  \item{colnames}{character vector}

  \item{sortby.col}{undocumented}
}
\value{
  \code{data.frame}
}
\description{
  General wrapper to import and NCBI file that has usually
  been mirrored to the local disc. It will copy the
  compressed or uncompressed file, unzip it in the tmp dir,
  grep out the taxon of interest, resort by a given column,
  and return the data
}
\author{
  Mark Cowley, 10 April 2006
}

