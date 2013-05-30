\name{import.ucsc.file}
\alias{import.ucsc.file}
\title{import UCSC-style file}
\usage{
  import.ucsc.file(file)
}
\arguments{
  \item{file}{the path to a file}
}
\value{
  a \code{data.frame}
}
\description{
  File from UCSC tend to have the column names in the first
  line preceeded with a hash character.
}
\details{
  This function imports the data and fixes the names to
  remove this hash char.
}
\author{
  Mark Cowley, 8 August 2006
}

