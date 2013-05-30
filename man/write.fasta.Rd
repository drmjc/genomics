\name{write.fasta}
\alias{write.fa}
\alias{write.fasta}
\title{write a FastA file}
\usage{
  write.fasta(sequences, names = names(sequences),
    file = "", subset = NULL)

  write.fa(sequences, names = names(sequences), file = "",
    subset = NULL)
}
\arguments{
  \item{sequences}{a character vector of sequences}

  \item{names}{the names of each fasta record}

  \item{file}{the output file name}

  \item{subset}{optional character vector of fasta records
  to export}
}
\value{
  nothing.
}
\description{
  write a FastA file
}
\author{
  Mark Cowley, 2012-07-06
}

