\name{miRbase.miR2mir}
\alias{miRbase.miR2mir}
\title{Convert mirbase database
Convert the miRbase db (which has one row per precursor miRNA) to have one
row per mature miRNA.}
\usage{
  miRbase.miR2mir(mir2miR, collapse = TRUE, sep = " // ")
}
\arguments{
  \item{mir2miR}{see \code{\link{import.miRbase}}}

  \item{collapse}{logical: collapse multiple rows into one?
  separated by \code{sep}}

  \item{sep}{see \code{collapse}}
}
\value{
  a \code{data.frame} with columns: \dQuote{Mature_Acc},
  \dQuote{Mature_ID}, \dQuote{Mature_Seq},
  \dQuote{Precursor_Acc}, \dQuote{Precursor_ID},
  \dQuote{Precursor_Status}, \dQuote{Precursor_Seq},
  \dQuote{miRBase.order}. This final column indicates if it
  was in the Mature1 or Mature2 columns from the original
  miRNA file
}
\description{
  Convert mirbase database Convert the miRbase db (which
  has one row per precursor miRNA) to have one row per
  mature miRNA.
}
\author{
  Mark Cowley, 2009-01-07
}

