\name{import.miRbase.miR2mir}
\alias{import.miRbase.miR2mir}
\title{import.miRbase.miR2mir}
\usage{
  import.miRbase.miR2mir(f = "~/data/miRBase/sequences/latest/miRNA.csv",
    species = "hsa", as.list = FALSE)
}
\arguments{
  \item{f}{the lastest miRBase miRNA.csv file (from Sanger
  FTP)}

  \item{species}{which 3 letter species code? (default
  hsa)}

  \item{as.list}{if FALSE, the results are a table, it
  TRUE, the result is a list with n elements, 1 per miR}
}
\description{
  Import the master miRBase miRNA table, but re-order the
  relationships to be from miR to mir
}
\author{
  Mark Cowley, 2008-09-22
}

