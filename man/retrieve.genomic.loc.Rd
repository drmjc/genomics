\name{retrieve.genomic.loc}
\alias{retrieve.genomic.loc}
\title{retrieve.genomic.loc}
\usage{
  retrieve.genomic.loc(geneID, override.cache = FALSE,
    cache.dir = "~/data/mirror/gene.searches", sleep = 1)
}
\arguments{
  \item{geneID}{an Entrez GeneID}

  \item{override.cache}{a ncbi copy of the webpage is
  saved. This can be overwritten to get the very latest
  webpage if you desire.}

  \item{sleep}{if length(x) > 1, how many seconds should we
  sleep between batch commands?}
}
\description{
  Lookup a GeneID at NCBI and see if there is a refseq
  record that has been mapped to the genome. Why? for some
  reason the gene2refseq file does not reflect the NCBI
  website all the time? strange.
}
\examples{
\dontrun{
retrieve.genomic.loc( 11643 )
}
}
\author{
  Mark Cowley, 10 April 2006
}

