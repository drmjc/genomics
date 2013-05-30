\name{accession2chr}
\alias{accession2chr}
\title{Change the GENOMEacc from gene2refseq into a chromosome name}
\usage{
  accession2chr(x, taxid = 10090)
}
\description{
  Currently supports mouse, human and rat. For new
  organisms, find out the chromosome accession numbers by
  searching for them in the tmp genome page: eg: try
  searching for "rattus AND norvegicus AND chromosome" in
  http://www.tmp.nlm.nih.gov/entrez/query.fcgi?CMD=search&DB=genome
}
\author{
  Mark Cowley, 14 April 2006
}

