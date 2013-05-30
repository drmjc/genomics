\name{make.unique.gene2accession}
\alias{make.unique.gene2accession}
\title{find a unique genomic mapping for each entrez gene (if known)}
\usage{
  make.unique.gene2accession(g2acc)
}
\arguments{
  \item{g2acc}{undocumented}
}
\value{
  undocumented
}
\description{
  The gene2refseq file maps many refseq records to one
  Entrez gene record. If the refseq has been mapped to the
  genome, then the GENOMEacc column contains a chromosomal
  accession no starting with NC. Note, NT records are
  contig mappings which do not tell us genomic mapping
  locations. Essentially, only lines that are mapped to NC
  are considered, and if there are >1 mappings, then the
  record with the stronger refseq status is considered (see
  rank.refseq). If there's still a tie then the row closest
  to the top of the file is used. NB: for mouse build 35.1,
  30172 out of 43981 GeneID records were mapped
  successfully (ie 69\%) using gene2refseq downloaded on
  2/3/06.
}
\author{
  Mark Cowley, 10 April 2006
}

