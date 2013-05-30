\name{import.miRbase.locations}
\alias{import.miRbase.locations}
\title{Import the genomic locations of the miRBase precursor miRNA's.}
\usage{
  import.miRbase.locations(gff.file)
}
\arguments{
  \item{gff.file}{the path to the gff file, avalable from
  references [1-3]}
}
\value{
  see \code{\link{import.gff}}, with 3 additional columns:
  \dQuote{ACC}, \dQuote{ID} and \dQuote{pos}, eg
  \dQuote{MI0006363}, \dQuote{hsa-mir-1302-2},
  \dQuote{chr1:20229-20366}
}
\description{
  Import the genomic locations of the miRBase precursor
  miRNA's.
}
\examples{
\dontrun{
tmp <- import.miRbase.locations("~/data/miRBase/sequences/12.0/genomes/hsa.gff")
head(tmp)
#         ACC           ID chr     start       end strand score                       pos
# 1 MI0000060 hsa-let-7a-1   9  95978060  95978139      +     .    chr9:95978060-95978139
# 2 MI0000061 hsa-let-7a-2  11 121522440 121522511      -     . chr11:121522440-121522511
# 3 MI0000062 hsa-let-7a-3  22  44887293  44887366      +     .   chr22:44887293-44887366
# 4 MI0000063   hsa-let-7b  22  44888230  44888312      +     .   chr22:44888230-44888312
# 5 MI0000064   hsa-let-7c  21  16834019  16834102      +     .   chr21:16834019-16834102
# 6 MI0000065   hsa-let-7d   9  95980937  95981023      +     .    chr9:95980937-95981023
}
}
\author{
  Mark Cowley, 2009-01-07
}
\references{
  [1]
  \url{ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/genomes/hsa.gff}
  [2]
  \url{ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/genomes/mmu.gff}
  [3]
  \url{ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/genomes/rno.gff}
}

