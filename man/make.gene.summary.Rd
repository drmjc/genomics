\name{make.gene.summary}
\alias{make.gene.summary}
\title{Make an Entrez Gene Summary table which includes the genes symbol,
description, genomic mapping info and synonyms}
\usage{
  make.gene.summary(outfile = "/biomirror/mouse/Gene/gene.summary.mm.txt",
    taxid = 10090)
}
\arguments{
  \item{outfile}{undocumented parameter}

  \item{taxid}{undocumented parameter}
}
\value{
  a \code{data.frame} with the following columns:
  c("GeneID", "Symbol", "Description", "Synonym",
  "GenbankID", "mgiID", "Chr", "Start", "Stop", "Strand",
  "cM") GeneID Symbol Description 1 11287 Pzp pregnancy
  zone protein 2 11298 Aanat arylalkylamine
  N-acetyltransferase 3 11302 Aatk apoptosis-associated
  tyrosine kinase 4 11303 Abca1 ATP-binding cassette,
  sub-family A (ABC1), member 1 5 11304 Abca4 ATP-binding
  cassette, sub-family A (ABC1), member 4 6 11305 Abca2
  ATP-binding cassette, sub-family A (ABC1), member 2
  Synonyms GenbankID mgiID 1 A1m /// A2m /// AI893533 ///
  MAM NM_007376 MGI:87854 2 MGC151344 /// Nat-2 /// Nat4
  /// SNAT NM_009591 MGI:1328365 3 AATYK /// mKIAA0641
  NM_007377 MGI:1197518 4 Abc1 NM_013454 MGI:99607 5
  AW050280 /// Abc10 /// Abcr /// D430003I15Rik /// RmP
  NM_007378 MGI:109424 6 AI413825 /// Abc2 /// D2H0S1474E
  /// mKIAA1062 NM_007379 MGI:99606 Chr Start Stop Strand
  cM 1 6 1.28e+08 1.28e+08 - 62.0 2 11 1.16e+08 1.16e+08 +
  70.0 3 11 1.20e+08 1.20e+08 - NA 4 4 5.31e+07 5.32e+07 -
  23.1 5 3 1.22e+08 1.22e+08 + 61.8 6 2 2.53e+07 2.53e+07 +
  12.6 checked with Mouse, Rat and Humans.
}
\description{
  Data is retreived from "gene2refseq" and "gene_info"
  available from:
  "ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz" and
  "ftp.ncbi.nih.gov/gene/DATA/gene_info.gz" Algorithm: The
  best genomic location is determined from gene2refseq file
  (see import.gene2refseq) The relevant columns from
  gene_info are kept (see import.gene_info)
}
\author{
  Mark Cowley, 19/3/07
}

