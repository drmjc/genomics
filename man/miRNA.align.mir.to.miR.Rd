\name{miRNA.align.mir.to.miR}
\alias{miRNA.align.mir.to.miR}
\title{Align miR to mir sequences
Sequence based alignment of mature miR sequences within precursor mir
sequences}
\usage{
  miRNA.align.mir.to.miR(mature.id, mature.seq,
    precursor.id, precursor.seq)
}
\arguments{
  \item{mature.id}{a vector of identifiers for the mature
  and precursor}

  \item{precursor.id}{a vector of identifiers for the
  mature and precursor miR/mir genes respectively.}

  \item{mature.seq}{a vector of sequences. same length as
  each other, and}

  \item{precursor.seq}{a vector of sequences. same length
  as each other, and same length as the names}
}
\value{
  a table of N rows, containing: Mature ID, length of
  sequence, Precursor ID, length of sequence, Alignment
  score (see ?pairwiseAlignment), start and end coordinates
  from within the precursor align.prop is the align start /
  precursor length
}
\description{
  Align miR to mir sequences Sequence based alignment of
  mature miR sequences within precursor mir sequences
}
\author{
  Mark Cowley, 2009-01-13
}

