\name{miRNA.which.arm}
\alias{miRNA.which.arm}
\title{Predict miRNA arm
Predict which arm (3' or 5') of the precursor miRNA that the mature miR is
on This uses purely the sequences of the precursor, and the mature miRNA.}
\usage{
  miRNA.which.arm(mature.seq, precursor.seq, THRESH = 0.44,
    detailed.alignment = NULL)
}
\arguments{
  \item{mature.seq}{a vector of N mature miRNA sequences}

  \item{precursor.seq}{a vector of N precursor miRNA
  sequences}

  \item{THRESH}{The alignment score (see notes) that
  distinguishes the 3'/5' arms}

  \item{detailed.alignment}{If you've already run
  \code{\link{miRNA.align.mir.to.miR}}, then you can
  provide this argument, instead of \code{mature.seq} and
  \code{precursor.seq}}
}
\value{
  a vector of N values that can be either "5p" or "3p"
}
\description{
  Predict miRNA arm Predict which arm (3' or 5') of the
  precursor miRNA that the mature miR is on This uses
  purely the sequences of the precursor, and the mature
  miRNA.
}
\note{
  Alignment score is the alignments' start index relative
  to the precursor divided by the length of the precursor
}
\author{
  Mark Cowley, 2009-01-08
}

