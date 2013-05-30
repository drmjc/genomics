\name{import.gff}
\alias{import.gff}
\title{Import a GFF file}
\usage{
  import.gff(file, row.pattern = NULL)
}
\arguments{
  \item{file}{the path to a gff file}

  \item{row.pattern}{currently unused}
}
\value{
  a \code{data.frame} with these columns: "chr", "source",
  "feature", "start", "end", "score", "strand", "frame",
  "group"
}
\description{
  from:
  \url{http://genome.ucsc.edu/FAQ/FAQformat#format3}\cr GFF
  (General Feature Format) lines are based on the GFF
  standard file format. GFF lines have nine required fields
  that must be tab-separated. If the fields are separated
  by spaces instead of tabs, the track will not display
  correctly. For more information on GFF format, refer to
  \url{http://www.sanger.ac.uk/Software/formats/GFF}. Here
  is a brief description of the GFF fields:\cr seqname -
  The name of the sequence. Must be a chromosome or
  scaffold.\cr source - The program that generated this
  feature.\cr feature - The name of this type of feature.
  Some examples of standard\cr feature types are "CDS",
  "start_codon", "stop_codon", and "exon".\cr start - The
  starting position of the feature in the sequence. The
  first base is numbered 1.\cr end - The ending position of
  the feature (inclusive).\cr score - A score between 0 and
  1000. If the track line useScore attribute is set to 1
  for this annotation data set, the score value will
  determine the level of gray in which this feature is
  displayed (higher numbers = darker gray). If there is no
  score value, enter ".".\cr strand - Valid entries include
  '+', '-', or '.' (for don't know/don't care).\cr frame -
  If the feature is a coding exon, frame should be a number
  between 0-2 that represents the reading frame of the
  first base. If the feature is not a coding exon, the
  value should be '.'.\cr group - All lines with the same
  group are linked together into a single item.\cr
}
\author{
  Mark Cowley, 2/1/08
}

