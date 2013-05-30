\name{int2char}
\alias{int2char}
\title{Convert integers to characters}
\usage{
  int2char(x)
}
\arguments{
  \item{x}{an integer vector}
}
\value{
  a character vector
}
\description{
  Really big integers like genomic locations can become
  scientific when converted to strings eg "1e09" instead of
  "1000000000"
}
\author{
  Mark Cowley, 10/7/07
}

