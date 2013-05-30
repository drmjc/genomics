\name{write.expression.data}
\alias{write.expression}
\alias{write.expression.data}
\title{write.expression.data}
\usage{
  write.expression.data(x, file, idtype = "GenBankID",
    digits = 4, na = "")

  write.expression(x, file, idtype = "GenBankID",
    digits = 4, na = "")
}
\description{
  wrapper to write expression data to a tab delimited,
  non-quoted table with the rownames in the first column
  with a column header which makes it suitable for opening
  in a spreadsheet (which would normally puts the column
  names starting from 1st col instead of the 2nd col)
}
\author{
  Mark Cowley, 25 May 2006
}

