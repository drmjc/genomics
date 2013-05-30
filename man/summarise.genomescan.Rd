\name{summarise.genomescan}
\alias{summarise.genomescan}
\title{summarise.genomescan}
\usage{
  summarise.genomescan(pvals,
    thresholds = c(0.05, 0.001, 1e-04, 1e-05, 1e-06, 1e-07),
    genes = NULL, file = NULL, which.genes = F,
    which.markers = F)
}
\arguments{
  \item{pvals}{a GxM matrix, where G is no. of genes, M is
  no of markers}

  \item{thresholds}{the various levels at which to
  threshold the data}

  \item{genes}{if NULL, then all genes are used;
  alternatively, it can be a list of numeric vectors which
  will be used to index the rows in the pvals table}

  \item{file}{the name of a text file for the results to be
  written to}
}
\value{
  if genes == NULL then a named vector of counts of genes
  that pass the thresholds is produced. If genes != NULL,
  then a table of counts is produced, with each set of
  genes in a different appropriately named column.
}
\description{
  Function to threshold a pvalue table at various levels,
  and report how many genes pass those linkage thresholds.
  You can look at various subsets of genes from within the
  pval table, by supplying a list argument to 'genes' where
  each element of the list is a vector of row indices.
}
\examples{
\dontrun{
	# summarise.genomescan(brain.M.spt$p, genes=genes, file="brain.M.spt.summary.txt")
}
}
\author{
  Mark Cowley, 22 Feb 2006
}

