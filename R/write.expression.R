# wrapper to write expression data to a tab delimited, non-quoted
# table with the rownames in the first column with a column header
# which makes it suitable for opening in a spreadsheet (which would
# normally puts the column names starting from 1st col instead of
# the 2nd col)
#
# Mark Cowley, 25 May 2006
#
write.expression.data <- function(x, file, idtype="GenBankID", digits=4, na="") {
	cc <- which(colclasses(x) == "numeric")
	x[,cc] <- round(x[,cc], digits)
    write.delim(x, file, row.names=idtype, na=na)
}

# alias for write.expression.data
#
# Mark Cowley, 25 May 2006
write.expression <- function(x, file, idtype="GenBankID", digits=4, na="") {
    write.expression.data(x, file, idtype, digits, na=na)
}
