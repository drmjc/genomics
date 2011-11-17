
write.fa <- function(sequences, names=names(sequences), file="", subset=NULL) {
    write.fasta(sequences=sequences, names=names, file=file, subset=subset)
}

write.fasta <- function(sequences, names=names(sequences), file="", subset=NULL) {
	if( is.null(subset) )
		subset <- c(1:length(sequences))

	if( is.character(subset) ) {
		subset <- match(intersect(subset, names), names)
	}
	else {
		subset <- rm.na(subset)
	}
	
	tmp <- p(">", names[subset], "\n", sequences[subset])
	write.delim(tmp, file, col.names=F, row.names=F)
}

