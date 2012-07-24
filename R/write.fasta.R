#' write a FastA file
#'
#' @param sequences a character vector of sequences
#' @param names the names of each fasta record
#' @param file the output file name
#' @param subset optional character vector of fasta records to export
#' 
#' @return nothing.
#' 
#' @author Mark Cowley, 2012-07-06
#' @export
write.fasta <- function(sequences, names=names(sequences), file="", subset=NULL) {
	if( is.null(subset) )
		subset <- c(1:length(sequences))

	if( is.character(subset) ) {
		subset <- match(intersect(subset, names), names)
	}
	else {
		# subset <- rm.na(subset)
		subset <- na.omit(subset)
	}
	
	tmp <- p(">", names[subset], "\n", sequences[subset])
	write.delim(tmp, file, col.names=F, row.names=F)
}


#' @export
#' @aliases write.fa
#' @rdname write.fasta
write.fa <- write.fasta
