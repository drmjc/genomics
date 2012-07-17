#' Convert an alignment to a ucsc location.
#' 
#' alignment looks like:
#'        name   chrom   chromStart  chromEnd   strand
#' 1 NM_011547   chr13   40727709    40745894   -
#' 2 NM_010791   chr11   101693600   101710444  -
#' 3 NM_013684   chr17   15204872    15222390   +
#' 
#' @param x an alignment object (see Details)
#' @return a ucsc string, eg \dQuote{chr13:40727709-40745894}
#' @author Mark Cowley
#' @export
alignment2ucsc <- function(x) {
	#
	# find the column indices in x that contain the
	# chromosome, start, stop.
	#
	# NB: chromosome is harder to find since most of the other columns
	# also contain the word "chrom", so find the other colnames that
	# relate to start, stop, then exclude these from the search
	# for the column that contains just "chrom", "chr", or "Chrom", "Chr".
	#
	if( !is.null(colnames(x)) ) {
		start.col <- grep("start", colnames(x), ignore.case=TRUE)
		stop.col <- grep("stop|end", colnames(x), ignore.case=TRUE)
		chr.col <- setdiff(grep("chr", colnames(x), ignore.case=TRUE), c(start.col, stop.col))
	}
	else {
		chr.col <- 1
		start.col <- 2
		stop.col <- 3
	}
	stopifnot( alleq( c(1, length(start.col), length(stop.col), length(chr.col)) ) )

	res <- apply(x[,c(chr.col, start.col, stop.col)], 1, make.ucsc.string)

	return( res )
}



#' Convert an alignment to a ucsc location.
#' 
#' alignment looks like:
#'        name   chrom   chromStart  chromEnd   strand
#' 1 NM_011547   chr13   40727709    40745894   -
#' 2 NM_010791   chr11   101693600   101710444  -
#' 3 NM_013684   chr17   15204872    15222390   +
#' 
#' @param chr the chromosome names
#' @param start the start locations
#' @param stop the stop/end locations
#' @return a ucsc string, eg \dQuote{chr13:40727709-40745894}
#' @author Mark Cowley
#' @seealso \code{\link{make.ucsc.string}}
#' @export
ucsc.pos <- function(chr, start=NULL, stop=NULL) {
	make.ucsc.string(chr=chr, start=start, stop=stop)
}
# # updated 14/4/07: fixed when there are NA's in the input.


#' Convert an alignment to a ucsc location.
#' 
#' alignment looks like:
#'        name   chrom   chromStart  chromEnd   strand
#' 1 NM_011547   chr13   40727709    40745894   -
#' 2 NM_010791   chr11   101693600   101710444  -
#' 3 NM_013684   chr17   15204872    15222390   +
#' 
#' @param chr the chromosome names
#' @param start the start locations
#' @param stop the stop/end locations
#' @return a ucsc string, eg \dQuote{chr13:40727709-40745894}
#' @author Mark Cowley
#' @export
make.ucsc.string <- function(chr, start=NULL, stop=NULL) {
#	  opd <- options()$digits
#	  on.exit(options(digits=opd))
#	  options(digits=11) # avoid wrapping long integers into scientific: 1.38e09

	if( is.matrix.like(chr) ) {
		stopifnot( ncol(chr) %in% c(3,4) )
		apply(chr, 1, make.ucsc.string)
	}
	else {
		if( is.null(start) && length(chr) %in% c(3,4) ) {
			start <- min(as.numeric(c(chr[2], chr[3])))
			stop <- max(as.numeric(c(chr[2], chr[3])))

			chr <- chr[1]
		}

		if( any(is.na(chr), is.na(start), is.na(stop)) ) {
			# remove the NA's, search for remaining "good" results,
			# then return results with NA's inserted where there was an NA
			# in original data.
			na.idx <- which(is.na(chr) | is.na(start) | is.na(stop))
			ok.idx <- setdiff(1:length(chr), na.idx)
			res <- rep(NA, length(chr))
			if( length(ok.idx) > 0 )
				res[ok.idx] <- make.ucsc.string(chr[ok.idx], start[ok.idx], stop[ok.idx])
			return( res )
		}
		else {
			if( any(grepl("///", chr)) ) {
				chr <- o2n( chr )[[1]]
				start <- o2n( start )[[1]]
				stop <- o2n( stop )[[1]]
				res <- make.ucsc.string(cbind(chr=chr, start=start, stop=stop))
				return( n2o(list(res)) )
#			  res <- character( length=length(chr) )
#			  for(i in 1:length(chr))
#				  res[i] <- make.ucsc.string(chr[i], start[i], stop[i])
#			  return( n2o(list(res)) )
			}
			else {
				chr <- sub("^chr", "", chr)
				pos <- sprintf("chr%s:%d-%d", chr, start, stop)
				return( pos )
				# #return( paste(chr, ":", format(start, scientific=FALSE), "-", format(stop, scientific=FALSE), sep="") )
				# chr <- paste(sep="", "chr", chr)
				# return( paste(chr, ":", int2char(start), "-", int2char(stop), sep="") )
			}
		}
	}
}


#' is a UCSC string?
#' 
#' 
#' @param x a vector of UCSC strings
#' @return logical vector if the elements look like valid UCSC locations
#' @author Mark Cowley
#' @seealso \code{\link{make.ucsc.string}}
#' @export
is.ucsc.string <- function(x) {
	tmp <- strsplit(x, ":")
	chr <- sapply(tmp, "[", 1)
	mod <- sapply(tmp, "[", 2)
	tmp <- strsplit(mod, "-")
	start <- sapply(tmp, "[", 1)
	stop <- sapply(tmp, "[", 2)
	
	res <- is.na(x) | (grepl("^chr[0-9XYM]", chr) & grepl("^[0-9]+$", start) & grepl("^[0-9]+$", stop))
	res[is.na(x)] <- NA
	
	return( res )
}


#' Convert integers to characters
#' Really big integers like genomic locations can become scientific when
#' converted to strings eg "1e09" instead of "1000000000"
#' @param x an integer vector
#' @return a character vector
#' @author Mark Cowley, 10/7/07
#' @export
int2char <- function(x) {
	.trim <- function(x) {
		sub("^[ \t]*([^ \t]+.*[^ \t]+)[ \t]+$", "\\1", x, perl=TRUE)
	}
	.trim( format(x, scientific=FALSE) )
}


#' interval sizes
#' Calculate the size of each genomic interval listed in a ucsc-style string
#'
#' @param x a vector of ucsc-style location strings
#' @return an integer vector of interval sizes
#' @author Mark Cowley, 2011-08-02
#' @export
ucsc.size <- function(x) {
	pos <- ucsc.pos.split(x)
	return( pos[,3] - pos[,2] )
}

#' union of 2 genomic intervals
#'
#' @param x a vector of ucsc-style location strings
#' @param y a vector of ucsc-style location strings
#' @return the union of the 2 intervals
#' @author Mark Cowley, 2011-08-02
#' @export
ucsc.union <- function(x,y) {
	x.pos <- ucsc.pos.split(x)
	y.pos <- ucsc.pos.split(y)

	if(x.pos[1,1] == y.pos[1,1] ) {
		return(
			make.ucsc.string(x.pos[1,1],
							min(x.pos[1,2], y.pos[1,2]),
							max(x.pos[1,3], y.pos[1,3])
							)
				)
	}
	else {
		warning( paste(x, "and", y, "aren't on the same chromosome\n") )#'
		return( c(x, y) )
	}
}


#' union of multiple genomic intervals
#'
#' @param \dots at least 2 vectors of ucsc-style location strings
#' @return the union of the genomic intervals
#' @author Mark Cowley, 2011-08-02
#' @export
ucsc.unionN <- function(...) {
	args <- as.list(...)
	if( length(args) == 1 ) {
		if( is.character(args[[1]]) )
			return( args[[1]] )
		else if( is.list(args[[1]]) && is.character(args[[1]][[1]]) )
			args <- args[[1]]
	}
	#
	# If any markers are on different chromosomes, this function will
	# NOT handle it!
	#
	warn <- options()$warn
	on.exit(options(warn=warn))
	options(warn=2)

	res <- ucsc.union(args[[1]], args[[2]])
	if( length(args) > 2 ) {
		for(i in 3:length(args)) {
			res <- ucsc.union(res, args[[i]])
		}
	}

	return( res )

}

#' chromosome position to genome position
#' Useful for Manhattan style plots where markers on chr 2 should be to the 
#' right of those on chr1
#'
#' @param x a vector of ucsc-style genomic locations
#' @param marker.info a table of marker locations, with at least these columns: Mb, GenoMb, chr
#' @return vector of genomic locations
#' @author Mark Cowley, 2011-08-02
#' @export
ucsc.pos2genoMb.pos <- function(x, marker.info) {
	chrs <- unique(marker.info$chr)
	chr.geno.start <- rep(0, length(chrs))
	names(chr.geno.start) <- as.character( chrs )

	for(chr in chrs) {
		chr.geno.start[as.character(chr)] <- (marker.info$GenoMb - marker.info$Mb)[match(chr, marker.info$chr)]
	}

	xpos <- ucsc.pos.split(x)
	xpos[,2:3] <- xpos[,2:3] / 1e06

	for(chr in chrs) {
		idx <- which((xpos[,1] == chr) | (xpos[,1] == paste(sep="", "chr",chr)))
		if( length(idx) > 0 )
			xpos[idx,2:3] <- xpos[idx,2:3] + chr.geno.start[as.character(chr)]
	}

	xpos <- xpos[,2:3]
	colnames(xpos) <- c("Start", "Stop")
	return( xpos )
}
