#' Method to determine the distance between two genomic regions in "ucsc.pos"
#' format.
#' 
#' @return an integer, or Inf
#' @author Mark Cowley, 18/4/07
#' @examples
#' ucsc.pos.diff("chr9:1000-2000", "chr9:1500-2100", "min")
#' # [1] 0
#' ucsc.pos.diff("chr9:1000-2000", "chr9:1500-2100", "start")
#' # [1] 500
#' ucsc.pos.diff("chr9:1000-2000", "chr9:1500-2100", "stop")
#' # [1] 100
#' ucsc.pos.diff("chr9:1000-2000", "chr9:1500-2100", "middle")
#' # [1] 300
#' ucsc.pos.diff("chr9:1000-2000", "chr9:2100-2200", "min")
#' # [1] 100
#' ucsc.pos.diff("chr9:1000-2000", "chr9:2100-2200", "start")
#' # [1] 1100
#' ucsc.pos.diff("chr9:1000-2000", "chr9:2100-2200", "stop")
#' # [1] 200
#' ucsc.pos.diff("chr9:1000-2000", "chrX:2100-2200", "stop")
#' # [1] Inf
#' @export
ucsc.pos.diff <- function(from, to, method=c("min", "start", "stop", "middle")) {
    method <- method[1]

    from <- ucsc.pos.split(from)
    to <- ucsc.pos.split(to)

    res <- rep(Inf, nrow(from))

    na.idx <- which( rowSums(is.na(from)) > 0 | rowSums(is.na(to)) > 0 )
    if( length(na.idx) > 0 )
        res[na.idx] <- NA

    samechr.idx <- which(from[, 1] == to[, 1])
    if (length(samechr.idx) == 0)
        return(res)

    for (i in samechr.idx) {
        if (method == "start")
            res[i] <- to[i, 2] - from[i, 2]
        else if (method == "stop")
            res[i] <- to[i, 3] - from[i, 3]
        else if (method == "middle")
            res[i] <- median(c(to[i, 2], to[i, 3])) - median(c(from[i,
                2], from[i, 3]))
        else if (method == "min") {
            # to ends before from starts; thus negative distance
            if( to[i,3] < from[i,2] ) {
                res[i] <- to[i,3] - from[i,2]
            }
            # to starts after from ends; thus positive distance
            else if ( to[i,2] > from[i,3] ) {
                res[i] <- to[i,2] - from[i,3]
            }
            # at-least a partial overlap
            else {
                res[i] <- 0
            }
        }
#             .ucsc.pos.partial.overlap <- function(from, to) { # had to fix this. if to was within from then FALSE was returned.
#                 return(
#                         in.range(from[1, 2], to[1, c(2, 3)]) |
#                         in.range(from[1, 3], to[1, c(2, 3)]) |
#                         in.range(to[1, 2], from[1, c(2, 3)]) |
#                         in.range(to[1, 3], from[1, c(2, 3)])
#                         )
#             }
#             if (.ucsc.pos.partial.overlap(from[i,], to[i,])) {
#                 res[i] <- 0
#             }
#             else {
#                 tmp <- c(to[i, 2] - from[i, 3], to[i, 3] - from[i,
#                   2])
#                 res[i] <- tmp[which.min(abs(tmp))]
#             }
#         }

#         # depreciated method:
#         else if (method == "max") {
#             .ucsc.pos.complete.overlap <- function(from, to) {
#                 widths <- c(from[1, 3] - from[1, 2] + 1, to[1, 3] - to[1, 2] + 1)
#                 if (widths[1] > widths[2]) {
#                     swap <- from
#                     from <- to
#                     to <- swap
#                 }
#                 res <- ((from[1, 2] >= to[1, 2]) & (from[1, 2] <= to[1, 3])) ||
#                        ((from[1, 3] >= to[1, 2]) & (from[1, 3] <= to[1, 3]))
#                 return(res)
#             }
#             if (.ucsc.pos.complete.overlap(from[i, ], to[i, ])) {
#                 res[i] <- 0
#             }
#             else {
#                 tmp <- c(to[i, 2] - from[i, 3], to[i, 3] - from[i, 3], to[i, 2] - from[i, 2], to[i, 3] - from[i, 2])
#                 res[i] <- tmp[which.max(abs(tmp))]
#             }
#         }
        else {
            stop("You have specified an undefined method\n")
        }
    }
    return(res)
}

#' ucsc.pos.split
#' 
#' split UCSC-style locations to \code{data.frame} style
#'
#' @param x character vector of UCSC-style locations
#' @return a 3 column \code{data.frame} 'chr', 'start', 'stop'
#' @author Mark Cowley, 2013-05-30
#' @export
ucsc.pos.split <- function(x) {
    x <- gsub(",","",x)
    tmp <- strsplit(x, ":")
    chr <- sapply(tmp, "[", 1)
    tmp <- sapply(tmp, "[", 2)
    tmp <- strsplit(tmp, "-")
    start <- as.numeric(sapply(tmp, "[", 1))
    stop <- as.numeric(sapply(tmp, "[", 2))
    res <- as.data.frame(cbind(chr = chr, start = start, stop = stop))
    colclasses(res) <- c("character", "numeric", "numeric")
    return(res)
}


#' ucsc.pos.add
#' 
#' Shift a set of UCSC-style intervals a certain number of bases along the
#' same chromosome. This can be + or -, but no error checking is done, only
#' warnings made if start locations become negative.
#' 
#' @author Mark Cowley, 2009-01-08
#' @export
ucsc.pos.add <- function(x, bases) {
	x.split <- ucsc.pos.split(x)
	x.split$start <- x.split$start + bases
	x.split$stop <- x.split$stop + bases
	if( any(x.split$start) < 0 )
		warning("Some start locations are now < 0.")

	res <- ucsc.pos(x.split)

	res
}

#' ucsc.pos.subset
#' 
#' For a set of UCSC style intervals, produce a new set of intervals that are
#' subsets of those intervals. Strand must be provided, because intervals are
#' strand-less, but the logic to subset an interval depends very much on the
#' strand.
#' 
#' eg: subset chr1:1000-1100 (+) to the interval that starts at 2-7 ->
#' chr1:1001-1006
#' eg: subset chr1:1000-1100 (-) to the interval that starts at 2-7 ->
#' chr1:1094-1099
#' 
#' @param intervals vector of UCSC style intevals, eg "chr1:100-300"
#' @param start vector of numerics defining the starting base within each
#' interval (1-based)
#' @param stop vector of numerics defining the ending base within each
#' interval (1-based). NB stop can be longer than the original interval.
#' @param strand either a charcter vector of "+"/"-", or {0,1}, or {-1,1}
#' @return a vector of UCSC style intervals, adjusted as detailed above
#' @author Mark Cowley, 2009-01-08
#' @examples
#' ucsc.pos.subset("chr1:1000-1100", 1,11, "+")
#' # [1] "chr1:1000-1010"
#' ucsc.pos.subset("chr1:1000-1100", 1,11, "-")
#' # [1] "chr1:1090-1100"
#' # (you can use > 1 entry too.)
#' @export
ucsc.pos.subset <- function(intervals, start, stop, strand) {
	if( is.numeric(strand) ) {
		if( all(strand) %in% c(-1,1) )
			strand <- ifelse(strand==-1, "-", "+")
		else if( all(strand) %in% c(0,1) )
			strand <- ifelse(strand==0, "-", "+")
		else
			stop("Can't work out what type of numeric Strand you've provided. Please try '+' or '-'.\n")
	}
	
	stopifnot( all(strand %in% c("+", "-")) )
	
	x.split <- ucsc.pos.split(intervals)
	res <- x.split
	
	# consider the sense and anti-sense strand independently
	SS <- which(strand == "+") # SS for "same strand"
	if( length(SS) > 0 ) {
		res$start[SS] <- x.split$start[SS] + start[SS] - 1
		res$stop[SS] <- x.split$start[SS] + stop[SS] - 1
	}
	SS <- which(strand == "-")
	if( length(SS) > 0 ) {
		res$start[SS] <- x.split$stop[SS] - stop[SS] + 1
		res$stop[SS] <- x.split$stop[SS] - start[SS] + 1
	}
	
	res <- ucsc.pos(res)
	res
}
