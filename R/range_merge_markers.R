#' merge a range of markers
#' 
#' Usful function to annotate the markers that are significant for a
#' particular gene in a genome scan. First work out which markers are
#' significant (either as indices, or as the marker names) then supply these
#' markers to this function, along with the names of all markers, and you will
#' be supplied with a char vector of these markers, and if a few markers in a
#' row are significant, then just the start and end marker with a "-" (or sep)
#' is returned.
#' 
#' @param x 1:length(marker.names)
#' @param marker.names The marker names (or words) to look up the indices
#' supplied in x.
#' @param gaplength allow N gaps between consecutive markers to be allowed. 0
#' implies no gaps; 1 implies if markers A and C are significant, then marker
#' B is expected to be significant, and included in the range.
#' @param sep Which character should be used to delimit a range of words?
#' @param single.string if NULL, then a character vector is returned. If not
#' NULL then it's value will be used to concatenate the ranges into a single
#' string, thus ", " comma separates the ranges. For clarity, avoid using the
#' same symbol as specified by sep.
#' @return A character vector with ranges identified by the first and last
#' marker, and singletons identified by their name.  
#' @author Mark Cowley, 22 march 2006
#' @export
#' @importFrom mjcbase range_merge
#' @examples
#' range_merge_markers(c(1:4,6), LETTERS)
#' # [1] "A-D" "F"
#' range_merge_markers(c(LETTERS[1:4], LETTERS[6]), LETTERS)
#' # [1] "A-D" "F"
range_merge_markers <- function(x, marker.names, gaplength=0, sep="-", single.string=", ") {

    if( missing(marker.names) )
        stop("must supply marker.names to work out if there are any gaps in 'x'\n")

    if( is.character(x) ) {
        #
        # convert the characters to indices into the marker.names vector
        # so the programming logic can work
        #

        x <- match(x, marker.names)
    }

    range <- range_merge(x, gaplength=gaplength)

    if( !is.null(range) && nrow(range) > 0 ) {
        res <- rep("", nrow(range))
        for(i in 1:nrow(range)) {
            if( range[i,1] == range[i,2] )
                res[i] <- marker.names[ range[i,1] ]
            else
                res[i] <- paste( marker.names[ range[i,1] ], marker.names[ range[i,2] ], sep=sep)
        }
    }

    if( !is.null(single.string) )
        res <- paste(res, collapse=single.string, sep="")

    return( res )
}
