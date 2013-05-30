#' rename.gnf.markers
#' 
#' Simple function to rename Gnf markers to their proper name: eg: rename
#' 02.118.650 -> S02Gnf118.650 It's kinda slow, so avoid calling it in a loop
#' 
#' @author Mark Cowley, 22 March 2006
#' @export
rename.gnf.markers <- function(x) {
    for(i in 1:length(x)) {
        if( grepl("\\..*\\.", x[i]) ) {
            # then the marker looks like: '02.118.650'
            # needs to look like: 'S02Gnf118.650'
            tmp <- strsplit(x[i], "\\.")[[1]]
            x[i] <- paste("S", tmp[1], "Gnf", tmp[2], ".", tmp[3], collapse="", sep="")
        }
    }

    return(x)
}
