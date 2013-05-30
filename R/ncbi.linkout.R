#' Generate a weblink for a given Entrez GeneID to the Entrez gene NCBI
#' database.
#' 
#' @param x vector of Entrez Gene ID's
#' @param g2r ignored.
#' 
#' @return character vector of URL's
#' 
#' @author Mark Cowley, 10 April 2006
#' @export
linkout.gene <- function(x, g2r) {
    pre <- "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=full_report&list_uids="

    return( paste(pre, x, sep="") )
}

#' generate a weblink for a given Entrez GeneID to the UCSC genome browser.
#' 
#' @param x vector of Entrez Gene ID's
#' @param g2r undocumented
#' 
#' @return character vector of URL's
#' 
#' @author Mark Cowley, 10 April 2006
#' @export
linkout.ucsc <- function(x, g2r) {

    row <- grep(x, g2r$GeneID)

    if( "-" %in% g2r[row, c("GENOMEacc", "Start", "Stop")] )
        return( NULL )

    chr <- accession2chr( g2r$GENOMEacc[row] )
    start <- g2r$Start[row]
    stop <- g2r$Stop[row]
    pre <- "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=vertebrate&org=Mouse&db=mm6&position="

    return( paste(pre, "chr", chr, ":", start, "-", stop, sep="") )
}

