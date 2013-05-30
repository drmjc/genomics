#' find a unique genomic mapping for each entrez gene (if known)
#' 
#' The gene2refseq file maps many refseq records to one Entrez gene record.
#' If the refseq has been mapped to the genome, then the GENOMEacc column
#' contains a chromosomal accession no starting with NC. Note, NT records
#' are contig mappings which do not tell us genomic mapping locations.
#' Essentially, only lines that are mapped to NC are considered, and if there
#' are >1 mappings, then the record with the stronger refseq status is
#' considered (see rank.refseq). If there's still a tie then the row closest
#' to the top of the file is used.
#' NB: for mouse build 35.1, 30172 out of 43981 GeneID records were mapped
#' successfully (ie 69%) using gene2refseq downloaded on 2/3/06.
#' 
#' @param g2acc undocumented
#' @return undocumented
#' 
#' @author Mark Cowley, 10 April 2006
#' @export
make.unique.gene2accession <- function( g2acc ) {
    make.unique.gene2refseq( g2acc )
}

#' find a unique genomic mapping for each entrez gene (if known)
#' 
#' The gene2refseq file maps many refseq records to one Entrez gene record.
#' If the refseq has been mapped to the genome, then the GENOMEacc column
#' contains a chromosomal accession no starting with NC. Note, NT records
#' are contig mappings which do not tell us genomic mapping locations.
#' Essentially, only lines that are mapped to NC are considered, and if there
#' are >1 mappings, then the record with the stronger refseq status is
#' considered (see rank.refseq). If there's still a tie then the row closest
#' to the top of the file is used.
#' NB: for mouse build 35.1, 30172 out of 43981 GeneID records were mapped
#' successfully (ie 69%) using gene2refseq downloaded on 2/3/06.
#' 
#' @param g2r undocumented
#' @param taxid undocumented
#' @return undocumented
#' 
#' @author Mark Cowley, 10 April 2006
#' @export
make.unique.gene2refseq <- function(g2r, taxid=TAXID.MOUSE) {
    g2r <- g2r[grep("^NC", g2r$GENOMEacc),]
    g2r <- g2r[order(g2r$GeneID),]

    #
    # since g2r is sorted via GeneID's, it'll be quicker to remember
    # which rows relate to which geneID in one go at the start that to
    # iterate through all rows searching for each geneID at a time.
    #
    genes <- unique(g2r$GeneID)
    indices <- match(genes, g2r$GeneID)
    names(indices) <- genes

    rm.idx <- NULL
    for(i in 1:length(indices)) {
        if(i == length(indices))
            rows <- c(indices[i]: nrow(g2r))
        else
            rows <- c(indices[i]:(indices[i+1]-1))

        stopifnot(alleq(g2r$GeneID[rows])) ## should NEVER fail if g2r is sorted by GeneID.

        if( length(rows) > 1 ) {
            if( !alleq(g2r$GENOMEacc[rows]) ) {
                cat( paste("GeneID: ", g2r$GeneID[rows[1]], " maps to >1 chr -- choosing the most validated one (or the first one)\n", sep="") )
            }

            rows <- rows[order.refseq(g2r$Status[rows])]
            rm.idx <- c(rm.idx, rows[2:length(rows)])
        }
    }

    if( length(rm.idx) > 0 )
        g2r <- g2r[-rm.idx, ]

    #
    # add the chromosomal info.
    #
    g2r$Chr <- accession2chr( g2r$GENOMEacc, taxid=taxid )

    return( g2r )
}



#' rank RefSeq identifiers
#' 
#' From the tmp Handbook:
#' http://www.tmp.nlm.nih.gov/books/bv.fcgi?rid=handbook.table.697
#' 
#' I've placed INFERRED as stronger evidence than PREDICTED, GA or WGS since
#' there is evidence
#' in other organisms that is exists.
#' VALIDATED >> REVIEWED >> PROVISIONAL >> INFERRED >> PREDICTED = GENOME
#' ANNOTATION = WGS
#' 
#' @param x undocumented
#' @return undocumented
#' 
#' @author Mark Cowley
#' @export
rank.refseq <- function(x) {
    ranks <- rep(5, length(x))

    ranks[grep("VALIDATED", x)] <- 1
    ranks[grep("REVIEWED", x)] <- 2
    ranks[grep("PROVISIONAL", x)] <- 3
    ranks[grep("INFERRED", x)] <- 4

    return( ranks )
}


#' order RefSeq ID's by validation status
#' 
#' return the order of a vector of RefSeq status's such that those that are
#' VALIDATED come before those that are PREDICTED.
#' 
#' Here are the status rankings:
#' VALIDATED >> REVIEWED >> PROVISIONAL >> INFERRED >> PREDICTED = GENOME
#' ANNOTATION = WGS
#' 1 2 3 4 5 5 5
#' see ?order for equivalent method using normal integers.
#' 
#' @param x undocumented
#' @param na.last see order
#' @param decreasing see order
#' 
#' @return Undocumented return value
#' 
#' @author Mark Cowley, 10 April 2006
#' @export
order.refseq <- function(x, na.last = TRUE, decreasing = FALSE) {
    return( order(rank.refseq(x), na.last=na.last, decreasing=decreasing) )
}

