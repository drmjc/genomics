#' Change the GENOMEacc from gene2refseq into a chromosome name
#' 
#' Currently supports mouse, human and rat. For new organisms, find
#' out the chromosome accession numbers by searching for them in the
#' tmp genome page:
#' eg: try searching for "rattus AND norvegicus AND chromosome" in
#' http://www.tmp.nlm.nih.gov/entrez/query.fcgi?CMD=search&DB=genome
#' 
#' @author Mark Cowley, 14 April 2006
#' @export
accession2chr <- function(x, taxid=10090) {
    if( taxid == TAXID.MOUSE ) {
        chr <- c( "NC_000067", # chr 1
                  "NC_000068", # chr 2
                  "NC_000069", # chr 3
                  "NC_000070", # chr 4
                  "NC_000071", # chr 5
                  "NC_000072", # chr 6
                  "NC_000073", # chr 7
                  "NC_000074", # chr 8
                  "NC_000075", # chr 9
                  "NC_000076", # chr 10
                  "NC_000077", # chr 11
                  "NC_000078", # chr 12
                  "NC_000079", # chr 13
                  "NC_000080", # chr 14
                  "NC_000081", # chr 15
                  "NC_000082", # chr 16
                  "NC_000083", # chr 17
                  "NC_000084", # chr 18
                  "NC_000085", # chr 19
                  "NC_000086", # chr X
                  "NC_000087", # chr Y
                  "NC_005089") # Mitochondrion

        names(chr) <- c(as.character(1:19), "X", "Y", "M")
    }
    else if( taxid == TAXID.HUMAN ) {
    # http://www.ncbi.nlm.nih.gov/genome/guide/human/release_notes.html
        chr <- c( "NC_000001", # chr 1
                  "NC_000002", # chr 2
                  "NC_000003", # chr 3
                  "NC_000004", # chr 4
                  "NC_000005", # chr 5
                  "NC_000006", # chr 6
                  "NC_000007", # chr 7
                  "NC_000008", # chr 8
                  "NC_000009", # chr 9
                  "NC_000010", # chr 10
                  "NC_000011", # chr 11
                  "NC_000012", # chr 12
                  "NC_000013", # chr 13
                  "NC_000014", # chr 14
                  "NC_000015", # chr 15
                  "NC_000016", # chr 16
                  "NC_000017", # chr 17
                  "NC_000018", # chr 18
                  "NC_000019", # chr 19
                  "NC_000020", # chr 20
                  "NC_000021", # chr 21
                  "NC_000022", # chr 22
                  "NC_000023", # chr X
                  "NC_000024", # chr Y
                  "NC_001807") # Mitochondrion

        names(chr) <- c(as.character(1:22), "X", "Y", "M")
    }
    else if( taxid == TAXID.RAT ) {
        chr <- c("NC_005100", # chr 1
                 "NC_005101", # chr 2
                 "NC_005102", # chr 3
                 "NC_005103", # chr 4
                 "NC_005104", # chr 5
                 "NC_005105", # chr 6
                 "NC_005106", # chr 7
                 "NC_005107", # chr 8
                 "NC_005108", # chr 9
                 "NC_005109", # chr 10
                 "NC_005110", # chr 11
                 "NC_005111", # chr 12
                 "NC_005112", # chr 13
                 "NC_005113", # chr 14
                 "NC_005114", # chr 15
                 "NC_005115", # chr 16
                 "NC_005116", # chr 17
                 "NC_005117", # chr 18
                 "NC_005118", # chr 19
                 "NC_005119", # chr 20
                 "NC_005120") # chr X

        names(chr) <- c(as.character(1:20), "X")
    }
    else {
        stop( "unsupported taxon. see function description for adding new taxa." )
    }

    res <- rep(NA, length(x))
    for(i in 1:length(chr)) {
##         res[grep(paste0(chr[i], "\\."), x)] <- names(chr)[i]
        res[grep(chr[i], x)] <- names(chr)[i]
    }

    return( res )
}
