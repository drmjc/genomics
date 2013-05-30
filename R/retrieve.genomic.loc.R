#' retrieve.genomic.loc
#' 
#' Lookup a GeneID at NCBI and see if there is a refseq record that has been
#' mapped to the genome. Why? for some reason the gene2refseq file does not
#' reflect the NCBI website all the time? strange.
#' 
#' @param geneID an Entrez GeneID
#' @param override.cache a ncbi copy of the webpage is saved. This can be
#' overwritten to get the very latest webpage if you desire.
#' @param sleep if length(x) > 1, how many seconds should we sleep between
#' batch commands?
#' @author Mark Cowley, 10 April 2006
#' @export
#' @importFrom mjcbase scan.text
#' @examples
#' \dontrun{
#' retrieve.genomic.loc( 11643 )
#' }
#' 
retrieve.genomic.loc <- function(geneID, override.cache=FALSE, cache.dir="~/data/mirror/gene.searches", sleep=1) {

    if( !file.exists(cache.dir) ) {
        stop("the cache.dir does not exist. recommend setting cache.dir=\"tempdir()\".\n")
    }

    COLNAMES <- c("GeneID", "GENOMEacc", "Start", "Stop", "Strand", "Chr")

    #
    # if > 1 gene has been specified, loop through them all, and download info if necessary.
    #
    if( length(geneID) > 1 ) {
        res <- matrix(NA, length(geneID), length(COLNAMES))
        colnames(res) <- COLNAMES
        for(i in 1:length(geneID)) {
            incache <- file.exists( file.path(cache.dir, paste0(as.character(geneID[i]), ".htm")) )

            tmp <- retrieve.genomic.loc( geneID[i], override.cache=override.cache,
                                            cache.dir=cache.dir )

            if( length(tmp) == 0 ) {
                res[i,1] <- geneID[i]
            }
            else {
                res[i, ] <- tmp
            }

            if( !incache && i < length(geneID) ) {
                cat(".")
                Sys.sleep(sleep)
            }
        }

        #
        # clean up res's formatting.
        #
        res <- as.data.frame(res)
        for(i in 1:ncol(res)) {
            res[,i] <- as.character(res[,i])
            if(i %in% c(1,3,4))
                res[,i] <- as.numeric(res[,i])
        }

        return( res )
    }


    #
    # NCBI Entrez Gene weblink in the "Gene Table" view which is the least amount of data
    # that will indicate whether the refseq record has been mapped or not.
    #
    pre <- "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=gene_table&list_uids="

    #
    # primite use of cache files
    #
    tmp.result.file <- file.path(cache.dir, paste0(as.character(geneID), ".htm"))

    if( override.cache || !file.exists(tmp.result.file ) ) {
        #
        # download and cache the weblink from NCBI
        #
        ncbi <- scan.text(paste0(pre, geneID))
        write(ncbi, tmp.result.file )
    }
    else {
##         cat("using cache\n")
##         cat("|")
        ncbi <- scan.text(tmp.result.file )
    }

    idx <- grep("NC_", ncbi)
    #
    # if "NC_" is found in the webpage then in all cases checked so far,
    # there is a RefSeq record that has been aligned to a chromosome.
    #
    if( length(idx) >= 1 ) {

        #
        # Choose the first time that this "NC_" pattern is spotted...
        # and parse out the chromosome, start, stop and strand information.
        #
        tmp <- ncbi[ idx[1] ]
        ## tmp="["FASTA","window.location='/entrez/viewer.fcgi?val=NC_000086.4&from=5501393&to=5511569&dopt=fasta'","",""],"

        tmp <- gsub("^.+\\?", "", tmp, perl=T)
        ## tmp="val=NC_000086.4&from=5501393&to=5511569&dopt=fasta'","",""],"

        tmp <- gsub("&dopt.*", "", tmp)
        ## tmp="val=NC_000086.4&from=5501393&to=5511569"

        parts <- strsplit(tmp, split="&")[[1]]
        ## parts=c("val=NC_000086.4", "from=5501393", "to=5511569")

        for(i in 1:3)
            parts[i] <- gsub(".*=", "", parts[i])

        GENOMEacc <- parts[1]
        start <- as.numeric( parts[2] )
        stop <-  as.numeric( parts[3] )
        strand <- ifelse( grepl("minus strand", ncbi), "-", "+" )
        chr <- accession2chr( parts[1] )

        res <- c(geneID, GENOMEacc, start, stop, strand, chr)
        names(res) <- COLNAMES

        return(res)
    }
    else {
        return( character(0) )
    }

}
