
## ## update.mapping <- function(genbankIDs, g2acc=NULL, ginfo=NULL, taxid=10090, mirror="/import/sparta/data/mark/mirror/ftp.ncbi.nih.gov") {
## ##
## ##     if( is.character(taxid) )
## ##         taxid <- get.taxonID(taxid)
## ##
## ##     if( is.null(g2acc) )
## ##         g2acc <- get.gene2accession(taxid=taxid, mirror=mirror)
## ##
## ##     if( is.null(ginfo) )
## ##         ginfo <- get.gene_info(taxid=taxid, mirror=mirror)
## ##
## ##     #
## ##     # get rid of the rows where there is no genbankID
## ##     #
## ##     g2acc <- g2acc[!is.na(g2acc$genbankID),]
## ##
## ##     g2acc$genbankVer <- sub("[^.]+", "", g2acc$genbankID)
## ##     g2acc$genbankID <- sub("\\.[0-9]+$", "", g2acc$genbankID)
## ##
## ##     # venn diagram style results
## ##     AnotB <- length(setdiff(genbankIDs, g2acc$genbankID))
## ##     AandB <- length(intersect(genbankIDs, g2acc$genbankID))
## ##     BnotA <- length(setdiff(g2acc$genbankID, genbankIDs))
## ##
## ##
## ## }
##
## get.gene2accession <- function(taxid=10090, mirror="/import/sparta/data/mark/mirror/ftp.ncbi.nih.gov") {
##
##     if( is.character(taxid) )
##         taxid <- get.taxonID(taxid)
##
##     gene2accession <- file.path(mirror, "gene/DATA", "gene2accession.gz")
##     if( !file.exists(gene2accession) )
##         stop(paste(gene2accession, "does not exist\n"))
##
##     #
##     # unzip the gene2accession file
##     #
##     tmp <- paste0(tempfile(pattern="g2acc"), ".gz")
##
##     file.copy(gene2accession, tmp)
##     file.gunzip(tmp)
##     tmp <- sub(".gz", "", tmp)
##
##     #
##     # extract out the specified organism
##     #
##     g2acc.file <- tempfile(pattern="g2acc")
##     system(paste0("cd ", tempdir(), " && ", "grep -w '^", taxid, "' ", tmp, " > ", g2acc.file ))
##     file.remove( tmp )
##
##     g2acc <- read.delim(g2acc.file, as.is=T, header=F, na.strings="-")
##     colnames( g2acc ) <- c("tax_id", "GeneID", "status",
##                            "genbankID", "genbankGI",
##                            "proteinID", "proteinGI",
##                            "genomeAcc", "genomeGI", "genomeStart", "genomeEnd", "genomeOri")
##
##     file.remove( g2acc.file )
##
##     return( g2acc )
## }
##
##
##
## get.gene_info <- function(taxid=10090, mirror="/import/sparta/data/mark/mirror/ftp.ncbi.nih.gov") {
##
##     if( is.character(taxid) )
##         taxid <- get.taxonID(taxid)
##
##     gene_info <- file.path(mirror, "gene/DATA", "gene_info.gz")
##     if( !file.exists(gene_info) )
##         stop(paste(gene_info, "does not exist\n"))
##
##     #
##     # unzip the gene_info file
##     #
##     tmp <- paste0(tempfile(pattern="ginfo"), ".gz")
##
##     file.copy(gene_info, tmp)
##     file.gunzip(tmp)
##     tmp <- sub(".gz", "", tmp)
##
##     #
##     # extract out the specified organism
##     #
##     ginfo.file <- tempfile(pattern="ginfo")
##     system(paste0("cd ", tempdir(), " && ", "grep -w '^", taxid, "' ", tmp, " > ", ginfo.file ))
##     file.remove( tmp )
##
##     ginfo <- read.delim(ginfo.file, as.is=T, header=F, na.strings="-")
##     colnames( ginfo ) <- c("tax_id", "GeneID", "symbol",
##                            "locusTag", "synonyms",
##                            "dbXrefs", "chr",
##                            "mapLoc", "description", "type",
##                            "officialSymbol", "officialFullName", "officialStatus")
##
##     file.remove( ginfo.file )
##
##     return( ginfo )
## }

#' get.taxonID
#' 
#' return taxon ID for a species name or short code
#'
#' @param taxid species name, some aliases, or 2 letter code
#' @return numeric taxon ID
#' @author Mark Cowley, 2013-05-30
#' @export
#' @examples
#' get.taxonID("mouse")
#' get.taxonID("mm")
#' get.taxonID("hs")
#' get.taxonID("rat")
#' get.taxonID("rn")
get.taxonID <- function(taxid) {
    if( grepl("mouse", taxid, ignore.case=T) || grepl("mm", taxid, ignore.case=T) )
        return( 10090 )
    else if( grepl("human", taxid, ignore.case=T) || grepl("hs", taxid, ignore.case=T) )
        return( 9606 )
    else if( grepl("rat", taxid, ignore.case=T) || grepl("rn", taxid, ignore.case=T) )
        return( 10116 )
    else if( grepl("fly", taxid, ignore.case=T) || grepl("fruit", taxid, ignore.case=T) ||
             grepl("drosophila", taxid, ignore.case=T) || grepl("dm", taxid, ignore.case=T) )
        return( 7227 )
    else if( grepl("cress", taxid, ignore.case=T) || grepl("arabidopsis", taxid, ignore.case=T) )
        return( 3702 )
    else if( grepl("yeast", taxid, ignore.case=T) || grepl("sacchar", taxid, ignore.case=T) )
        return( 4932 )
    else {
        stop(paste0("unsupported organism -- suggest looking up Taxa ID at\n",
                "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Taxonomy\n"))
        return( 0 )
    }
}
