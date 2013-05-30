TAXID.RAT <- 10116
TAXID.MOUSE <- 10090
TAXID.HUMAN <- 9606

#' import.ncbi.file
#' 
#' General wrapper to import and NCBI file that has usually been mirrored to the
#' local disc. It will copy the compressed or uncompressed file, unzip it in the
#' tmp dir, grep out the taxon of interest, resort by a given column, and return
#' the data
#' 
#' @param file path to ncbi-formatted file
#' @param taxid numeric taxon id
#' @param colnames character vector
#' @param sortby.col undocumented
#' 
#' @return \code{data.frame}
#' 
#' @author Mark Cowley, 10 April 2006
#' @export
import.ncbi.file <- function(file, taxid, colnames, sortby.col) {

    file.type <- gsub("\\.gz", "", basename(file))

    stopifnot( file.exists(file) )

    tmp <- tempfile(pattern="ncbi") # this will contain the full-sized uncompressed NCBI file

    if( grepl("\\.gz$", file) ) {
        tmp <- paste0(tmp, ".gz")
        file.copy(file, tmp)
        file.gunzip(tmp)

        tmp <- sub(".gz$", "", tmp)
    }
    else {
        file.copy(file, tmp)
    }

    tmp2 <- tempfile(pattern="ncbi")  # this will contain NCBI file for organism of interest.

    cat( paste0("subsetting to taxon ID: ", taxid, "\n") )
    system(paste("grep -w '^", taxid, "' ", tmp, " > ", tmp2, sep=""))

    #
    # import the subsetted NCBI file
    #
    cat( paste0("importing ", file.type, " for ", taxid, "\n") )
    ncbi <- read.delim( tmp2, as.is=T, header=F )
    colnames(ncbi) <- colnames[1:ncol(ncbi)]

    #
    # cleanup /tmp
    #
    unlink(c(tmp, tmp2))

    if( !is.null(sortby.col) ) {
        if( is.character(sortby.col) )
            idx <- match(sortby.col, colnames)
        else
            idx <- sortby.col

        ncbi <- ncbi[order(ncbi[,idx]), ]
    }

    return( ncbi )
}


#' import a local copy of the gene_info file.
#' 
#' @author Mark Cowley, 10 April 2006
#' @export
#' @seealso \code{\link{import.ncbi.file}}
import.gene_info <- function(file="/biomirror/ftp.ncbi.nih.gov/gene/DATA/gene_info.gz", taxid=10090, add.cM.column=T) {
    colnames <- c("TaxonID", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs", "Chr", "MapLoc", "Description", "Type", "SymbolNomenc", "FullNameNomenc", "StatusNomenc", "Other")
    sortby.col <- "GeneID"

    gi <- import.ncbi.file(file, taxid, colnames, sortby.col)

    if( add.cM.column ) {
        gi$cM <- sub("^.*\\|[^ ]+ ", "", gi$MapLoc)
        gi$cM[-grep("cM", gi$cM)] <- NA
        # this WILL spit out errors.
        suppressWarnings( gi$cM <- as.numeric(sub(" cM", "", gi$cM)) )
    }

    return( gi )
}
# CHANGELOG
# 8/8/2006 robustified the colnaming. new 13th column indicates
# the assembly used. added extra "unknown" column names and make sure
# that import.ncbi.file only tries to give as many column names to the
# table as there are in that table.


#' import a local copy of the gene2refseq file.
#' 
#' @author Mark Cowley, 10 April 2006
#' @export
#' @seealso \code{\link{import.ncbi.file}}
#' @importFrom mjcbase colclasses "colclasses<-"
import.gene2refseq <- function(file="/biomirror/ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz", taxid=10090,
                               make.unique=T) {
    colnames <- c("TaxonID", "GeneID", "Status",
                  "RNAacc", "RNAgi", "PROTacc", "PROTgi", "GENOMEacc", "GENOMEgi",
                  "Start", "Stop", "Strand", "Assembly",
                  "Unknown1", "Unknown2", "Unknown3")

    sortby.col <- "GeneID"

    g2r <- import.ncbi.file(file, taxid, colnames, sortby.col)

    g2r$Status <- toupper(g2r$Status) # there are some entries for which the Reviewed status is in title-case

    if( make.unique )
        g2r <- make.unique.gene2refseq(g2r, taxid=taxid)

    #
    # convert the numeric columns into numerics.
    # this WILL spit out errors when replacing "-" to a numeric.
    #
    suppressWarnings(colclasses(g2r)[match(c("RNAgi", "PROTgi", "GENOMEgi", "Start", "Stop"), colnames(g2r))] <- "numeric")

    return( g2r )
}

#' import.gene2accession
#'
#' @param file path to gene2accession.gz
#' @param taxid numeric taxon id
#' @param make.unique logical
#' 
#' @return undocumented
#' 
#' @author Mark Cowley, 2013-05-30
#' @export
import.gene2accession <- function(file="/biomirror/ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz", taxid=10090,
                               make.unique=F) {
    g2acc <- import.gene2refseq(file, taxid, make.unique=make.unique)

    return( g2acc )
}

#' Import the gene2unigene file.
#' 
#' @param taxid numeric taxon id
#' @param file path to gene2unigene file
#' 
#' @return \code{data.frame}
#' 
#' @author Mark Cowley, 30/3/07
#' @export
import.unigene2gene <- function(taxid=TAXID.MOUSE,
                                file="/biomirror/ftp.ncbi.nih.gov/gene/DATA/gene2unigene") {
    if(taxid==TAXID.MOUSE)
        prefix <- "Mm"
    else if(taxid==TAXID.HUMAN)
        prefix <- "Hs"
    else if(taxid==TAXID.RAT)
        prefix <- "Rn"

    tmp <- tempfile()
    cmd <- paste("grep",sQuote(paste0(prefix, "\\.")),file.path(file),">",tmp)
    system( cmd )

    res <- read.delim(tmp, as.is=T, header=F)[, c(2,1)]
    unlink(tmp)
    colnames(res) <- c("UniGeneID", "GeneID")
    res <- res[order(as.numeric(sub("^.*\\.","",res$UniGeneID))),]
    # res <- rm.duplicate.rows(res, issorted=T)
    res <- res[!duplicated(res), ]
    rownames(res) <- 1:nrow(res)

    res
}

#' Import the Mm.data.gz file file.
#' 
#' @param taxid numeric taxon id
#' @param files list of named paths to Hs.data.gz, Mm.data.gz, Rn.data.gz
#' @param genbank.ids undocumented
#' 
#' return undocumented
#' 
#' @author Mark Cowley, 30/3/07
#' @export
import.unigene.data <- function( taxid=TAXID.MOUSE,
                                 files=list(mouse="/biomirror/bio-mirror.grangenet.net/biomirror/unigene/Mus_musculus/Mm.data.gz",
                                            human="/biomirror/bio-mirror.grangenet.net/biomirror/unigene/Homo_sapiens/Hs.data.gz",
                                            rat="/biomirror/bio-mirror.grangenet.net/biomirror/unigene/Rattus_norvegicus/Rn.data.gz"),
#                                  file="/biomirror/bio-mirror.grangenet.net/biomirror/unigene/Mus_musculus/Mm.data.gz",
                                 genbank.ids=NULL ) {
#                                  fields=c("ID", "TITLE", "GENE", "CYTOBAND", "GENE_ID", "LOCUSLINK", "EXPRESS",
#                                           "RESTR_EXPR", "CHROMOSOME", "STS", "PROTSIM", "SCOUNT", "SEQUENCE" )

    if( taxid == TAXID.MOUSE )
        file <- files$"mouse"
    else if( taxid == TAXID.HUMAN )
        file <- files$"human"
    else if( taxid == TAXID.RAT )
        file <- files$"rat"
    else {
        if( (is.list(files) && length(files) == 1) && file.exists(files[1]) )
            file <- files[[1]]
        else if( is.character(files) && length(files) == 1 && file.exists(files[1]) )
            file <- files[1]
        else
            stop('unsupported taxid; try setting taxid and files=list("/path/to/your/taxid/XX.data.gz")')
    }

    file.name <- gsub("\\.gz", "", basename(file))

    stopifnot( file.exists(file) )

    tmp <- tempfile(pattern="ncbi") # this will contain the full-sized uncompressed NCBI file

    if( grepl("\\.gz$", file) ) {
        cat("unzipping")
        tmp <- paste0(tmp, ".gz")
        file.copy(file, tmp)
        file.gunzip(tmp)

        tmp <- sub(".gz$", "", tmp)
        cat(".\n")
    }
    else {
        file.copy(file, tmp)
    }

    tmp2 <- tempfile(pattern="ncbi")  # this will contain NCBI file for organism of interest.

    cat("parsing")
    system(paste("parse.unigene.data.pl", tmp, tmp2, sep=" "))
    unlink( tmp )
    cat(".\n")

    if( !is.null(genbank.ids) ) {
        #
        # Since the parsed Mm.data becomes very large, it is unwise to import all
        # of this data, especially as much of it will not relate to a particular
        # set of genbank ID's of interest.
        # Thus, use grep -f to subset the parsed unigene data to just those that are
        # in our genbank IDs vector,
        # tmp2 contains the parsed unigene data
        # tmp3 contains the genbank ids
        # tmp4 will contain the lines from tmp2 that are also in tmp3
        # move tmp4 back to tmp2 so that next command can just import tmp2.
        #
        cat("subsetting")
        genbank.ids <- sort( genbank.ids )

        tmp3 <- tempfile(pattern="ncbi")
        write(genbank.ids, tmp3)

        tmp4 <- tempfile(pattern="ncbi")
        system( paste( "grep -f", tmp3, tmp2, ">", tmp4, sep=" ")  )
        unlink( tmp3 ) # rm the genbank ID's

        unlink( tmp2 )
        file.rename(tmp4, tmp2)
        cat(".\n")
    }
    #
    # import the subsetted NCBI file
    #
    ncbi <- read.delim( tmp2, as.is=T, header=T )

    #
    # cleanup /tmp
    #
    unlink( tmp2 )

    return( ncbi )
}


# import.gene2refseq <- function(file="/biomirror/ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz", taxid=10090) {
#
#     stopifnot( file.exists(file) )
#
#     tmp <- file.path(tempdir(), "gene2refseq") # this will contain the full-sized uncompressed gene2refseq
#
#     if( grepl(".gz$", file) ) {
#         tmp <- paste0(tmp, ".gz")
#         file.copy(file, tmp)
#         file.gunzip(tmp)
#
#         tmp <- sub(".gz$", "", tmp)
#     }
#     else {
#         file.copy(file, tmp)
#     }
#
#     tmp2 <- tempfile() # this will contain gene2refseq for organism of interest.
#
#     cat( paste0("subsetting to taxon ID: ", taxid, "\n") )
#     system(paste("grep -w '^", taxid, "' ", tmp, " > ", tmp2, sep=""))
#
#     #
#     # import the subsetted gene2refseq (approx 15 Mb for mouse)
#     #
#     cat( paste0("importing gene2refseq for ", taxid, "\n") )
#     gr <- read.delim( tmp2, as.is=T, header=F )
#     colnames(gr) <- c("TaxonID", "GeneID", "Status", "RNAacc", "RNAgi", "PROTacc", "PROTgi", "GENOMEacc", "GENOMEgi", "Start", "Stop", "Strand")
#
#     #
#     # cleanup /tmp
#     #
#     unlink(c(tmp, tmp2))
#
#     gr <- gr[order(gr$GeneID), ]
#
#     return( gr )
# }

# import.gene_info <- function(file="/biomirror/ftp.ncbi.nih.gov/gene/DATA/gene_info.gz", taxid=10090) {
#
#     stopifnot( file.exists(file) )
#
#     tmp <- file.path(tempdir(), "gene_info") # this will contain the full-sized uncompressed gene_info
#
#     if( grepl(".gz$", file) ) {
#         tmp <- paste0(tmp, ".gz")
#         file.copy(file, tmp)
#         file.gunzip(tmp)
#
#         tmp <- sub(".gz$", "", tmp)
#     }
#     else {
#         file.copy(file, tmp)
#     }
#
#     tmp2 <- tempfile() # this will contain gene_info for organism of interest.
#
#     cat( paste0("subsetting to taxon ID: ", taxid, "\n") )
#     system(paste("grep -w '^", taxid, "' ", tmp, " > ", tmp2, sep=""))
#
#     #
#     # import the subsetted gene_info (approx 15 Mb for mouse)
#     #
#     cat( paste0("importing gene_info for ", taxid, "\n") )
#     ginfo <- read.delim( tmp2, as.is=T, header=F )
#     colnames(ginfo) <- c("TaxonID", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs", "Chr", "MapLoc", "Description", "Type", "SymbolNomenc", "FullNameNomenc", "StatusNomenc")
#
#     #
#     # cleanup /tmp
#     #
#     unlink(c(tmp, tmp2))
#
#     ginfo <- ginfo[order(ginfo$GeneID), ]
#
#     return( ginfo )
# }
