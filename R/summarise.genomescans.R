## Function to threshold a pvalue table at various levels, and report how many genes
## pass those linkage thresholds.
## You can look at various subsets of genes from within the pval table, by supplying a
## list argument to 'genes' where each element of the list is a vector of row indices.
##
## Parameters:
##     pvals: a GxM matrix, where G is no. of genes, M is no of markers
##     thresholds: the various levels at which to threshold the data
##     genes: if NULL, then all genes are used;
##            alternatively, it can be a list of numeric vectors which will be used
##            to index the rows in the pvals table
##     file: the name of a text file for the results to be written to
##
## Value:
##     if genes == NULL then a named vector of counts of genes that pass the thresholds
##     is produced. If genes != NULL, then a table of counts is produced, with each set
##     of genes in a different appropriately named column.
##
## Mark Cowley, 22 Feb 2006
##
summarise.genomescan <- function( pvals, thresholds=c(0.05, 0.001, 1e-04, 1e-05, 1e-06, 1e-07),
                                  genes=NULL, file=NULL, which.genes=F, which.markers=F ) {

##     if( is.null(genes) && nrow(pvals) == 23040 ) {
##         de.ids <- scan.text("~/data/bxd3/de755.ids")
##         exp.ids <- scan.text("~/data/bxd3/exp6k.ids")
##
##         genes <- list(match(de.ids, rownames(pvals)), match(exp.ids, rownames(pvals)), 1:23040)
##     }

    if( !which.genes && which.markers )
        stop( "If you specify which.markers, you must also set which.genes=T\n" )

    res <- NULL

    if( !is.null(genes) && is.list(genes) ) {
        res <- matrix(0, length(thresholds), length(genes))

        for(g in 1:length(genes)) {
            tmp <- summarise.genomescan(pvals[genes[[g]], ], thresholds=thresholds, genes=NULL, which.genes=F)
            res[,g] <- tmp
        }
        colnames(res) <- paste(sapply(genes, length), "genes")
        rownames(res) <- names(tmp)
    }
    else {
        minp <- apply(pvals, 1, min, na.rm=T)

        if( !which.genes ) {
            res <- rep(0, length(thresholds))
            for(thresh in 1:length(thresholds)) {
                res[thresh] <- sum(minp < thresholds[thresh])
            }
        }
        else {
            res <- list()
            for(thresh in 1:length(thresholds)) {
                res[[thresh]] <- rownames(pvals)[which(minp < thresholds[thresh])]
            }
        }
        names(res) <- paste("P <", as.character( thresholds ))

        gc( verbose=F )
    }

    if( !is.null(file) && !which.genes )
        write.delim(res, file, row.names=T, row.names="threshold")


    return( res )
}

## summarise.genomescan(brain.M.spt$p, genes=genes, file="brain.M.spt.summary.txt")
