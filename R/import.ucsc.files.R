# import the knownGene file downloaded manually from UCSC genome browser.
# By default it looks in the /biomirror/hgdownload.cse.ucsc.edu/tables for
# a file called knownGene
#
# 	#name          chrom   strand   txStart txEnd   cdsStart cdsEnd <snip>
# 	NM_001011874    chr1    -       3174829 3630690 3176369 3630540 <snip>
# 	AK035064        chr1    +       3709446 3783878 3709490 3783206 <snip>
# 	NM_011283       chr1    -       4304522 4320770 4304977 4313203 <snip>
# 	NM_011441       chr1    -       4451317 4456803 4452105 4453796 <snip>
# 	<snip>
#
# Mark Cowley, 8 August 2006
#
import.knowngene <- function(file="/biomirror/hgdownload.cse.ucsc.edu/tables/knownGene") {
	import.ucsc.file(file)
}


# import the kgXref file downloaded manually from UCSC genome browser.
# By default it looks in the /biomirror/hgdownload.cse.ucsc.edu/tables for
# a file called kgXref.gz
#
# 	#kgID      mRNA      spID      spDisplayID   geneSymbol  <snip>
# 	AB001750   AB001750  O35551    RABE1_MOUSE   AB001750    <snip>
# 	AB001926   AB001926  P48722-2  P48722-2      AB001926    <snip>
# 	AB003147   AB003147  O54770    O54770_MOUSE  AB003147    <snip>
# 	AB003503   AB003503  O88180    O88180_MOUSE  AB003503    <snip>
# 	<snip>
#
# Mark Cowley, 8 August 2006
#
import.kgXref <- function(file="/biomirror/hgdownload.cse.ucsc.edu/tables/kgXref.gz") {
	tmp <- p(tempfile(pattern="ucsc"), ".gz")
	file.copy(file, tmp)
	file.gunzip(tmp)
	tmp <- sub("\\.gz", "", tmp)

	res <- import.ucsc.file(tmp)

	unlink(tmp)

	return( res )
}


# File from UCSC tend to have the column names in the first line
# preceeded with a hash character.
#
# This function imports the data and fixes the names to remove
# this hash char.
#
# Mark Cowley, 8 August 2006
#
import.ucsc.file <- function(file) {
	res <- read.delim(file, as.is=T, comment.char="")
	colnames(res)[1] <- "name"

	return( res )
}



# Download knownGeneToLocusLink from UCSC Table Browser to:
# Import the knownToLocusLink file, downloaded from the UCSC
# Table Browser.
# taxid: the numerical taxonID
# files: a list with your file name in it.
#
# Mark Cowley, 30/3/07
#
import.knownGene2gene <- function(taxid=TAXID.MOUSE,
                                  files=list(mouse="/biomirror/ucsc/tables/mm8/knownToLocusLink.gz",
                                             human="/biomirror/ucsc/tables/hg18/knownToLocusLink.gz",
                                             rat="/biomirror/ucsc/tables/rn4/knownToLocusLink.gz")
                                  ) {

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
            stop('unsupported taxid; try setting taxid and files=list("/path/to/your/taxid/knownToLocusLink.gz")')
    }


    tmp <- tempfile()
    file <- file.gunzip(file, tmp)
    res <- import.ucsc.file(tmp)
    unlink(tmp)
    colnames(res) <- c("knownGeneID", "GeneID")
    res <- res[order(res$knownGeneID),]
    res <- rm.duplicate.rows(res, issorted=T)
    rownames(res) <- 1:nrow(res)

    res
}
