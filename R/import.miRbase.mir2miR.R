#' Import the master csv file from miRbase.
#' 
#' Download from the Sanger FTP site, and this file should be in a directory
#' tree like:
#' miRBase > genoloc | GFF | targets | sequences > 11.0 > miRNA.csv
#' This function returns a 3 column table of mir ID, then miR ID 1, miR ID 2
#' 
#' @param f the lastest miRBase miRNA.csv file (from Sanger FTP)
#' @param species which 3 letter species code? (default hsa)
#' @param as.list if FALSE, the results are a table, it TRUE, the result is a
#'   list with n elements, 1 per miR
#' @author Mark Cowley, 2008-09-22
#' @export
import.miRbase.mir2miR <- function(f="~/data/miRBase/sequences/latest/miRNA.csv", species="hsa", as.list=FALSE) {
	stopifnot(file.exists(f))
	raw <- read.csv(f)
	idx <- grep(p("^", species), raw$ID)
	mir2miRs <- raw[idx, c("ID", "Mature1_ID", "Mature2_ID")]
	colnames(mir2miRs) <- c("mir", "miR1", "miR2")
	
	if( as.list ) {
		tmp <- list()
		# tmp[nrow(mir2miRs)] <- NA
		for(i in 1:nrow(mir2miRs))
			tmp[[i]] <- setdiff(unlist(mir2miRs[i,2:3]), c("", NA))
		# tmp <- apply(mir2miR[,2:3], 1, as.list)
		names(tmp) <- mir2miRs[,1]
		mir2miRs <- tmp
	}
	
	mir2miRs
}

#' Import the master miRBase miRNA table, but re-order the relationships to be
#' from miR to mir
#' 
#' @param f the lastest miRBase miRNA.csv file (from Sanger FTP)
#' @param species which 3 letter species code? (default hsa)
#' @param as.list if FALSE, the results are a table, it TRUE, the result is a
#'   list with n elements, 1 per miR
#' @author Mark Cowley, 2008-09-22
#' @export
import.miRbase.miR2mir <- function(f="~/data/miRBase/sequences/latest/miRNA.csv", species="hsa", as.list=FALSE) {
	mir2miR <- import.miRbase.mir2miR(f, species, FALSE)
	# alter the direction
	major <- mir2miR[,c(2,1)]; colnames(major) <- c("miR", "mir")
	minor <- mir2miR[,c(3,1)]; colnames(minor) <- c("miR", "mir")
	rm.idx <- which(is.na(minor$miR) | minor$miR == "")
	minor <- minor[-rm.idx, ]
	
	merged <- rbind(major, minor)
	merged <- merged[order(merged$miR, merged$mir), ]
	merged <- merged[mirorder(merged$miR), ]
	
	# if( merge ) {
	# 	for(miR in unique(merged$miR)) {
	# 		idx <- which(merged$miR == miR)
	# 		if( length(idx) > 1 ) {
	# 			merged$mir[idx[1]] <- paste(merged$mir[idx], collapse=" /// ")
	# 		}
	# 	}
	# 	merged <- merged[match(unique(merged$miR), merged$miR),]
	# }
	if( as.list ) {
		res <- list()
		for(miR in unique(merged$miR)) {
			idx <- which(merged$miR == miR)
			res[[miR]] <- merged$mir[idx]
		}
		names(res) <- unique(merged$miR)
		merged <- res
	}
	
	merged
}


#' Obtain the set of all immature 'mirs' from miRBase = the 'mir-ome'
#' 
#' @param f the lastest miRBase miRNA.csv file (from Sanger FTP)
#' @param species which 3 letter species code? (default hsa)
#' @author Mark Cowley, 2008-09-22
#' @export
import.miRbase.mirome <- function(f="~/data/miRBase/sequences/latest/miRNA.csv", species="hsa") {
	mir2miRs <- import.miRbase.mir2miR(f, species=species)
	mirs <- unique(mir2miRs$mir)
	mirs <- mirsort(mirs)
	mirs
}

#' Obtain the set of all mature 'miRs' from miRBase = the 'miR-ome'
#' 
#' @param f the lastest miRBase miRNA.csv file (from Sanger FTP)
#' @param species which 3 letter species code? (default hsa)
#' @author Mark Cowley, 2008-09-22
#' @export
import.miRbase.miRome <- function(f="~/data/miRBase/sequences/latest/miRNA.csv", species="hsa") {
	mir2miRs <- import.miRbase.mir2miR(f, species=species)
	miRs <- c(mir2miRs$miR1, mir2miRs$miR2)
	miRs <- setdiff(miRs, c("", NA))
	miRs <- mirsort(miRs)
	miRs
}
