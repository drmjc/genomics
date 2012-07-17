#
# miRbase utilities
#


#' Import miRBase records.
#' From the miRBase FTP site ([1]), there is a file called
#' \dQuote{miRNA.xls}, which contains lots of the useful information about each miRNA.
#' 
#' This function imports the \code{miRNA.xls} file, subsets the mir's to those from a
#' particular species, and re-sorts the resulting file.
#' It can also obtain the locations of the miRNA's, either automatically (ie by
#' mirroring the entire \sQuote{sequences/CURRENT} directory from above FTP site), or by
#' specifying the gff file.
#' 
#' @param miRNA.xls.file the path to the miRNA.xls file
#' @param species.code 3 letter code such as \dQuote{hsa}, \dQuote{mmu}, \dQuote{rno}
#' @param gff if \code{NULL} then method will look for a file called \dQuote{genomes/<species code>.gff}
#'  in the same directory as the \code{miRNA.xls.file}; if \code{NULL}, then no
#'   location info will be inferred; otherwise provide the path to the 
#' \dQuote{<species code>.gff} file
#' @param predict.arm logical: If \code{TRUE}, then predict which arm the mature miR exists on
#'   (3'/5'). Otherwise set to \code{FALSE}.
#' @author Mark Cowley, 2009-01-07
#' @references 
#' [1] \url{ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT}
#' @export
import.miRbase <- function(miRNA.xls.file, species.code="hsa", gff=NULL, predict.arm=TRUE) {
	#
	# import the miRNA.xls file, and clean it up.
	#
	mirbaseall <- read.xls(miRNA.xls.file)
	mirbase <- mirbaseall[grep(species.code, mirbaseall$ID),]

	mirbase <- mirbase[mirorder(mirbase[,2]), ]
	mirbase <- mirbase[,c(1:grep("Mature2_Seq", colnames(mirbase)))]

	mirbase[mirbase == ""] <- NA
	
	cat(nrow(mirbase), "\n")
	
	#
	# determine the locations of the mir's?
	#
	mir.loc <- NULL
	if( is.null(gff) ) { # try to find the file automatically
		gff <- file.path(dirname(miRNA.xls.file), "genomes", paste(species.code, "gff", sep="."))
		if( file.exists(gff) )
			cat("Found the miRNA locations file automatically.\n")
	}
	if( !is.na(gff) && file.exists(gff) ) {
		cat("Grabbing miRNA locations.\n")
		mir.loc <- import.miRbase.locations(gff)
		mir.loc <- vecmerge(mirbase$ID, mir.loc, "ID")
		
		mirbase <- merge(mirbase, mir.loc[,c("ACC", "pos", "strand")], by.x="Accession", by.y="ACC", all.x=TRUE, all.y=FALSE, sort=FALSE)
		mirbase <- mirbase[mirorder(mirbase$ID), ]
	}
	else {
		mirbase$pos <- NA
	}

	#
	# Determine which arm (3' or 5') each miR is on.
	#
	mirbase$Mature1_Arm <- NA
	mirbase$Mature1_Pos <- NA
	mirbase$Mature2_Arm <- NA
	mirbase$Mature2_Pos <- NA
	mirbase <- mirbase[, c('Accession', 'ID', 'Status', 'Sequence', 'pos', 'strand', 'Mature1_Acc', 'Mature1_ID', 'Mature1_Seq', 'Mature1_Arm', 'Mature1_Pos', 'Mature2_Acc', 'Mature2_ID', 'Mature2_Seq', 'Mature2_Arm', 'Mature2_Pos')]

	if( predict.arm ) {
		cat("Predicting which arm the miRNA belongs to.\n")
		#
		# align the Mature1 sequences first
		#
		mirbase$Mature1_Arm <- NA
		idx <- which( !is.na(mirbase$Mature1_Seq) )
		alignment <- miRNA.align.mir.to.miR(
			mirbase$Mature1_ID[idx], mirbase$Mature1_Seq[idx], 
			mirbase$ID[idx], mirbase$Sequence[idx])
		mirbase$Mature1_Arm[idx] <- miRNA.which.arm(detailed.alignment=alignment)
		
		# given this detailed alignment, and the mir locations, determine the miR locations.
		if( !is.null(mir.loc) ) {
			idx <- intersect(idx, which( !is.na(mirbase$pos) ))
			mirbase$Mature1_Pos[idx] <- ucsc.pos.subset(mirbase$pos[idx], alignment$Start, alignment$End, mirbase$strand[idx])
		}
		
		#
		# align the Mature2 sequences next
		#
		mirbase$Mature2_Arm <- NA
		idx <- which( !is.na(mirbase$Mature2_Seq) )
		alignment <- miRNA.align.mir.to.miR(
			mirbase$Mature2_ID[idx], mirbase$Mature2_Seq[idx], 
			mirbase$ID[idx], mirbase$Sequence[idx])
		mirbase$Mature2_Arm[idx] <- miRNA.which.arm(detailed.alignment=alignment)
		
		# given this detailed alignment, and the mir locations, determine the miR locations.
		if( !is.null(mir.loc) ) {
			idx <- intersect(idx, which( !is.na(mirbase$pos) ))
			mirbase$Mature2_Pos[idx] <- ucsc.pos.subset(mirbase$pos[idx], alignment$Start, alignment$End, mirbase$strand[idx])
		}
		# mirbase$Mature2_Arm <- ""
		# idx <- which( mirbase$Mature2_Seq != "" )
		# alignment <- miRNA.align.mir.to.miR(
		# 	mirbase$Mature2_ID[idx], mirbase$Mature2_Seq[idx], 
		# 	mirbase$ID[idx], mirbase$Sequence[idx])
		# mirbase$Mature2_Arm[idx] <- miRNA.which.arm(detailed.alignment=alignment)
		# 
		# # given this detailed alignment, and the mir locations, determine the miR locations.
		# if( !is.null(mir.loc) ) {
		# 	mirbase$Mature2_Pos[idx] <- ucsc.pos.subset(mirbase$pos[idx], alignment$Start, alignment$End, mir.loc$strand[idx])
		# }
	}
	else {
		# mirbase$Mature1_Arm <- NULL
		# mirbase$Mature1_Pos <- NULL
		# mirbase$Mature2_Arm <- NULL
		# mirbase$Mature2_Pos <- NULL
	}
	
	mirbase
}



#' Import the genomic locations of the miRBase precursor miRNA's.
#' 
#' @param gff.file the path to the gff file, avalable from references [1-3]
#' @return see \code{\link{import.gff}}, with 3 additional columns: \dQuote{ACC}, 
#' \dQuote{ID} and \dQuote{pos}, eg \dQuote{MI0006363}, \dQuote{hsa-mir-1302-2}, 
#' \dQuote{chr1:20229-20366}
#' @author Mark Cowley, 2009-01-07
#' @references 
#' [1] \url{ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/genomes/hsa.gff}
#' [2] \url{ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/genomes/mmu.gff}
#' [3] \url{ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/genomes/rno.gff}
#' @examples
#' \dontrun{
#' tmp <- import.miRbase.locations("~/data/miRBase/sequences/12.0/genomes/hsa.gff")
#' head(tmp)
#' #         ACC           ID chr     start       end strand score                       pos
#' # 1 MI0000060 hsa-let-7a-1   9  95978060  95978139      +     .    chr9:95978060-95978139
#' # 2 MI0000061 hsa-let-7a-2  11 121522440 121522511      -     . chr11:121522440-121522511
#' # 3 MI0000062 hsa-let-7a-3  22  44887293  44887366      +     .   chr22:44887293-44887366
#' # 4 MI0000063   hsa-let-7b  22  44888230  44888312      +     .   chr22:44888230-44888312
#' # 5 MI0000064   hsa-let-7c  21  16834019  16834102      +     .   chr21:16834019-16834102
#' # 6 MI0000065   hsa-let-7d   9  95980937  95981023      +     .    chr9:95980937-95981023
#' }
#' @export
import.miRbase.locations <- function(gff.file) {
	mir.loc <- import.gff(gff.file)
	mir.loc$ACC <- sub("ACC=", "", sub(";.*", "", mir.loc$group))
	mir.loc$ID <- sub("^.*; ID=", "", sub(";$", "", mir.loc$group))
	mir.loc$pos <- ucsc.pos(mir.loc$chr, mir.loc$start, mir.loc$end)

	# rownames(mir.loc) <- make.unique(paste("chr", mir.loc$chr, sep=""))
	mir.loc <- mir.loc[mirorder(mir.loc$ID), ]
	rownames(mir.loc) <- 1:nrow(mir.loc)
	
	mir.loc <- mir.loc[,c('ACC', 'ID', 'chr', 'start', 'end', 'strand', 'score', 'pos')]
	mir.loc
}


#' Convert mirbase database
#' Convert the miRbase db (which has one row per precursor miRNA) to have one
#' row per mature miRNA.
#' 
#' @param mir2miR see \code{\link{import.miRbase}}
#' @param collapse logical: collapse multiple rows into one? separated by \code{sep}
#' @param sep see \code{collapse}
#' @return a \code{data.frame} with columns: \dQuote{Mature_Acc}, \dQuote{Mature_ID}, \dQuote{Mature_Seq},
#'   \dQuote{Precursor_Acc}, \dQuote{Precursor_ID}, \dQuote{Precursor_Status}, \dQuote{Precursor_Seq},
#'   \dQuote{miRBase.order}. This final column indicates if it was in the Mature1 or
#'   Mature2 columns from the original miRNA file
#' @author Mark Cowley, 2009-01-07
#' @export
miRbase.miR2mir <- function(mir2miR, collapse=TRUE, sep=" // ") {
	generic.columns <- intersect(c('Accession', 'ID', 'Status', 'Sequence', 'pos'), colnames(mir2miR))
	
	specific.columns <- intersect(c('Mature1_Acc', 'Mature1_ID', 'Mature1_Seq', 'Mature1_Arm', 'Mature1_Pos'), colnames(mir2miR))
	tmp1 <- mir2miR[,c(specific.columns, generic.columns)]
	colnames(tmp1) <- sub("1", "", colnames(tmp1))
	tmp1$miRBase.order <- "Mature1"
	
	specific.columns <- intersect(c('Mature2_Acc', 'Mature2_ID', 'Mature2_Seq', 'Mature2_Arm', 'Mature2_Pos'), colnames(mir2miR))
	tmp2 <- mir2miR[,c(specific.columns, generic.columns)]
	miRuniverse <- mirsort(setdiff(union(mir2miR$Mature1_ID, mir2miR$Mature2_ID), c("", NA)))
	colnames(tmp2) <- sub("2", "", colnames(tmp2))
	tmp2$miRBase.order <- "Mature2"
	
	if( any(is.na(tmp1$Mature_ID)) )
		tmp1 <- tmp1[!is.na(tmp1$Mature_ID), ]
	if( any(is.na(tmp2$Mature_ID)) )
		tmp2 <- tmp2[!is.na(tmp2$Mature_ID), ]
	
	miR2mir <- rbind(tmp1, tmp2)
	
	miR2mir <- miR2mir[mirorder(miR2mir$Mature_ID), ] # since miR2mir should already have 1's before the 2's, no need to addtionally sort by this argument - plus mirorder only handle's one argument, unlike the more flexible order(...)
	rownames(miR2mir) <- 1:nrow(miR2mir)

	# fix the precursor column names
	colnames(miR2mir) <- sub("^Accession", "Precursor_Acc", colnames(miR2mir))
	colnames(miR2mir) <- sub("^ID", "Precursor_ID", colnames(miR2mir))
	colnames(miR2mir) <- sub("^Status", "Precursor_Status", colnames(miR2mir))
	colnames(miR2mir) <- sub("^Sequence", "Precursor_Seq", colnames(miR2mir))
	colnames(miR2mir) <- sub("^pos", "Precursor_pos", colnames(miR2mir))
	
	# now fix the mature miRNA column names.
	colnames(miR2mir) <- sub("Mature_", "", colnames(miR2mir))
	
	if( collapse ) {
		miR2mir <- collapse.rows(miR2mir, key.column="ID", grep("Precursor|miRBase.order", colnames(miR2mir)), sep)
		miR2mir <- miR2mir[mirorder(miR2mir$ID), ]
	}
	
	miR2mir
}


#' Predict miRNA arm
#' Predict which arm (3' or 5') of the precursor miRNA that the mature miR is
#' on This uses purely the sequences of the precursor, and the mature miRNA.
#' 
#' @param mature.seq a vector of N mature miRNA sequences
#' @param precursor.seq a vector of N precursor miRNA sequences
#' @param THRESH The alignment score (see notes) that distinguishes the 3'/5'
#'   arms
#' @param detailed.alignment If you've already run \code{\link{miRNA.align.mir.to.miR}}, then
#'   you can provide this argument, instead of \code{mature.seq} and \code{precursor.seq}
#' @return a vector of N values that can be either "5p" or "3p"
#' @note Alignment score is the alignments' start index relative to the
#'   precursor divided by the length of the precursor
#' @author Mark Cowley, 2009-01-08
#' @export
miRNA.which.arm <- function(mature.seq, precursor.seq, THRESH=0.44, detailed.alignment=NULL) {

	if( is.null(detailed.alignment) ) {
		stopifnot( length(mature.seq) == length(precursor.seq) )
		alignments <- miRNA.align.mir.to.miR(rep("miR", length(mature.seq)), mature.seq, rep("mir", length(mature.seq)), precursor.seq)
	}
	else {
		alignments <- detailed.alignment
	}

	res <- ifelse(alignments$Align.Prop < THRESH, "5p", "3p")
	
	# N <- length(mature.seq)
	# res <- rep(NA, N)
	# 
	# for(i in 1:N) {
	# 	tmp <- pairwiseAlignment(pattern=mature.seq[i], 
	# 							 subject=precursor.seq[i], 
	# 							gapOpening=0, gapExtension=-100, type="local")
	# 	start.prop <- start(subject(tmp)) / nchar(precursor.seq[i])
	# 	
	# 	if( start.prop < THRESH )
	# 		res[i] <- "5p"
	# 	else
	# 		res[i] <- "3p"
	# }
	res
}

#' Align miR to mir sequences
#' Sequence based alignment of mature miR sequences within precursor mir
#' sequences
#' 
#' @param mature.id a vector of identifiers for the mature and precursor
#' @param precursor.id a vector of identifiers for the mature and precursor
#'   miR/mir genes respectively.
#' @param mature.seq a vector of sequences. same length as each other, and
#' @param precursor.seq a vector of sequences. same length as each other, and
#'   same length as the names
#' @return a table of N rows, containing: Mature ID, length of sequence,
#'   Precursor ID, length of sequence, Alignment score (see
#'   ?pairwiseAlignment), start and end coordinates from within the precursor
#'   align.prop is the align start / precursor length
#' @author Mark Cowley, 2009-01-13
#' @export
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom IRanges subject score
miRNA.align.mir.to.miR <- function(mature.id, mature.seq, precursor.id, precursor.seq) {
	stopifnot( length(mature.seq) == length(precursor.seq) )

	N <- length(mature.seq)
	res <- as.data.frame(matrix(c("miR-1", 0, "mir-1", 100, 0, 1, 21, 0.1), nrow=N, ncol=8, byrow=TRUE), stringsAsFactors=FALSE)
	colnames(res) <- c("Mature_ID", "Mature_Len", "Precursor_ID", "Precursor_Len", "Align.Score", "Start", "End", "Align.Prop")
	colclasses(res) <- c("character", "numeric", "character", "numeric", rep("numeric", 4))
	
	for(i in 1:N) {
		tmp <- pairwiseAlignment(pattern=mature.seq[i], 
								 subject=precursor.seq[i], 
								 gapOpening=0, gapExtension=-100, type="local")
		res$Mature_ID[i] <- mature.id[i]
		res$Mature_Len[i] <- nchar(mature.seq[i])
		res$Precursor_ID[i] <- precursor.id[i]
		res$Precursor_Len[i] <- nchar(precursor.seq[i])
		res$Align.Score[i] <- score(tmp)
		res$Start[i] <- start(subject(tmp))
		res$End[i] <- end(subject(tmp))
		res$Align.Prop[i] <- res$Start[i] / res$Precursor_Len[i]
	}
	res
}

#' Fix the names of microRNA's
#' 
#' If there are any systematic naming problems with miRs, this method can fix
#' them.
#' Currently implemented rules:\cr
#' replace -star -> *
#' 
#' @param miRids a vector of miR ID's that need checking for validity.
#' @return a vector of miR ID's that have been fixed according to a set of
#'   simple rules.
#' @author Mark Cowley, 2009-01-13
#' @export
miRNA.fix.miR.ids <- function(miRids) {
	miRids <- sub("-star", "*", miRids)
	miRids
}


#' Import the miRNA.dead file
#' Import the miRNA.dead file, and return the mir ID's that have been dropped.
#' 
#' @param file the path to \dQuote{miRNA.dead.gz}. It can be unzipped or gzipped
#' @param species \code{NULL} means import all species, otherwise use a 3 letter
#'   species code, eg \dQuote{hsa}
#' @author Mark Cowley, 2009-01-13
#' @export
import.mirbase.dropped <- function(file, species="hsa") {
	if( file.isgzip(file) )
		file <- file.gunzip(file, tempfile())
	all.lines <- readLines(file)
	id.lines <- all.lines[grep("^ID", all.lines)]
	ids <- sub("ID[ \t]+", "", id.lines)

	if( !is.null(species) )
		ids <- grep(species, ids, value=TRUE)
	
	ids
}

#' Validate miR ID's
#' 
#' EXPERIMENTAL: Take a vector of miR ID's, and the formatted miRbase table, and determine
#' the latest miR ID for each supplied miRNA
#' ALGORITHM
#' for each miR id\cr
#' if miR ID found in miRBase, then no change\cr
#' else if miR ID ends in -5p or -3p, and there's only one precursor mir and
#' that mir has a miR on the correct arm (5' or 3') then the correct miR is
#' returned.
#' 
#' This handles cases where miRBase uses */<> whereas query uses
#' -5p/-3p
#' else if miR ID ends in *, and there's only one precursor mir, and that mir
#' has two miRs, then use the newest miR ID. This usually corresponds to
#' records in miRbase that used to end in */<>, but now end in -5p/-3p
#' else if miR ID is plain (eg miR-131) and there's only one precursor mir,
#' then choose the oldest miR record. If the mir has only one miR then simple.
#' if mir has 2 miRs then choose the oldest one. Historically mir-123 and
#' miR-123 would have been deposited into db when discovered, and later a
#' (usually) secondard (and thus miR-123*) would have been discovered and
#' deposited in db with a higher MIMAT accession.
#' 
#' @author Mark Cowley, 2009-01-13
#' @export
miRNA.validate.miR.ID <- function(miRids, mirbase5p3p, dropped.mirs=NULL) {
	if( is.null(dropped.mirs) )
		dropped.mirs <- NULL

	miRids.original <- miRids
	miRids <- miRNA.fix.miR.ids(miRids.original)
	
	miRuniverse <- miRuniverse(mirbase5p3p)
	
	res <- data.frame(user.supplied.id=miRids.original, cleaned.id=miRids, miRbase.id=NA, stringsAsFactors=FALSE)
	res$explanation <- ""
	res$ambiguities <- NA
	
	# optional speedup step. - map the easy ones first.
	
	#
	# for each miR, determine the miRBase ID for it.
	#
	for(i in 1:length(miRids)) {
		miRid <- miRids[i]
		
		#
		# 1. is the miRid in the miRuniverse once?
		#
		if( miRid %in% miRuniverse ) {
			res$miRbase.id[i] <- miRid
			res$explanation[i] <- "NOCHANGE: miR ID found"
		}
		#
		# 2. does it end in -5p or -3p?
		#
		else if( str.right(miRid, 3) %in% c("-5p", "-3p") ) {
			mirid <- miR2mir(miRid)
			mir.row <- match(mirid, mirbase5p3p$ID)
			if( mirid %in% mirbase5p3p$ID ) {
				tmp <- unlist(mirbase5p3p[mir.row, c("Mature5p_ID", "Mature3p_ID")])
				if( str.right(miRid, 3) == "-5p" )
					tmp <- tmp[1]
				else
					tmp <- tmp[2]
				res$miRbase.id[i] <- tmp
				res$explanation[i] <- "CHANGED: choice of */<no star> in miRBase. Used sequence alignment to determine match to query-miR-3p|5p"
			}
			else {
				res$miRbase.id[i] <- NA
				res$explanation[i] <- "ERROR: miR ends in 5p/3p BUT mir not found"
			}
			
		}
		#
		# 3. does it end in *?
		#
		else if( str.right(miRid, 1) == "*" ) {
			# if mir exists, then miR* are generally newer than the miR - choose the one with the newer MIMAT id
			mirid <- miR2mir(miRid)
			mir.row <- match(mirid, mirbase5p3p$ID)
			if( mirid %in% mirbase5p3p$ID ) {
				tmp <- unlist(mirbase5p3p[mir.row, c("Mature5p_Acc", "Mature3p_Acc")])
				tmp <- sub("MIMAT", "", tmp)
				tmp <- as.numeric(tmp)
				tmp.col <- c("Mature5p_ID", "Mature3p_ID")[which.min(tmp)]
				tmp <- mirbase5p3p[mir.row, tmp.col]
				res$miRbase.id[i] <- tmp
				res$explanation[i] <- "CHANGED: choice of two 5p/3p miR's. Chose the newest MIMAT Acc to match the miR*"
			}
			else {
				res$miRbase.id[i] <- NA
				if( any(grepl(mirid, mirbase5p3p$ID)) )
					res$explanation[i] <- "ERROR: Multiple precursor mir genes found (from query miR*)"
				else
					res$explanation[i] <- "ERROR: precursor not in miRBase (from query miR*)"
			}
		}
		#
		# 4. The miR was not found, and doesn't end in -3p, -5p, or *
		# eg "hsa-miR-125a"
		#
		else {
			# the miR could be ambiguous.
			ambig.pattern <- paste(miRid, "[a-z]", sep="")
			ambig.results <- grep(ambig.pattern, miRuniverse, value=TRUE)
			if( length(ambig.results) > 0 ) {
				res$miRbase.id[i] <- NA
				res$explanation[i] <- "ERROR: miR ID is ambiguous"
				res$ambiguities[i] <- paste(ambig.results, collapse=" // ")
				next
			}
			# the miR could have been the major miR for a mir that now has 
			# two miRs, labelled with -5p or -3p. 
			# Most of the time, the oldest miR of the 2 would be the major one.
			# Choose the oldest MIMAT id.
			mirid <- miR2mir(miRid)
			mir.row <- match(mirid, mirbase5p3p$ID)
			if( mirid %in% mirbase5p3p$ID ) {
				tmp <- unlist(mirbase5p3p[mir.row, c("Mature5p_Acc", "Mature3p_Acc")])
				tmp <- sub("MIMAT", "", tmp)
				tmp <- as.numeric(tmp)
				tmp.col <- c("Mature5p_ID", "Mature3p_ID")[which.min(tmp)]
				tmp <- mirbase5p3p[mir.row, tmp.col]
				res$miRbase.id[i] <- tmp
				res$explanation[i] <- "CHANGED: choice of two 5p/3p miR's. Chose the oldest MIMAT Acc to match the miR$"
			}
			else if( length(grep(p(mirid, "-[1-9]"), mirbase5p3p$ID)) > 1 ) {
				res$miRbase.id[i] <- NA
				res$explanation[i] <- "ERROR: Multiple precursor mir genes found (from query miR)"
			}
			else if( mirid %in% dropped.mirs ) {
				res$miRbase.id[i] <- NA
				res$explanation[i] <- "ERROR: precursor mir dropped from miRBase"				
			}
			else {
				res$miRbase.id[i] <- NA
				res$explanation[i] <- "ERROR: precursor mir not in miRbase (from query miR)"
			}
		}
	}
	#
	# DEBUG:
	# Use grep to find miR matches in the miRuniverse
	#
	grep.patterns <- sub("\\*", "\\\\*", miRids)
	grep.patterns <- sub("([0-9a-z])$", "\\1[^0-9]", miRids)
	grep.results <- mgrep(grep.patterns, miRuniverse, value=TRUE, perl=TRUE)
	res$grep.patterns <- grep.patterns
	res$grepID <- sapply(grep.results, paste, collapse=" // ")
	
	res
}



#' Convert miR ID's to 3p/5p
#' The miRNA.xls file from miRbase uses the notation Mature1 and Mature2.
#' 
#' If a mir has 2 miRs then the Mature 1 and 2 are always the 5' and 3'
#' respectively. BUT if there's only one known miR for a mir, then you don't
#' know if it's a 5' or 3'.
#' This function converts the mirbase file from Mature1/2 to Mature5p/3p
#' @param mirbase see \code{\link{import.miRbase}}. You must specify \code{predict.arm=TRUE}
#' 
#' @return a data.frame same dim as mirbase, but with mature1 and 2 changed for
#'   5p and 3p, and miR data shuttled into the correct columns depending on
#'   whether it's on the 5' or 3' arm. The Mature1_Arm/Mature2_Arm columns are
#'   dropped though.
#' @author Mark Cowley, 2009-01-09
#' @export
mirbase.to.3p.5p <- function(mirbase) {
	res <- mirbase
	colnames(res) <- sub("Mature1", "Mature5p", colnames(res))
	colnames(res) <- sub("Mature2", "Mature3p", colnames(res))
	res[,grep("Mature", colnames(res))] <- NA
	
	idx <- which(mirbase$Mature1_Arm == "5p")
	if( length(idx) > 0 )
		res[idx, grep("Mature5p", colnames(res))] <- mirbase[idx, grep("Mature1", colnames(mirbase))]
	# Mature2 shouldn't contain any 5p miRs
	idx <- which(mirbase$Mature2_Arm == "5p")
	if( length(idx) > 0 )
		res[idx, grep("Mature5p", colnames(res))] <- mirbase[idx, grep("Mature2", colnames(mirbase))]

	idx <- which(mirbase$Mature1_Arm == "3p")
	if( length(idx) > 0 )
		res[idx, grep("Mature3p", colnames(res))] <- mirbase[idx, grep("Mature1", colnames(mirbase))]
	idx <- which(mirbase$Mature2_Arm == "3p")
	if( length(idx) > 0 )
		res[idx, grep("Mature3p", colnames(res))] <- mirbase[idx, grep("Mature2", colnames(mirbase))]

	res$Mature5p_Arm <- NULL
	res$Mature3p_Arm <- NULL

	res
}


#' Obtain the miR universe
#' What are the known miR's in a given mirbase table.
#' 
#' the mirbase table can be native (contains Mature1 or Mature2 columns), or
#' updated to Mature5p/Mature3p style, otherwise it can be a miRbase where the
#' miR's
#' are in the ID column
#' @param mirbase a mirbase table
#' @return a vector of all known miR's
#' @author Mark Cowley, 2009-01-09
#' @export
miRuniverse <- function(mirbase) {
	if( any(grepl("Mature1", colnames(mirbase))) )
		res <- setdiff(union(mirbase$Mature1_ID, mirbase$Mature2_ID), c("", NA))
	else if( any(grepl("Mature5p", colnames(mirbase))) )
		res <- setdiff(union(mirbase$Mature5p_ID, mirbase$Mature3p_ID), c("", NA))
	else if( any(grepl("Precursor_ID", colnames(mirbase))) ) {
		# then this is a miRbase
		res <- mirbase$ID
	}
	else
		stop("Unsupported mirbase style.\n")

	res <- mirsort(res)
	res
}





########################################################
# miRNA utility functions.                            ##
# - conversion from precursor to mature, sorting...   ##
########################################################


#' Convert pre to mature microRNA's
#' Convert a list of precursor mir ID's to all possible mature
#' sequences/synonyms
#' 
#' @param mirs a character vector of precursor mirbase ID's
#' @param patterns the possible patterns that make up a mature miRNA ID
#' @param major.only return only the major form of the mature miRNA? (ie only 1
#'   value)
#' @return a list with N elements, one for each mir in 'mirs'. Each element is
#'   a character vector of possible miRNA ID's.
#' @author Mark Cowley, isoD
#' @export
mirna.pre2mature <- function(mirs, patterns=c("*", "-5p", "-3p"), major.only=FALSE) {
	# strip the -1, -2, -3 from the end of the precursor
	miRs <- sub("-[0-9]$", "", mirs)
	miRs <- sub("mir", "miR", miRs)
	# work out the possible synonyms for the mature mir's
	res <- list()
	res[length(miRs)] <- NA # sets the whole list to NA's

	if( ! "" %in% patterns)
		patterns <- c("", patterns)

	for(i in 1:length(miRs)) {
		res[[i]] <- paste(miRs[i], patterns, sep="")
	}
	
	if( major.only )
		res <- sapply(res, "[", which(patterns==""))

	names(res) <- mirs
	
	return( res )
}

#' Convert pre to mature microRNA's
#' Wrapper around \code{\link{mirna.pre2mature}}
#' 
#' @param mirs a character vector of precursor mirbase ID's
#' @param patterns the possible patterns that make up a mature miRNA ID
#' @param major.only return only the major form of the mature miRNA? (ie only 1
#'   value)
#' @return a list with N elements, one for each mir in 'mirs'. Each element is
#'   a character vector of possible miRNA ID's.
#' @author Mark Cowley, 2008-06-24
#' @export
mir2miR <- function(mirs, patterns=c("*", "-5p", "-3p"), major.only=FALSE) {
	mirna.pre2mature(mirs, patterns, major.only)
}



#' Convert mature miR ID's into the immature, precursor form.
#' 
#' NB, this has no way of knowing which mir gene if there are more than one...
#' @param miRs a character vector of mature miR ID's
#' @param patterns the possible patterns on the end of the miR ID
#' @return a character vector of the mir ID corresponding to each miR ID, with
#'   the caveat that there may be multiple mir genes corresponding to each
#'   mature miR     
#' @author Mark Cowley, 2008-06-25
#' @export
miR2mir <- function(miRs, patterns=c("\\*", "-5p", "-3p")) {
	mirs <- sub("miR", "mir", miRs)
	for( pattern in patterns ) {
		idx <- grep(pattern, mirs)
		if( length(idx) > 0 ) {
			mirs[idx] <- sub(paste(pattern, "$", sep=""), "", mirs[idx] )
		}
	}
	mirs
}





#' Split the parts of a mirana miRNA ID into the consituent parts
#'
#' @param mirs a character vector of mirbase ID's (either precursor or mature), duplicates allowed, NA's will NOT be removed.
#' @return a 4 column \code{data.frame} with colnames: \dQuote{species}, \dQuote{accession}, \dQuote{ids}, \dQuote{isoform} 
#' @export
#' @author Mark Cowley, 2008-06-24
#' @examples
#' mirIDsplit("hsa-let-7a-3")
#' # species   accession   ids  isoform 
#' #   "hsa"       "let"  "7a"      "3"
#'
mirsplit <- function(mirs) {
	mirIDsplit(mirs)
}



#' Split the parts of a miRanda miRNA ID into the constituent parts
#' 
#' @param mirs a character vector of mirbase ID's (either precursor or mature),
#'   duplicates allowed, NA's will NOT be removed.
#' @return a 4 column data.frame with colnames: species, accession, ids,
#'   isoform
#' @author Mark Cowley, 2008-06-24
#' @examples
#' mirIDsplit("hsa-let-7a-3")
#' # [1] "hsa" "let" "7a" "3"
#' @export
mirIDsplit <- function(mirs) {
	.mirsplit <- function(x) {
		havedash <- grep("-", x)
		if( length(havedash) == 0 ) {
			LHS <- x
			RHS <- rep("", length(x))
		}
		else {
			LHS <- RHS <- x
			LHS[havedash] <- sub("-.*", "", x[havedash])
			RHS[havedash] <- sub("^[^-]+-", "", x[havedash])
			if( length(havedash) < length(x) ) {
				nodash <- setdiff(1:length(x), havedash)
				LHS[nodash] <- x[nodash]
				RHS[nodash] <- ""
			}
		}
		# if( grepT("-", x) ) {
		# 	LHS <- sub("-.*", "", x)
		# 	RHS <- sub("^[^-]+-", "", x)
		# }
		# else {
		# 	LHS <- x
		# 	RHS <- ""
		# }
		
		return( list(LHS, RHS) )
	}
	x <- .mirsplit(mirs) # -> c("hsa", "let-7a-3")
	species <- x[[1]] # hsa
	x <- .mirsplit(x[[2]]) # -> c("let", "7a-3")
	accession <- x[[1]]
	x <- .mirsplit(x[[2]]) # -> c("7a", "3")
	ids <- x[[1]]
	isoform <- x[[2]]
	
	res <- data.frame(species, accession, ids, isoform, stringsAsFactors=FALSE)
	rownames(res) <- make.unique(mirs)
	
	if( length(mirs) == 1 )
	 	res <- unlist(res[1,,TRUE])
		
	return( res )
}



#' Sort a vector of mirbase ID's
#' 
#' @param mirs a character vector of mirbase ID's (either precursor or mature),
#'   duplicates allowed, NA's will be removed.
#' @return the mirs that have been sorted alpha-numerically, such that
#'   hsa-mir-2 comes before hsa-mir-199
#' @author Mark Cowley, 2008-06-24
#' @export
#' @examples
#' # not run
#' # mirsort(affymir)
mirsort <- function(mirs) {
	o <- mirorder(mirs)
	res <- mirs[o]
	return( res )
	# mirs <- na.rm( mirs )
	# mirsplit <- mirIDsplit( mirs )
	# 
	# # then by 1st, 2nd, 3th fields.
	# numerics <- as.numeric(sub("[a-zA-Z*]+$", "", mirsplit[,3]))
	# alphatail <- sub("^[0-9]+", "", mirsplit[,3])
	# o <- order(mirsplit[,1], mirsplit[,2], numerics, alphatail, mirsplit[,4])
	# 
	# res <- mirs[o]
	# return( res )
}


#' Ordering permutation of mirbase ID's
#' Get the index that each ID within a vector of mirbase ID's would be at if \code{mirs}
#' was sorted/
#' 
#' @param mirs a character vector of mirbase ID's (either precursor or mature),
#'   duplicates allowed, NA's will be removed.
#' @return the numeric index that each mir \emph{would} be in if mirs was sorted.
#' @author Mark Cowley, 2008-06-24
#' @export
#' @examples
#' # not run
#' # mirsort(affymir)
mirorder <- function(mirs) {
	mirs <- na.rm( mirs )
	mirsplit <- mirIDsplit( mirs )
	
	# then by 1st, 2nd, 3th fields.
	numerics <- as.numeric(sub("[a-zA-Z*]+$", "", mirsplit[,3]))
	alphatail <- sub("^[0-9]+", "", mirsplit[,3])
	o <- order(mirsplit[,1], mirsplit[,2], numerics, alphatail, mirsplit[,4])
	
	return( o )
}
