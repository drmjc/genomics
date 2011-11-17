# Make an Entrez Gene Summary table which includes the genes symbol, description,
# genomic mapping info and synonyms
#
# Data is retreived from "gene2refseq" and "gene_info" available from:
#     "ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz" and
#     "ftp.ncbi.nih.gov/gene/DATA/gene_info.gz"
#
# Algorithm:
#   The best genomic location is determined from gene2refseq file
#       (see import.gene2refseq)
#   The relevant columns from gene_info are kept (see import.gene_info)
#
# Value:
#  a data.frame with the following columns:
#    c("GeneID", "Symbol", "Description", "Synonym", "GenbankID", "mgiID", "Chr", "Start", "Stop", "Strand", "cM")
#
#
#       GeneID Symbol                                         Description
#     1  11287    Pzp                              pregnancy zone protein
#     2  11298  Aanat                  arylalkylamine N-acetyltransferase
#     3  11302   Aatk                apoptosis-associated tyrosine kinase
#     4  11303  Abca1 ATP-binding cassette, sub-family A (ABC1), member 1
#     5  11304  Abca4 ATP-binding cassette, sub-family A (ABC1), member 4
#     6  11305  Abca2 ATP-binding cassette, sub-family A (ABC1), member 2
#                                                    Synonyms GenbankID       mgiID
#     1                      A1m /// A2m /// AI893533 /// MAM NM_007376   MGI:87854
#     2                 MGC151344 /// Nat-2 /// Nat4 /// SNAT NM_009591 MGI:1328365
#     3                                   AATYK /// mKIAA0641 NM_007377 MGI:1197518
#     4                                                  Abc1 NM_013454   MGI:99607
#     5 AW050280 /// Abc10 /// Abcr /// D430003I15Rik /// RmP NM_007378  MGI:109424
#     6        AI413825 /// Abc2 /// D2H0S1474E /// mKIAA1062 NM_007379   MGI:99606
#       Chr    Start     Stop Strand   cM
#     1   6 1.28e+08 1.28e+08      - 62.0
#     2  11 1.16e+08 1.16e+08      + 70.0
#     3  11 1.20e+08 1.20e+08      -   NA
#     4   4 5.31e+07 5.32e+07      - 23.1
#     5   3 1.22e+08 1.22e+08      + 61.8
#     6   2 2.53e+07 2.53e+07      + 12.6
#
# checked with Mouse, Rat and Humans.
#
# Mark Cowley, 19/3/07
#
make.gene.summary <- function(outfile="/biomirror/mouse/Gene/gene.summary.mm.txt",
                              taxid=10090) {

    gene_info <- import.gene_info(taxid=taxid)
    gene_info$Chr[!is.na(gene_info$Chr) & gene_info$Chr == "MT"] <- "M"
    gene_info$Synonyms <- gsub("\\|", " /// ", gene_info$Synonyms)

    gene2refseq <- import.gene2refseq(taxid=taxid, make.unique=T)
    gene2refseq$GenbankID <- tsub(gene2refseq$RNAacc)

    gene <- merge(gene_info[,c("GeneID", "Symbol", "Description", "Synonyms", "dbXrefs", "Chr", "cM")],
                  gene2refseq[,c("GeneID", "GenbankID", "Chr", "Start", "Stop", "Strand")],
                  by="GeneID", all=T, sort=T)


    ## remove the cM info for all those that have different chr info, then remove the Chr.x column.
    gene[!is.na(gene$Chr.x) | !is.na(gene$Chr.y),][which(gene$Chr.x[!is.na(gene$Chr.x) | !is.na(gene$Chr.y)]
                                                        !=
                                                        gene$Chr.y[!is.na(gene$Chr.x) | !is.na(gene$Chr.y)]), "cM"] <- NA

    gene$Chr.x <- NULL

    rownames(gene) <- 1:nrow(gene)

    ## rename some columns
    colnames(gene)[match("Chr.y", colnames(gene))] <- "Chr"

    #
    # reorder the columns.
    #
    # c("GeneID", "Symbol", "Description", "mgiID", "chr", "cM", "GenbankID", "Start", "Stop", "Strand")
    gene <- gene[, c("GeneID", "Symbol", "Description", "Synonyms", "GenbankID", "dbXrefs", "Chr", "Start", "Stop", "Strand", "cM")]

    if( taxid == TAXID.MOUSE )
        colnames(gene)[match("dbXrefs", colnames(gene))] <- "mgiID"
    else if( taxid == TAXID.RAT )
        colnames(gene)[match("dbXrefs", colnames(gene))] <- "rgdID"


    ## write out the results.
    if( !is.null(outfile) )
        write.delim(gene, outfile)

    return( gene )
}

#
# all the other NCBI files have an import function, so why not this one?
#
# Mark Cowley, 19/3/07
#
import.gene.summary <- function(file="/biomirror/mouse/Gene/gene.summary.mm.txt") {
    read.delim(file, as.is=T)
}
