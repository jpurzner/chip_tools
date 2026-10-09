#' Average RPKM across a set of columns
#'
#' Reads per kilobase per million, averaged over the columns given by
#' `avg_ind`. Rows are sorted by name first, so `counts` and `ex_length` must
#' already correspond by row name.
#'
#' @param counts Features x samples count data frame.
#' @param ex_length Vector or single-column matrix of exonic lengths in bases,
#'   in the row order of the *sorted* `counts`.
#' @param avg_ind Column indices to average over.
#'
#' @return A one-column data frame `avg_RPKM`, row names from `counts`.
#'
#' @section Note:
#' `counts` is reordered by row name but `ex_length` is not, so pass an
#' `ex_length` that is already in sorted row-name order or the lengths will be
#' applied to the wrong features.
#'
#' @seealso [counts2tpm()] and [counts_to_tpm()], which are preferable for
#'   between-sample comparison.
#'
#' @examples
#' \dontrun{
#' count_table2rpkm(counts, exon_lengths, avg_ind = 1:3)
#' }
#' @export
# takes count data table and produces a table with average RPKM values
# need exon length array and what cols of count table to avg
# can only do one average at a time

count_table2rpkm <- function(counts, ex_length, avg_ind) {
	
    #library(GenomicFeatures)
    #txdb=makeTranscriptDbFromUCSC(genome='mm9',tablename='ensGene')
    #ex_by_gene=exonsBy(txdb,'gene')
    
counts <- counts[order(row.names(counts)),]
counts <- as.matrix(counts)
geneLengthsInKB <- as.matrix(ex_length) / 1000
millionsMapped <- colSums(counts) / 1000000
millionsMapped <- as.matrix(millionsMapped)
rpm = sweep(counts, MARGIN = 2, millionsMapped, "/")
rpkm = sweep(rpm, MARGIN = 1, geneLengthsInKB, "/")
avg_RPKM <- rowMeans(rpkm[,avg_ind])
avg_RPKM <- as.data.frame(avg_RPKM)
row.names(avg_RPKM) <- row.names(counts)
return(avg_RPKM)
}